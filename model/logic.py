import math
import random
import numpy as np


import logging


class Subregion():
    """ Subregion class.
    
    Each instance of Subregion contains
    the name, and the left and right 
    boundaries of a subregion. 
    
    Name and boundaries are set during 
    initualization and can be retrieved
    afterwards using helper methods.
    """
    def __init__(self, name, left_boundary, right_boundary, iterations):
        self.name = name
        self.left_boundary = left_boundary
        self.right_boundary = right_boundary
        self.flux_list = np.zeros(iterations, dtype=np.int64)
        
    def leftBoundary(self):
        return self.left_boundary
    
    def rightBoundary(self):
        return self.right_boundary
    
    def getName(self):
        return self.name

    def incrementFlux(self, iteration):
        self.flux_list[iteration] += 1
    
    def getFluxList(self):
        return self.flux_list

def get_event_particles(e_events, subregions, model_particles, level_limit, height_dependant=False):
    """ Find and return list of particles to be entrained

    Keyword arguments:
    e_events -- number of events requested per subregion 
    subregions -- array of Subregion objects
    model_particles -- array of all model particles

    Returns:
    event_particles -- List of particles to be entrained

    """
    if e_events == 0:
        e_events = 1 #???
    
    event_particles = []
    for subregion in subregions:
        # Take only particles in the boundaries of the current subregion
        subregion_particles = model_particles[
                (model_particles[:,0] >= subregion.leftBoundary())
              & (model_particles[:,0] <= subregion.rightBoundary())]
        # Take only particles that are in-stream (not ghost)
        in_stream_particles = subregion_particles[
                                                subregion_particles[:,0] != -1]
        # Take only particles that are 'active' 
        active_particles =  in_stream_particles[
                                                in_stream_particles[:,4] != 0]

        # Do not take any particles that have been selected for entrainment (i.e do not double select)
        # This only happens when particles rest on the boundary. 
        active_event, active_idx, event_idx = np.intersect1d(active_particles[:,3], event_particles, return_indices=True)
        active_particles = np.delete(active_particles, active_idx, axis=0)

        subregion_event_ids = []  
        if height_dependant: # any particle at the level limit must be entrained
            levels = elevation_list(subregion_particles[:,2], desc=False)
            tip_particles = []
            # find the tip particles -- these are the particles being entrained
            if len(levels) == level_limit: 
                tip_particles = active_particles[active_particles[:,2] == levels[level_limit-1]]
            for particle in tip_particles:
                subregion_event_ids.append(particle[3])
                active_particles = active_particles[active_particles[:,2] != particle[2]]
        # If there are not enough particles in the subregion to sample from, alter the sample size
        if e_events > len(active_particles):
            random_sample = random.sample(range(len(active_particles)), 
                                        len(active_particles))
        else: 
            random_sample = random.sample(range(len(active_particles)), 
                                        e_events)
        # TODO: change so that we don't rely on loop index to grab particle
        for index in random_sample:
            subregion_event_ids.append(int(active_particles[index][3])  )
        
        ghost_particles = np.where(model_particles[:,0] == -1)[0]
        for index in ghost_particles:
            model_particles[index][0] = 0 
            subregion_event_ids.append(index)
        
        if e_events != len(subregion_event_ids):
            msg = (
                     f'Requested {e_events} events in {subregion.getName()} ' 
                     f'but {len(subregion_event_ids)} are occuring'
            )
            logging.warning(msg)
        event_particles = event_particles + subregion_event_ids
    event_particles = np.array(event_particles, dtype=np.intp)

    return event_particles

def define_subregions(bed_length, num_subregions, iterations):
    """ Define subregion list for model stream.
    

    Keyword arguments:
    bed_length -- length of the model bed.
    num_subregions -- number of subregions to create
    iterations -- number of iterations for this model simulation

    Returns:
    subregions_arr -- array of initialized subregion objects

    """
    subregion_length = bed_length/num_subregions
    left_boundary = 0.0
    subregions_arr = []
    for region in range(num_subregions):  
        right_boundary = left_boundary + subregion_length   
        subregion = Subregion(f'subregion-{region}', left_boundary, right_boundary, iterations)
        left_boundary = right_boundary
        
        subregions_arr.append(subregion)
    
    return subregions_arr
    
def build_streambed(x_max, set_diam):
    """ Build the bed particle list.
    
    
    Handles calls to add_bed_particle, checks for 
    completness of bed and updates the x-extent
    of stream when the packing exceeds/under packs 
    within 8mm range.
    
    Note: the updates to x-extent are only required 
    when variable particle diameter is being used. 
    
    Return values:
    bed_particles -- array of all bed particles

    """
    max_particles = int(math.ceil( x_max / set_diam ))
    bed_particles = np.zeros([max_particles, 7],dtype=float)
    
    particle_id = -1
    centre = (set_diam/2)  
    state = 0
    age = 0
    loop_age = 0
    elevation = 0
    while not bed_complete(centre, x_max):  
        # index with negative indices... bed particles are built from the final element to the first
        bed_particles[particle_id] = [centre, set_diam, elevation, particle_id, state, age, loop_age]
        centre += set_diam
        particle_id += -1 # Bed particles get negative IDs
    
    return bed_particles

def bed_complete(pack_idx, x_max):
    """Check to see if bed is complete based on model params.""" 
    # similarly, if np.count_nonzero(bed_space) == x_max
    if pack_idx >= x_max:
        return 1
    else: return 0


def determine_num_particles(pack_frac, num_vertices):
    """Return the number of model particles to be used, based on 
    the packing fraction"""
    
    num_particles = num_vertices * pack_frac
    num_particles = int(math.ceil(num_particles))
    
    return num_particles

# Trig from: https://math.stackexchange.com/questions/2293201/
def place_particle(particle, model_particles, bed_particles, h):
    """ Calculate new Y of particle based on location in stream.
    
    
    Provided a particle p's location X in stream, 
    search for 2 supporting particles that p will rest
    on when placed at X.
    
    Calculate the Y position of p based on the elevation 
    and horizontal placement of both supporting particles.
    The calcuated X might be different up to some decimal
    point, so both X and Y are rounded to 2 decimal places.

    
    Keyword arguments:
    particle -- array representing the model particle that is being placed
    model_particles -- array of all model particles
    bed_particles -- array of all bed particles

    Return values:
    rounded_x -- rounded float of particle's new x location
    rounded_y -- rounded float of particle's new y location
    left_support -- id of the left support for the placed particle
    right_support -- id of the right support for the placed particle
    
    """
    left_support, right_support = find_supports(particle, model_particles, bed_particles)
    rounded_x = round(particle[0], 2)
    rounded_y = round(np.add(h, left_support[2]), 2)
    return rounded_x, rounded_y, left_support[3], right_support[3]


def update_particle_states(model_particles, model_supports):
    """ Set/update each model particle's state.
    
    
    If any model particle p has a particle 
    resting on it in the stream then p must 
    be set to inactive indicated by a boolean 0.
    
    If p does not have any particles resting
    on top of it then it is considered active 
    indicated by a boolean 1.
    
    Note that bed particles are always considered
    inactive.
    
    
    Keyword arguments:
    model_particles -- array of all model particles
    model_supports -- array of left/right supporting particles 
                                        for each model particle
    bed_particles -- array of all bed particles

    Return values:
    model_particles -- array of all model particles with updated states
    
    """
    # Start by setting all model particles to active then 
    # only set to inactive if there is a particle sitting on top
    model_particles[:,4] = 1
    in_stream_particles = model_particles[model_particles[:,0] != -1]
    inactive_left = np.intersect1d(in_stream_particles[:,3], model_supports[:,0])
    inactive_right = np.intersect1d(in_stream_particles[:,3], model_supports[:,1])

    if inactive_left.size != 0:
        model_particles[inactive_left.astype(int), 4] = 0
    if inactive_right.size != 0:
        model_particles[inactive_right.astype(int), 4] = 0

    return model_particles


def find_supports(particle, model_particles, bed_particles):
    """ Find the 2 supporting particles for a given particle.
    


    Provided a particle p at location x, this function 
    will search the stream for particles that p would
    rest on, if being dropped at x. 
    
    More generally, supporting particles are those 
    particles which directly hold up a particle. Supporting
    particles will always have a centre location that is 
    exactly a radius length away from p's centre, since 
    every particle has the same diameter.
    
    Keyword arguments:   
    particle -- array representing a particle 
    model_particles -- model particle list
    bed_particles -- bed particle list

    Returns:
    left_support -- the left supporting particle
    right_support -- the right supporting particle
    """  
    all_particles = np.concatenate((model_particles, bed_particles), axis=0)
    # Define location where left and right supporting particles must sit
    # in order to be considered a supporting particle.
    # Note: This limits the model to using same-sized grains.
    left_center = particle[0] - (particle[1] / 2)
    right_center = particle[0] + (particle[1] / 2)
     
    l_candidates = all_particles[all_particles[:,0] == left_center]
    try:
        left_support = l_candidates[l_candidates[:,2] 
                                    == np.max(l_candidates[:,2])]
    except ValueError:
        error_msg = (
                     f'No left supporting particle at {left_center}' 
                     f'for particle {particle[3]}'
        )
        logging.error(error_msg)
        raise   
        
    r_candidates = all_particles[all_particles[:,0] == right_center] 
    try:
        right_support = r_candidates[r_candidates[:,2] == np.max(r_candidates[:,2])]
    except ValueError:
        error_msg = (
                     f'No right supporting particle at {right_center}' 
                     f'for particle {particle[3]}'
        )
        logging.error(error_msg)
        raise
    return left_support[0], right_support[0]


def set_model_particles(bed_particles, available_vertices, set_diam, pack_fraction, h):
    """ Create array of model particles set each particle in stream.
    
    Locations of model particles are randomly
    assigned based on the vertices that are available
    from the bare bed. 

    As a reminder the structure of a model particle is:
        [0] = horizontal location (centre),
        [1] = diameter,
        [2] = elevation (centre),
        [3] = uid,
        [4] = active (boolean)
        [5] = age counter
        [6] = loop age counter
    This should be the exact same as a bed particle.
    
    Keyword arguments:
    bed_particles -- array of all bed particles
    available_vertices -- array of all available vertices in the stream
    set_diam -- diameter of all particles
    pack_fraction -- float packing density value
    h -- 
    
    Return values:
    model_partilces -- array of all model particles
    model_supp -- array of left/right supporting particle ids
                                        for each model particle

     """ 
    num_placement_loc = np.size(available_vertices)
    # determine the number of model particles that should be introduced into the stream bed
    num_particles = determine_num_particles(pack_fraction, num_placement_loc)
    # create an empty n-6 array to store model particle information
    model_particles = np.zeros([num_particles, 7], dtype='float')
    model_supp = np.zeros([num_particles, 2], dtype='float')
  
    for particle in range(num_particles):  
        # the following lines select a vertex to place the current particle at, 
        # and ensure that it is not already occupied by another particle
        random_idx = random.randint(0, np.size(available_vertices)-1)
        vertex = available_vertices[random_idx]
        available_vertices = available_vertices[available_vertices != vertex]

        # intialize the particle information
        model_particles[particle][0] = vertex 
        model_particles[particle][1] = set_diam
        
        model_particles[particle][3] = particle # id number for each particle
        model_particles[particle][4] = 1 # each particle begins as active
        
        # place particle at the chosen vertex
        p_x, p_y, left_supp, right_supp  = place_particle(model_particles[particle], 
                                                            model_particles, 
                                                            bed_particles, 
                                                            h)
        model_particles[particle][0] = p_x
        model_particles[particle][2] = p_y
        model_particles[particle][5] = 0
        model_particles[particle][6] = 0

        model_supp[particle][0] = left_supp
        model_supp[particle][1] = right_supp
    
    return model_particles, model_supp


def compute_available_vertices(model_particles, bed_particles, set_diam, level_limit,
                               lifted_particles=None):
    """ Compute the avaliable vertices in the model stream.

    
    An available vertex is an x location that a 
    model particle is able to be entrained to.
    This function identifies the distinct elevations 
    present in the stream then looks at subsets of 
    particles in decesnding order of their elevation,
    in order to compute available vertices.
    
    For each elevation group, if a particle is sitting on a vertex
    x, then x cannot be available (it is occupied) and it
    is considered nulled. Then vertices v created by two 
    particles touching are considered. If v is not already 
    nulled by a particle occupying the location at a higher
    level, then it is considered an available vertex.
    
    This ends once the bed particles (the lowest elevation)
    have been considered.
    
    
    Keyword arguments: 
        model_particles -- list of model particles
        bed_particles -- list of bed particles
        lifted_particles  -- idx of the 'lifted' particles. Default None
    
    Returns:
        available_vertices -- the set of available vertices
    """
    nulled_vertices = []
    avail_vertices = []
    
    # If we are lifting particles, we need to consider the subset of particles
    # that includes every particles _except_ the particles being 
    if lifted_particles is not None:  
        model_particles_lifted = np.delete(model_particles, 
                                           lifted_particles, 0)
        all_particles = np.concatenate((model_particles_lifted, 
                                        bed_particles), axis=0)
    else:
        all_particles = np.concatenate((model_particles, 
                                        bed_particles), axis=0)
    # Get unique model particle elevations in stream (descending)
    elevations = elevation_list(all_particles[:,2])
    
    for idx, elevation in enumerate(elevations):
        tmp_particles = all_particles[all_particles[:,2] == elevation]
        
        for particle in tmp_particles:    
            nulled_vertices.append(particle[0])
        
        right_vertices = tmp_particles[:,0] + (set_diam / 2)
        left_vertices = tmp_particles[:,0] - (set_diam / 2)
        tmp_shared_vertices = np.intersect1d(left_vertices, right_vertices)
        
        # Enforce level limit by nulling any vertex above limit:
        if len(elevations) == level_limit+1 and idx==0: 
            for vertex in tmp_shared_vertices:
                nulled_vertices.append(vertex)
        
        for vertex in tmp_shared_vertices:
            if vertex not in nulled_vertices:
                avail_vertices.append(vertex)
                
        del(tmp_shared_vertices)
        del(tmp_particles)
        
    available_vertices = np.array(avail_vertices)
    
    return available_vertices


def elevation_list(elevations, desc=True):
    """ Return a sorted list of unique elevation values """
    ue = np.unique(elevations)
    if desc:
           ue = ue[::-1]
    return ue
 
def compute_hops(event_particle_ids, model_particles, mu, sigma, normal=False):
    """ Given a list of (event) paritcles, this function will 
    add a hop distance to current x locations of all event particles. 
    
    Current + hop = desired hop distance.
    
    Hop distances are randomly selected from a log-normal or normal
    distribution. 
    
    Keyword arguments:
        event_particle_ids -- list of event particle ids
        model_particles -- the model's np arry of model_particles
        normal -- boolean flag for sampling from Normal (default Flase)
    
    Returns:
        event_particles -- list of event particles with 'hopped' x-locations
    
    """
    event_particles = model_particles[event_particle_ids]
    if normal:
        s = np.random.normal(mu, sigma, len(event_particle_ids))
    else:
        s = np.random.lognormal(mu, sigma, len(event_particle_ids))
    s_hop = np.round(s, 1)
    s_hop = list(s_hop)
    event_particles[:,0] = event_particles[:,0] + s_hop
    
    return event_particles
 
def move_model_particles(event_particles, model_particles, model_supp, bed_particles, available_vertices, h):
    """ Given an array of event particles and their desired hops, move each
    event particle to the closest valid vertex if its desired hop is not a vertex.  
    Update the model particle and support arrays accordingly.

    Keyword arguments:
        event_particles -- array of particles (full 1-7 struct) to be entrained
        model_particles -- array of all model particles 
        model_supp -- array of left/right supporting particles for each 
                                                            model particle
        bed_particles -- array of all bed particles
        available_particles -- array of available vertices in the stream
    
    Returns:
        model_particles -- array of model particles with event particle updates
        model_supports -- array of left/right supporting particles for each 
                                    model particle, with event particles updated

    """
    # Randomly iterate over event particles
    for particle in np.random.permutation(event_particles):
        orig_x = model_particles[model_particles[:,3] == particle[3]][0][0]
        verified_hop = find_closest_vertex(particle[0], available_vertices)
        
        if verified_hop == -1:
            exceed_msg = (
                f'Particle {int(particle[3])} exceeded stream...'
                f'sending to -1 axis'
            )
            logging.info(exceed_msg) 
            particle[6] = particle[6] + 1
            particle[0] = verified_hop

            model_supp[int(particle[3])][0] = np.nan
            model_supp[int(particle[3])][1] = np.nan
        else:
            hop_msg = (
                f'Particle {int(particle[3])} entrained from {orig_x} '
                f'to {verified_hop}. Desired hop was: {particle[0]}'
            )
            logging.info(hop_msg)
            particle[0] = verified_hop
            available_vertices = available_vertices[available_vertices != verified_hop]

            placed_x, placed_y, left_supp, right_supp = place_particle(particle, model_particles, bed_particles, h)
            particle[0] = placed_x
            particle[2] = placed_y

            model_supp[int(particle[3])][0] = left_supp
            model_supp[int(particle[3])][1] = right_supp

        model_particles[model_particles[:,3] == particle[3]] = particle
    return model_particles, model_supp


def update_flux(initial_positions, final_positions, iteration, subregions):
    """ Given arrays of initial and final positions, this function 
    will update each subregion s to indicate how many paticles crossed
    s's downstream boundary.

    Keyword arguments:
    initial_positions -- array of initial x locations
    final_positions -- array of final (verified) x locations
    iteration -- the iteration the that the crossing should be recorded under
    subregions -- array of subregion objects

    Return values:
    subregions -- array of subregion objects with updated flux values
    """
    # This can _most definitely_ be made quicker but for now, it works
    if len(initial_positions) != len(final_positions):
        raise ValueError(f'Initial_positions and final_positions do not contain the same # of elements')
    
    for position in range(0, len(initial_positions)):

        initial_pos = initial_positions[position]
        final_pos = final_positions[position]

        for idx, subregion in enumerate(subregions):
            if (initial_pos >= subregion.leftBoundary()) and (subregion.rightBoundary() > initial_pos):
                start_idx = idx

        for subregion_idx in range(start_idx, len(subregions)):
            if final_pos >= subregions[subregion_idx].rightBoundary():
                subregions[subregion_idx].incrementFlux(iteration)
            elif final_pos == -1 and subregion_idx == len(subregions)-1:
                subregions[subregion_idx].incrementFlux(iteration)

    return subregions

    
def find_closest_vertex(desired_hop, available_vertices):
    """ Find the closest downstream (greater than or equal) vertex
    in availbale vertices. If nothing exists, then return -1.
    
    Keyword arguments:
    desired_hop -- float representing the desired hop location
    available_location -- np array of available vertices in model
    
    Returns:
    vertex -- the closest available vertex that is >= desired_hop
    """    
    if available_vertices.size == 0:
        raise ValueError('Available vertices array is empty, cannot find closest vertex')
    if desired_hop < 0:
        raise ValueError('Desired hop is negative (invalid)')

    available_vertices = np.sort(available_vertices)
    forward_vertices = available_vertices[available_vertices >= desired_hop]
    
    if forward_vertices.size < 1:
        vertex = -1
    else:
        vertex = forward_vertices[0]
    return vertex    


def increment_age(model_particles, e_event_ids):
    """"Increment model particles' age, set event particles to age 0"""
    
    model_particles[:,5] = model_particles[:,5] + 1 
    model_particles[e_event_ids, 5] = 0
    
    return model_particles
