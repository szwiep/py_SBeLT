""" 
This module contains all the logic required to execute
sbelt simulations/runs. All functions/classes are designed 
for internal use and may change without note.

A primary purpose of this module is to manipulate n-7 NumPy 
arrays which represent multiple stream particles. In these
n-7 arrays, a _single_ particle is represented by a NumPy array 
with 7 attributes::

    [x_location, diam, y_location, UID, active state, age, loop age]

In this general example, x_location and y_location define the centre point of 
spherical particle whose diameter is defined by diam. See the project docs for
more information on the other attributes. 
"""

import math
import random
import numpy as np


import logging
logging.getLogger(__name__)


class Subregion():
    """ A subregion in the stream.
    
    A Subregion is defined by left (upstream)
    and right (downstream) boundaries. Each
    subregion maintains a NumPy list which is 
    used to record the number of model particles that 
    pass the downstream boundary in a given iteration. 
    For example::

        flux_list = [0,3,2,1]

    Means that 0 crossings happened in the first iteration, 
    3 happened in the 2nd, and so on. The list has 
    length equal to the number of iterations for a model run.
    
    Attributes:
        name: Name of the subregion.
        left_boundary: Location of the left boundary (float).
        right_boundary: Location of the right boundary (float).
        iterations: The number of iterations for the model run.
    """
    def __init__(self, name, left_boundary, right_boundary, iterations):
        self.name = name
        self.left_boundary = left_boundary
        self.right_boundary = right_boundary
        self.flux_list = np.zeros(iterations, dtype=np.int64)
        
    def leftBoundary(self):
        """Returns subregion's left boundary"""
        return self.left_boundary
    
    def rightBoundary(self):
        """Returns subregion's right boundary"""
        return self.right_boundary
    
    def getName(self):
        """Returns subregion's name"""
        return self.name

    def incrementFlux(self, iteration):
        """Increments flux list by 1.

        Args:
            iteration: The iteration/index to increment by 1
        """
        self.flux_list[iteration] += 1
    
    def getFluxList(self):
        """Returns subregion's flux list"""
        return self.flux_list

def get_event_particles(e_events, subregions, model_particles, level_limit, height_dependant=False):
    """ Find and return list of particles to be entrained

    Will loop through each subregion and select n = e_events
    model particles (within a subregion boundaries) at random 
    to be entrained. No particle will be selected twice.

    Args:
        e_events: The number of events requested per subregion (int).
        subregions: Python array of initialized Subregion objects. 
        model_particles: An n-7 NumPy array representing the stream's n 
            model particles.

    Returns:
        event_particles: A NumPy array of k uids representing the model particles
            that have been selected for entrainment. For example::

                [2.0, 5.0, 25.0]

            Will represent that model particles with uids 2.0, 5.0 and
            25.0 have been selected for entrainment. 
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
            logging.info(msg)
        event_particles = event_particles + subregion_event_ids
    event_particles = np.array(event_particles, dtype=np.intp)

    return event_particles

def define_subregions(bed_length, num_subregions, iterations):
    """ Define subregion list for model stream.
    
    Args:
        bed_length: The length of the stream (int). 
        num_subregions: The number of subregions (int).
        iterations: The number of iterations for the model run (int).

    Returns:
        subregions_arr: Python array of initialized Subregion objects. 
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
    
def build_streambed(bed_length, particle_diam):
    """ Builds the array of bed particles.
    
    Args:
        bed_length: The length of the stream (int).
        particle_diam: The diameter of all particles (float).
    
    Returns:
        bed_particles: An m-7 NumPy array representing the stream's m 
            bed particles. For example::

                [[x1, diam, y, uid1, active, age, loops], ... 
                            ,[xM, diam, y, uidM, active, age, loops]]
            
            Where all bed particles share the same diam (diam1=...=diamM) and 
            y (y1=...=yM), all uids are negative, and to represent 'static-ness'
            active = 0, age = 0, and loops = 0. 

    """
    max_particles = int(math.ceil( bed_length / particle_diam ))
    bed_particles = np.zeros([max_particles, 7],dtype=float)
    
    particle_id = -1
    centre = (particle_diam/2)  
    state = 0
    age = 0
    loop_age = 0
    elevation = 0
    while not bed_complete(centre, bed_length):  
        # index with negative indices... bed particles are built from the final element to the first
        bed_particles[particle_id] = [centre, particle_diam, elevation, particle_id, state, age, loop_age]
        centre += particle_diam
        particle_id += -1 # Bed particles get negative IDs
    
    return bed_particles

def bed_complete(centre, bed_length):
    if centre >= bed_length:
        return 1
    else: return 0


def determine_num_particles(pack_frac, num_vertices):
    """Return the number of model particles to be created
    based on the packing fraction"""
    
    num_particles = num_vertices * pack_frac
    num_particles = int(math.ceil(num_particles))
    
    return num_particles

# Trig from: https://math.stackexchange.com/questions/2293201/
def place_particle(particle, model_particles, bed_particles, h):
    """ Calculate new y (elevation) of particle based on it's 
        x (horizontal) location in stream.
    
    Provided a particle p's location x in the stream, 
    search for 2 supporting particles (s1, s2) that p 
    will rest on when placed at x.
    
    Calculate the y position of p based on the (x,y) of both 
    s1 and s2. The computed x for p might be different up to some 
    decimal point, so both x and y are rounded to 2 decimal places.

    
    Args:
        particle: NumPy array representing the model particle being placed.
        model_particles: An n-7 NumPy array representing the stream's n model particles.
        bed_particles: An m-7 NumPy array representing the stream's m bed particles.

    Returns:
        rounded_x: Rounded float of particle's new x location.
        rounded_y: Rounded float of particle's new y location.
        left_support: UID of the left support for the placed particle.
        right_support: UID of the right support for the placed particle.
    
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
    
    Args:
        model_particles: An n-7 NumPy array representing the stream's 
            n model particles.
        model_supports: An n-2 NumPy array with the uids of the two 
            particles supporting each model particle.
        bed_particles: An m-7 NumPy array representing the stream's m bed
            particles.

    Return values:
        model_particles: The provided model_particles array (Args) 
            but with updated active (attribute 4) values.
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
    rest on, if being dropped at x. More generally, supporting particles 
    are those particles which directly hold up a particle. Supporting
    particles will always have a centre location that is 
    exactly a radius length away from p's centre.

    Args:
        particle: 1-7 NumPy array representing a model particle.
        model_particles: An n-7 NumPy array representing the stream's 
            n model particles.
        bed_particles: An m-7 NumPy array representing the stream's m 
            bed particles.

    Returns:
        left_support: NumPy array representing the left supporting particle
        right_support: NumPy array representing the right supporting particle
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


def set_model_particles(bed_particles, available_vertices, particle_diam, pack_fraction, h):
    """ Create array of n model particles and set each particle in-stream.
    
    Model particles are randomly placed at available vertex
    locations (x,y) across the bed. Location and initial attribute
    values are stored in the returned NumPy array. 
    
    Args:
        bed_particles: An m-7 NumPy array representing the stream's m bed particles.
        available_vertices: A NumPy array with all available vertices in the stream. 
        particle_diam: The diameter of all particles (float).
        pack_fraction: Packing density value (float). See THEORY.md in project
            repo for more information.
        h: Geometric value used in calculations of particle placement (float). See
            in-line and project documentation for further explanation.
    
    Returns:
        model_particles: An n-7 NumPy array representing the stream's n model particles and their 
            initial placement in the stream. For example::
            
                [[x1, diam, y1, uid1, active, age, loops], ... ,[xN, diam, yN, uidN, active, age, loops]]

            Where (xi, yi) pairs will define the centre location of each particle and no two particles
            will have the same (xi, yi) pair values. All uids are unique and are positive whole numbers.
            All particles will start with active = 1, age = 0, and loops = 0.

        model_supp: An n-2 NumPy array with the uids of the two particles supporting each 
            model particle. For example::
        
                [[[-1,-2]], ... ,[-3,-4]]

            The above example states that model particle with uid 0 (model_supp[0]) is supported
            by bed particles with uids -1 and -2. Similarly, the model particle with uid n 
            (model_supp[n]) is supported by bed particles with uids -3 and -4. 
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
        model_particles[particle][1] = particle_diam
        
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


def compute_available_vertices(model_particles, bed_particles, particle_diam, level_limit,
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
    level, then it is considered an available vertex. This ends once the bed 
    particles (the lowest elevation) have been considered.
    

    Args:
        model_particles: An n-7 NumPy array representing the stream's 
            n model particles.
        bed_particles: An m-7 NumPy array representing the stream's m 
            bed particles. 
        lifted_particles: UID of particles that are 'lifted'. Lifted 
            particles will not be considered as present in the stream
            when the available vertices are being calculated; their
            (x,y) location will not be considered as occupied.
    
    Returns:
        available_vertices: A NumPy array with all available vertices in the stream. 
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
        
        right_vertices = tmp_particles[:,0] + (particle_diam / 2)
        left_vertices = tmp_particles[:,0] - (particle_diam / 2)
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
    """ Returns a sorted list of unique elevation values """
    ue = np.unique(elevations)
    if desc:
           ue = ue[::-1]
    return ue
 
def compute_hops(event_particle_ids, model_particles, mu, sigma, normal=False):
    """ Given a list of event paritcles, this function will 
    add a hop distance to current x locations of all event particles. 
    
    Current + Hop = desired hop distance. Hop values are randomly 
    selected from a log-normal or normal distribution. For more
    information on the use of these distributions see THEORY.md 
    in the project repository.
    
    Args:
        event_particle_ids: A NumPy array of k uids representing the model particles
            that have been selected for entrainment.
        model_particles: An n-7 NumPy array representing the stream's 
            n model particles.
        normal (default = False): Boolean flag for which distribution to sample
            Hop values from. True = sample from Normal, False = sample from log-Normal. 
    
    Returns:
        event_particles: A k-7 Numpy array representing each event particle
            with updated x locations (x=deried hop location).
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
    """ Move model particles in the stream.
    
    Given an array of event particles and their desired hops, move each
    event particle in-stream based on the event particle's desired hop. Sometimes
    the desired hop will equal an available vertex and sometimes it will not.
    When it does not, move the particle to the closest vertex to the 
    desired hop in the downstream direction.

    Args:
        event_particles: A k-7 Numpy array representing event particle.
        model_particles: An n-7 NumPy array representing the stream's 
            n model particles.
        model_supp: An n-2 NumPy array with the uids of the two particles supporting each 
            model particle (e.g model_supp[j] = supports for model particle j).  
        bed_particles: An m-7 NumPy array representing the stream's m 
            bed particles. 
        available_vertices: A NumPy array with all available vertices in the stream. 
        h: Geometric value used in calculations of particle placement (float). See
            in-line and project documentation for further explanation.
    
    Returns:
        model_particles: The provided model_particles array (Args) 
            but with updated (x,y) attributes based on placements. 
        model_supports: An updated model_supports (Args) based on 
            placements. Note that only event particles will ever 
            have their model supports updated.
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
    """ Update flux lists in each subregion.

    Given arrays of initial and final positions, this function 
    will update each subregion's flux list to indicate how many 
    paticles crossed the subregions's downstream boundary in 
    a given iteration.

    Args:
        initial_positions: NumPy array of initial x locations.
        final_positions: NumPy array of final (verified) x locations.
        iteration: The iteration for the flux to be updated for (int).
        subregions: Python array of Subregion objects.

    Returns:
        subregions: Python array of Subregion objects with updated flux lists.
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
    in availbale vertices. If nothing exists, return -1.
    
    Args:
        desired_hop: Float representing the desired hop location
        available_location: A NumPy array with all available vertices in the stream. 
    
    Returns:
        vertex: The location of the closest available vertex that 
            is >= desired_hop (float).
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
    """Increment model particles' age, set event particles to age 0"""
    
    model_particles[:,5] = model_particles[:,5] + 1 
    model_particles[e_event_ids, 5] = 0
    
    return model_particles
