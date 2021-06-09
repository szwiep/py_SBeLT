import numpy as np
from tqdm import tqdm
from codetiming import Timer
import yaml
import json
import logging
import logging.config
from datetime import datetime
from uuid import uuid4
import shelve

import logic
import util
import sys
import os

ITERATION_HEADER = ("""\n
    --------------------- Iteration {iteration} ---------------------
    """)
   
ITERATION_TEMPLATE = ("""\n
    # of entrainment events: {e_events}\n
    Particles to be entrained: {particles}\n                          
                      """)

def main(run_id, pid):
    
    #############################################################################
    # Set up logging
    #############################################################################

    with open('logs/conf.yaml', 'r') as f:
        config = yaml.safe_load(f.read())
        config['handlers']['file']['filename'] = f'logs/{run_id}.log'
        logging.config.dictConfig(config)
    
    #############################################################################
    # Get and validate parameters
    #############################################################################

    with open('param.yaml', 'r') as p:
        parameters = yaml.safe_load(p.read())
    # TODO: update validation to take dictionary of the parameters
    # util.validate_parameters(parameters)   

    #############################################################################
    #  Create model data and data structures
    #############################################################################

    # Pre-compute d and h values for particle elevation placement
    # see d and h here: https://math.stackexchange.com/questions/2293201/
    d = np.divide(np.multiply(np.divide(parameters['set_diam'], 2), 
                                        parameters['set_diam']), 
                                        parameters['set_diam'])
    h = np.sqrt(np.square(parameters['set_diam']) - np.square(d))


    print(f'[{pid}] Building Bed and  Model particle arrays...')
    # Create bed particle array and compute corresponding available vertices
    bed_particles, bed_length = logic.build_streambed(parameters['x_max'], parameters['set_diam'])   
    available_vertices = logic.compute_available_vertices([], bed_particles, parameters['set_diam'],
                                                        parameters['level_limit'], just_bed=True)    
    # Create model particle array and set on top of bed particles
    model_particles = logic.set_model_particles(bed_particles, available_vertices, parameters['set_diam'], 
                                                        parameters['pack_density'],  h)
    # Define stream's subregions
    subregions = logic.define_subregions(bed_length, parameters['num_subregions'])

    #############################################################################
    #  Create data structures for entrainment iteration information
    #############################################################################

    particle_flux_list = []
    particle_age_list = []
    particle_range_list = []
    snapshot_dict = {}
    snapshot_counter = 0
    milestones = [10, 20, 30, 40, 50, 60, 70, 80, 90, 100]

    snapshot_shelve = shelve.open(f"../plots/run-info-{run_id}")
    # Write parameters and bed (All are static values) to shelf
    snapshot_shelve['param'] = parameters 
    snapshot_shelve['bed'] = np.ndarray.tolist(bed_particles)

    #############################################################################
    #  Entrainment iterations
    #############################################################################

    print(f'[{pid}] Bed and Model particles built. Beginning entrainments...')
    for iteration in range(parameters['n_iterations']):
        snapshot_counter += 1
        # Calculate number of entrainment events iteration
        e_events = np.random.poisson(parameters['lambda_1'], None)
        # Select n (= e_events) particles, per-subregion, to be entrained
        event_particle_ids = logic.get_event_particles(e_events, subregions,
                                                    model_particles, 
                                                    parameters['level_limit'], 
                                                    parameters['height_dependancy'])
        # Determine hop distances of all event particles
        unverified_e = logic.fathel_furbish_hops(event_particle_ids, model_particles, parameters['mu'],
                                                parameters['sigma'], normal=parameters['normal_dist'])
        # Compute available vertices based on current model_particles state
        avail_vertices = logic.compute_available_vertices(model_particles, 
                                                    bed_particles,
                                                    parameters['set_diam'],
                                                    parameters['level_limit'],
                                                    just_bed=False, 
                                                    lifted_particles=event_particle_ids)
        # Run entrainment event                    
        model_particles, particle_flux = logic.run_entrainments(model_particles, 
                                                                bed_particles, 
                                                                event_particle_ids,
                                                                avail_vertices, 
                                                                unverified_e, 
                                                                h)
        # Record number of particles to cross downstream boundary per-iteration                                                        
        particle_flux_list.append(particle_flux)

        # Compute age range and average age, store in lists
        age_range = np.max(model_particles[:,5]) - np.min(model_particles[:,5])
        particle_range_list.append(age_range)

        avg_age = np.average(model_particles[:,5]) 
        particle_age_list.append(avg_age)

        # Record snapshot of relevant iteration information 
        # Currently recording iteration's: 
        #               1) model_particles array 
        #               2) available_vertices used for run_entrainments call
        #               3) event_particle_ids used for run_entrainments cal
        # snapshot_interval determines how often snapshots are stored
        if (snapshot_counter == parameters['snapshot_interval']):
            snapshot_dict[str(iteration)] = [ 
                                        np.ndarray.tolist(model_particles), 
                                        np.ndarray.tolist(available_vertices),
                                        np.ndarray.tolist(event_particle_ids)]
            snapshot_counter = 0

        # Incrementally write snapshot dictionary to file to avoid overwhelming memory
        if(iteration != 0 and iteration % 100000 == 0):
            print(f'[{pid}] Writing chunk of dictionary to shelf...')
            snapshot_shelve.update(snapshot_dict)
            snapshot_dict.clear()
            print(f'[{pid}] Writing of chunk complete. Continuing with entrainments...')

        # Display run progress for users using milestones list
        percentage_complete = (100.0 * (iteration+1) / parameters['n_iterations'])
        while len(milestones) > 0 and percentage_complete >= milestones[0]:
            print(f'[{pid}] {milestones[0]}% complete')
            #remove that milestone from the list
            milestones = milestones[1:]
    
    #############################################################################
    # Store entrainment iteration information
    #############################################################################

    # Once all entrainment events complete, store relevant information to shelf
    # TODO: might restructure to add flux info per-teration like above; plot_snapshot[iteration]
    print(f'[{pid}] Writing dictionary to shelf...')
    snapshot_shelve.update(snapshot_dict)
    print(f'[{pid}] Writting of dictionary complete...')
    
    print(f'[{pid}] Writting flux and age information to shelf...')
    snapshot_shelve['flux'] = particle_flux_list
    snapshot_shelve['avg_age'] = particle_age_list
    snapshot_shelve['age_range'] = particle_range_list
    print(f'[{pid}] Writting of flux and age information complete...')

    # TODO: close needs to go in a finally
    snapshot_shelve.close()
    print(f'[{pid}] Model run complete.')

    #############################################################################

if __name__ == '__main__':

    pid = os.getpid()
    print(f'Process [{pid}] running 1 instance of BeRCM...')
    run_id = datetime.now().strftime('%Y%m-%d%H-%M%S-') + str(uuid4())
    print(f'[{pid}] Using UUID: {run_id}')
    # TODO: make sure UUID/filename has not already been used
    main(run_id, pid)

    
    
