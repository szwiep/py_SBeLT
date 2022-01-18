import numpy as np
import yaml
import logging
import logging.config
from datetime import datetime
from pathlib import Path 
from shortuuid import uuid
import h5py
import cProfile
import time

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

def main(run_id, pid, param_path):

    logConf_path, log_path, output_path = get_relative_paths()

    #############################################################################
    # Set up logging
    #############################################################################
    
    configure_logging(run_id, logConf_path, log_path)
    
    #############################################################################
    # Get and validate parameters
    # TODO: update validation to work with yaml
    #############################################################################
    
    with open(param_path, 'r') as p:
        parameters = yaml.safe_load(p.read())
    # util.validate_parameters(parameters)   

    #############################################################################
    #  Create model data and data structures
    # TODO: Better names for d, h variables
    #############################################################################

    # Pre-compute d and h values for particle elevation placement
    # see d and h here: https://math.stackexchange.com/questions/2293201/
    d = np.divide(np.multiply(np.divide(parameters['set_diam'], 2), 
                                        parameters['set_diam']), 
                                        parameters['set_diam'])
    h = np.sqrt(np.square(parameters['set_diam']) - np.square(d))


    print(f'[{pid}] Building Bed and Model particle arrays...')
    bed_particles, model_particles, model_supp, subregions = build_stream(parameters, h)

    #############################################################################
    #  Create entrainment data and data structures
    #############################################################################

    particle_age_list = []
    particle_range_list = []
    snapshot_dict = {}
    snapshot_counter = 0
    milestones = [10, 20, 30, 40, 50, 60, 70, 80, 90, 100]

    # snapshot_shelve = prepare_output_shelve(run_id, output_path, param_path, parameters)
    h5py_filename = f'{parameters["filename_prefix"]}-{run_id}.hdf5'
    hdf5_path = f'{output_path}/{h5py_filename}'
    with h5py.File(hdf5_path, "a") as f: 
        
        grp_iv = f.create_group(f'initial_values')
        for key, value in parameters.items():
            grp_iv[key] = value
        grp_iv.create_dataset('bed', data=bed_particles)
        grp_iv.create_dataset('model', data=model_particles)

        #############################################################################
        #  Entrainment iterations
        #############################################################################

        print(f'[{pid}] Bed and Model particles built. Beginning entrainments...')
        for iteration in range(parameters['n_iterations']):
            logging.info(ITERATION_HEADER.format(iteration=iteration))
            snapshot_counter += 1
            # Calculate number of entrainment events iteration
            e_events = np.random.poisson(parameters['lambda_1'], None)
            # Select n (= e_events) particles, per-subregion, to be entrained
            event_particle_ids = logic.get_event_particles(e_events, subregions,
                                                        model_particles, 
                                                        parameters['level_limit'], 
                                                        parameters['height_dependancy'])
            # Determine hop distances of all event particles
            unverified_e = logic.compute_hops(event_particle_ids, model_particles, parameters['mu'],
                                                    parameters['sigma'], normal=parameters['normal_dist'])
            # Compute available vertices based on current model_particles state
            avail_vertices = logic.compute_available_vertices(model_particles, 
                                                        bed_particles,
                                                        parameters['set_diam'],
                                                        parameters['level_limit'],
                                                        lifted_particles=event_particle_ids)
            # Run entrainment event                    
            model_particles, model_supp, subregions = run_entrainments(model_particles, 
                                                                    model_supp,
                                                                    bed_particles, 
                                                                    event_particle_ids,
                                                                    avail_vertices, 
                                                                    unverified_e,
                                                                    subregions,
                                                                    iteration,  
                                                                    h)
            # Record number of particles to cross downstream boundary per-iteration                                                        
            # particle_flux_list.append(particle_flux)

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
            if (snapshot_counter == parameters['snapshot_interval']):
                # print(event_particle_ids)
                grp_i = f.create_group(f'iteration_{iteration}')
                grp_i.create_dataset("model", data=model_particles, compression="gzip")
                # grp.create_dataset("vertices", data=vertices_array, compression="gzip")
                grp_i.create_dataset("event_ids", data=event_particle_ids, compression="gzip")
                snapshot_counter = 0

            # Incrementally write snapshot dictionary to file to avoid overwhelming memory
            # if(iteration != 0 and iteration % 100000 == 0):
            #     print(f'[{pid}] Writing chunk of dictionary to shelf...')
            #     snapshot_shelve.update(snapshot_dict)
            #     snapshot_dict.clear()
            #     print(f'[{pid}] Finished writing chunk. Continuing with entrainments...')

            # Display run progress for users
            percentage_complete = (100.0 * (iteration+1) / parameters['n_iterations'])
            while len(milestones) > 0 and percentage_complete >= milestones[0]:
                print(f'[{pid}] {milestones[0]}% complete')
                #remove that milestone from the list
                milestones = milestones[1:]

        #############################################################################
        # Store final entrainment iteration information
        #############################################################################
        
        print(f'[{pid}] Writting flux and age information to shelf...')

        grp_final = f.create_group(f'final_metrics')
        grp_sub = grp_final.create_group(f'subregions')
        for subregion in subregions:
            name = f'{subregion.getName()}-flux'
            flux_list = subregion.getFluxList()
            grp_sub.create_dataset(name, data=flux_list, compression="gzip")

        # snapshot_shelve['flux'] = particle_flux_list
        grp_final.create_dataset('avg_age', data=particle_age_list, compression="gzip")
        grp_final.create_dataset('age_range', data=particle_range_list, compression="gzip")
        print(f'[{pid}] Finished writing flux and age information.')
        print(f'[{pid}] Model run finished successfully.')
    
    return

#############################################################################
# Helper functions
#############################################################################

def build_stream(parameters, h):
    bed_particles, bed_length = logic.build_streambed(parameters['x_max'], parameters['set_diam'])
    empty_model = np.empty((0, 7))      
    available_vertices = logic.compute_available_vertices(empty_model, bed_particles, parameters['set_diam'],
                                                        parameters['level_limit'])    
    # Create model particle array and set on top of bed particles
    model_particles, model_supp = logic.set_model_particles(bed_particles, available_vertices, parameters['set_diam'], 
                                                        parameters['pack_density'],  h)
    # Define stream's subregions
    subregions = logic.define_subregions(bed_length, parameters['num_subregions'], parameters['n_iterations'])
    return bed_particles,model_particles, model_supp, subregions

# So, so many parameters.
def run_entrainments(model_particles, model_supp, bed_particles, event_particle_ids, avail_vertices, unverified_e, subregions, iteration, h):
    """ This function mimics an 'entrainment event' through
    calls to the entrainment-related functions. 
    
    Uniqueness of entrainments is forced post-event. Particles
    which select non-unique entrainments will be forced
    to re-entrain at another vertex.
    
    Keyword arguments:
        model_particles -- model's model particles np array
        bed_particles -- model's bed particle np array
        event_particle_ids -- list of ids of event particles
        
    Returns:
        model_particles -- updated model_particles array
        particle_flux -- number (int) of particles which 
                            passed the downstream boundary
    """

    initial_x = model_particles[event_particle_ids][:,0]

    e_dict, model_particles, model_supp, avail_vertices = logic.move_model_particles(
                                                unverified_e, 
                                                model_particles,
                                                model_supp, 
                                                bed_particles, 
                                                avail_vertices,
                                                h)
    unique_entrainments, redo_ids = logic.check_unique_entrainments(e_dict)
     
    while not unique_entrainments:
        redo_entrainments = model_particles[np.searchsorted(model_particles[:,3], 
                                                            redo_ids)]
        e_dict, model_particles, model_supp, avail_vertices = logic.move_model_particles(
                                                            redo_entrainments, 
                                                            model_particles, 
                                                            model_supp,
                                                            bed_particles, 
                                                            avail_vertices,
                                                            h)
        unique_entrainments, redo_ids = logic.check_unique_entrainments(e_dict)

    final_x = model_particles[event_particle_ids][:,0]

    subregions = logic.update_flux(initial_x, final_x, iteration, subregions)
    model_particles = logic.update_particle_states(model_particles, model_supp, bed_particles)
    # Increment age at the end of each entrainment
    model_particles = logic.increment_age(model_particles, event_particle_ids)

    return model_particles, model_supp, subregions

def prepare_output_shelve(run_id, output_path, param_path, parameters):
    # Temporary logic to add 'simX' prefix to outputs
    filename = f'sim{parameters["filename_prefix"]}-{run_id}'
    if not os.path.exists(output_path):
        os.makedirs(output_path)
    snapshot_shelve = shelve.open(f"{output_path}/{filename}")
    return snapshot_shelve


def configure_logging(run_id, logConf_path, log_path):
    with open(logConf_path, 'r') as f:
        config = yaml.safe_load(f.read())
        config['handlers']['file']['filename'] = f'{log_path}/{run_id}.log'
        logging.config.dictConfig(config)


def get_relative_paths():
    logConf_path = Path(__file__).parent / 'logs/conf.yaml'
    log_path =  Path(__file__).parent / 'logs/'
    output_path = Path(__file__).parent / 'output/'
    return logConf_path,log_path,output_path


if __name__ == '__main__':

    # pr = cProfile.Profile()
    # pr.enable()
    tic = time.perf_counter()
    uid = uuid()
    pid = os.getpid()
    run_id = datetime.now().strftime('%y%m-%d%H-') + uid
    print(f'Process [{pid}] run ID: {run_id}')
    
    main(run_id, pid, sys.argv[1])
    toc = time.perf_counter()
    print(f"Completed in {toc - tic:0.4f} seconds")

    # pr.disable()
    # pr.dump_stats('profile_dump_nov_24')

    
    
