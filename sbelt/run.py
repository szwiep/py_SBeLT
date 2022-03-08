import numpy as np
import yaml
import logging
import logging.config
from datetime import datetime
from pathlib import Path 
import h5py
import time
from tqdm import tqdm
from jsonschema import validate, exceptions

import sbelt.logic as logic
import sys
import os

ITERATION_HEADER = ('Beginning iteration {iteration}...')
ENTRAINMENT_HEADER = ('Entraining particles {event_particles}')

def main():

    # pid = os.getpid()
    # run_id = datetime.now().strftime('%y%m-%d%H')
    # # TODO: send to default config, or get user to edit config, or ask for path to config
    # param_path = sys.argv[1]

    logConf_path, log_path, schema_path, output_path = get_relative_paths()

    #############################################################################
    # Set up logging
    #############################################################################
    
    configure_logging(run_id, logConf_path, log_path)
    
    #############################################################################
    # Get and validate parameters
    #############################################################################
    
    with open(schema_path, 'r') as s:
        schema = yaml.safe_load(s.read())
    with open(param_path, 'r') as p:
        parameters = yaml.safe_load(p.read())
    try:
        validate(parameters, schema)
    except exceptions.ValidationError as e:
        print("Invalid configuration of param file at {param_path}. See the exception below:\n" )
        raise e
    if parameters['x_max'] % parameters['set_diam'] != 0:
        print("Invalid configuration of param file at {param_path}: x_max must be divisible by set_diam.")
        raise ValueError("x_max must be divisible by set_diam")
    if parameters['x_max'] % parameters['num_subregions'] != 0:
        print("Invalid configuration of param file at {param_path}: x_max must be divisible by num_subregions.")
        raise ValueError("x_max must be divisible by num_subregions")

    #############################################################################
    #  Create model data and data structures
    # TODO: Better names for d, h variables
    #############################################################################

    print(f'[{pid}] Building Bed and Model particle arrays...')
    # Pre-compute d and h values for particle elevation placement
        # see d and h here: https://math.stackexchange.com/questions/2293201/
    d = np.divide(np.multiply(np.divide(parameters['set_diam'], 2), 
                                        parameters['set_diam']), 
                                        parameters['set_diam'])
    h = np.sqrt(np.square(parameters['set_diam']) - np.square(d))
    # Build the required structures for entrainment events
    bed_particles, model_particles, model_supp, subregions = build_stream(parameters, h)

    #############################################################################
    #  Create entrainment data and data structures
    #############################################################################

    particle_age_array = np.ones(parameters['n_iterations'])*(-1)
    particle_range_array = np.ones(parameters['n_iterations'])*(-1)
    snapshot_counter = 0

    #############################################################################

    h5py_filename = f'{parameters["filename_prefix"]}-{run_id}.hdf5'
    hdf5_path = f'{output_path}/{h5py_filename}'

    with h5py.File(hdf5_path, "a") as f: 
        
        grp_p = f.create_group(f'params')
        for key, value in parameters.items():
            grp_p[key] = value

        grp_iv = f.create_group(f'initial_values')
        grp_iv.create_dataset('bed', data=bed_particles)
        grp_iv.create_dataset('model', data=model_particles)

        #############################################################################
        #  Entrainment iterations
        #############################################################################

        print(f'[{pid}] Bed and Model particles built. Beginning entrainments...')
        for iteration in tqdm(range(parameters['n_iterations']), leave=False):
            logging.info(ITERATION_HEADER.format(iteration=iteration))
            snapshot_counter += 1

            # Calculate number of entrainment events iteration
            e_events = np.random.poisson(parameters['lambda_1'], None)
            # Select n (= e_events) particles, per-subregion, to be entrained
            event_particle_ids = logic.get_event_particles(e_events, subregions,
                                                        model_particles, 
                                                        parameters['level_limit'], 
                                                        parameters['height_dependancy'])
            logging.info(ENTRAINMENT_HEADER.format(event_particles=event_particle_ids))
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
            # Compute age range and average age, store in np arrays
            age_range = np.max(model_particles[:,5]) - np.min(model_particles[:,5])
            particle_range_array[iteration] = age_range

            avg_age = np.average(model_particles[:,5]) 
            particle_age_array[iteration] = avg_age

            # Record per-iteration information 
            if (snapshot_counter == parameters['data_save_interval']):
                grp_i = f.create_group(f'iteration_{iteration}')
                grp_i.create_dataset("model", data=model_particles, compression="gzip")
                grp_i.create_dataset("event_ids", data=event_particle_ids, compression="gzip")
                snapshot_counter = 0

        #############################################################################
        # Store flux and age information
        #############################################################################
        
        print(f'[{pid}] Writting flux and age information to shelf...')
        grp_final = f.create_group(f'final_metrics')
        grp_sub = grp_final.create_group(f'subregions')
        for subregion in subregions:
            name = f'{subregion.getName()}-flux'
            flux_list = subregion.getFluxList()
            grp_sub.create_dataset(name, data=flux_list, compression="gzip")

        grp_final.create_dataset('avg_age', data=particle_age_array, compression="gzip")
        grp_final.create_dataset('age_range', data=particle_range_array, compression="gzip")
        print(f'[{pid}] Finished writing flux and age information.')

        print(f'[{pid}] Model run finished successfully.')
    return

#############################################################################
# Helper functions
#############################################################################

def build_stream(parameters, h):
    """ Build the data structures which define a stream.       

    Build array of m bed particles and array of n model particles. 
    Model particles will be assigned (x,y) positions which represent 
    resting on top of the bed particles. At the end of the build, each 
    model particle will have 2 support particles from the bed recorded 
    in an array. Finally, define an array of subregion objects. 

    Args:
        parameters: 
            A dictionary (hopefully validated) of the 13 parameters 
            required by the model. For example:

            {'pack_density': (0.78),
            'x_max': (100),
                ...
            'filename_prefix': ('hello-river')}
        h: 
            Geometric value used in calculations of particle placement. See
            in-line and project documentation for further explanation.

    Returns:
        bed_particles: An m-7 NumPy array representing the stream's m bed
            particles. For example:

        [[x, diam, y, uid, active, age, loops], ... ,[x, diam, y, uid, active, age, loops]]
        
        All bed particles should share the same diameter and y (elevation valu), all uids 
        must be negative, and to represent 'static' bed particles active = 0 and loops = 0. 

        model_particles: An n-7 NumPy array representing the stream's n model 
            particles in their initial placement. For example:
        
        [[x, diam, y, uid, active, age, loops], ... ,[x, diam, y, uid, active, age, loops]]

        model_supp: An n-2 NumPy array with the uids of the two particles supporting each 
            model particle. For example:
        
        model_supp = [[[-1,-2]], ... ,[-3,-4]]

        The above example states that model particle with uid 0 (model_supp[0]) is supported
        by bed particles with uids -1 and -2. Similarly, the model particle with uid n 
        (model_supp[n]) is supported by bed particles with uids -3 and -4.

        subregions: An array of Subregion objects
    
    """
    bed_particles = logic.build_streambed(parameters['x_max'], parameters['set_diam'])
    empty_model = np.empty((0, 7))      
    available_vertices = logic.compute_available_vertices(empty_model, bed_particles, parameters['set_diam'],
                                                        parameters['level_limit'])    
    # Create model particle array and set on top of bed particles
    model_particles, model_supp = logic.set_model_particles(bed_particles, available_vertices, parameters['set_diam'], 
                                                        parameters['pack_density'],  h)
    # Define stream's subregions
    subregions = logic.define_subregions(parameters['x_max'], parameters['num_subregions'], parameters['n_iterations'])
    return bed_particles,model_particles, model_supp, subregions


def run_entrainments(model_particles, model_supp, bed_particles, event_particle_ids, avail_vertices, 
                                                                    unverified_e, subregions, iteration, h):
    """ This function mimics a single entrainment event through
    calls to the entrainment-related logic functions. 
    
    An entrainment event consists of:
        (1) Moving a set of the model particles 
            (event particles) downstream.
        (2) Recording crossings/flux across each downstream 
            boundary of the subregions.
        (3) Update model particles states (can it be
            selected for entrainment next iter?) and
            ages are updated 
    
    Args:
        model_particles: An n-7 NumPy array representing the stream's 
            n model particles. 
        model_supp: An n-2 NumPy array with the uids of the two particles supporting each 
            model particle (e.g model_supp[j] = supports for model particle j). 
        bed_particles: An m-7 NumPy array representing the stream's m 
            bed particles.
        event_particle_ids: A NumPy array of k uids representing the model particles
            that have been selected for entrainment.
        
    Returns:
        model_particles: Updated model_particles (Args) with updated age, location, 
            loops, and states, based on entrainment event placements.
        model_supp: An updated model_supports (Args) based on placements.
        subregions: Python array of Subregion objects with updated flux lists.
    """

    initial_x = model_particles[event_particle_ids][:,0]
    model_particles, model_supp = logic.move_model_particles(unverified_e, 
                                                                model_particles,
                                                                model_supp, 
                                                                bed_particles, 
                                                                avail_vertices,
                                                                h)
    final_x = model_particles[event_particle_ids][:,0]
    subregions = logic.update_flux(initial_x, final_x, iteration, subregions)
    model_particles = logic.update_particle_states(model_particles, model_supp)
    # Increment age at the end of each entrainment
    model_particles = logic.increment_age(model_particles, event_particle_ids)

    return model_particles, model_supp, subregions


def configure_logging(run_id, logConf_path, log_path):
    """Configure logging procedure using conf.yaml"""
    with open(logConf_path, 'r') as f:
        config = yaml.safe_load(f.read())
        config['handlers']['file']['filename'] = f'{log_path}/{run_id}.log'
        logging.config.dictConfig(config)


def get_relative_paths():
    """Return relavent relative paths in the project"""
    logConf_path = Path(__file__).parent / 'logs/conf.yaml'
    log_path =  Path(__file__).parent / 'logs/'
    schema_path = Path(__file__).parent / 'parameters/schema.yaml'
    output_path = Path(__file__).parent / 'output/'
    return logConf_path,log_path, schema_path, output_path

    
    
