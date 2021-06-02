import numpy as np
from tqdm import tqdm
from codetiming import Timer
import yaml
import json
import logging
import logging.config
from datetime import datetime
from uuid import uuid4

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

def main(run_id):
    
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

    # Pre-compute d and h values for particle elevation placement
    # see d and h here: https://math.stackexchange.com/questions/2293201/
    d = np.divide(np.multiply(np.divide(parameters['set_diam'], 2), 
                                        parameters['set_diam']), 
                                        parameters['set_diam'])
    h = np.sqrt(np.square(parameters['set_diam']) - np.square(d))


    print('Building Bed and  Model particle arrays...')
    # Create bed particle array and compute corresponding available vertices
    bed_particles, bed_length = logic.build_streambed(parameters['x_max'], parameters['set_diam'])   
    available_vertices = logic.compute_available_vertices([], bed_particles, parameters['set_diam'],
                                                        parameters['level_limit'], just_bed=True)    
    # Create model particle array and set on top of bed particles
    model_particles = logic.set_model_particles(bed_particles, available_vertices, parameters['set_diam'], 
                                                        parameters['pack_density'],  h)
    # Define stream's subregions
    subregions = logic.define_subregions(bed_length, parameters['num_subregions'])

    particle_flux_list = []
    plot_snapshot = {}
    plot_snapshot['param'] = parameters 
    plot_snapshot['bed'] = np.ndarray.tolist(bed_particles)
    snapshot_counter = 0

    print('Bed and Model particles built. Beginning entrainments...')
    for iteration in tqdm(range(parameters['n_iterations'])):
        snapshot_counter += 1
        # Calculate number of entrainment events iteration
        e_events = np.random.poisson(parameters['lambda_1'], None)
        # Randomly select n (= e_events) particles, per-subregion, to be entrained
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

        # Record snapshot of relevant iteration information 
        # Currently recording iteration's: 
        #               1) model_particles array 
        #               2) available_vertices used for run_entrainments call
        #               3) event_particle_ids used for run_entrainments call
        if (snapshot_counter == parameters['snapshot_interval']):
            plot_snapshot[iteration] = [ 
                                        np.ndarray.tolist(model_particles), 
                                        np.ndarray.tolist(available_vertices),
                                        np.ndarray.tolist(event_particle_ids)]
            snapshot_counter = 0
    
    # Store final list of flux information
    # TODO: might restructure to add flux info per-teration like above; plot_snapshot[iteration]
    plot_snapshot['flux'] = particle_flux_list
    print('Model run complete...')

    #############################################################################
    # Time profiling
    #############################################################################

    print(Timer.timers) 

    #############################################################################

    print('Writing iteration snapshots to file...')
    plot_snapshot_jsn = json.dumps(plot_snapshot)

    print("Estimated size: " + str(sys.getsizeof(plot_snapshot_jsn) / 1024) + "KB")

    outfilename = "../plots/run-info-" + run_id + ".json"
    with open(outfilename, 'w') as outfile:
        outfile.write(plot_snapshot_jsn)
  
    print("Actual size: " + str(os.path.getsize(outfilename) / 1024) + "KB")

    # print('Plotting flux information...')
    # plot.flux_info(particle_flux_list, parameters['n_iterations'], to_file=True)
    print('Model execution complete.')

    #############################################################################

if __name__ == '__main__':

    run_id = datetime.now().strftime('%Y%m-%d%H-%M%S-') + str(uuid4())

    main(run_id)

    
    
