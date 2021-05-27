import numpy as np
from tqdm import tqdm
from codetiming import Timer
import yaml
import json
import logging
import time

import logic as ml 
import plotting as plot
import util

ITERATION_HEADER = ("""\n
    --------------------- Iteration {iteration} ---------------------
    """)
   
ITERATION_TEMPLATE = ("""\n
    # of entrainment events: {e_events}\n
    Particles to be entrained: {particles}\n                          
                      """)

def main():
    plots_flag, hp_flag = util.parse_arguments()

    #############################################################################
    # Set up logging
    #############################################################################

    import logging
    import logging.config

    with open('logs/conf.yaml', 'r') as f:
        config = yaml.safe_load(f.read())
        logging.config.dictConfig(config)
    
    #############################################################################
    # Get and validate parameters
    #############################################################################

    with open('param.yaml', 'r') as p:
        parameters = yaml.safe_load(p.read())
    # util.validate_parameters(parameters)   

    #############################################################################

    # parameters = (
    #     f"Pack = {pm.Pack}, x_max = {pm.x_max}, set_diam = {pm.set_diam},"
    #     f"num_subregions = {pm.num_subregions}, level_limit = {pm.level_limit},"
    #     f"n_iterations = {pm.n_iterations}, lambda_1 = {pm.lambda_1},"
    #     f"normal_dist = {pm.normal_dist}, mu = {pm.mu}, sigma = {pm.sigma} \n" 
    # )  

    #############################################################################

    logging.info('Beginning model run...')
    # Pre-compute d and h values for particle elevation placement
    d = np.divide(np.multiply(np.divide(parameters['set_diam'], 2), 
                                        parameters['set_diam']), 
                                        parameters['set_diam'])
    h = np.sqrt(np.square(parameters['set_diam']) - np.square(d))

    bed_particles, bed_length = ml.build_streambed()   
    model_particles = ml.set_model_particles(bed_particles, h)

    subregions = ml.define_subregions(bed_length, parameters['num_subregions'])

    particle_flux_list = []
    plot_snapshot = {}
    plot_snapshot['param'] = parameters 
    snapshot_counter = 0
    for iteration in tqdm(range(parameters['n_iterations'])):
        snapshot_counter += 1
        # print(ITERATION_HEADER.format(iteration=iteration))  
        # Calculate number of entrainment events per-subregion for iteration
        e_events = np.random.poisson(parameters['lambda_1'], None)
        # Retrieve the Event Particles
        event_particles = ml.get_event_particles(e_events, subregions,
                                                model_particles, parameters['level_limit'], hp_flag)
        # print(ITERATION_TEMPLATE.format(
        #                             e_events=len(event_particles), 
        #                             particles=event_particles))     
        
        model_particles, particle_flux = ml.run_entrainments(model_particles, 
                                                            bed_particles, 
                                                            event_particles, 
                                                            parameters['normal_dist'], 
                                                            h)
        particle_flux_list.append(particle_flux)

        
        available_vertices = ml.compute_available_vertices(model_particles, 
                                                        bed_particles)
        if (snapshot_counter == parameters['snapshot_interval']):
            plot_snapshot[iteration] = [np.ndarray.tolist(bed_particles), 
                                        np.ndarray.tolist(model_particles), 
                                        np.ndarray.tolist(available_vertices)]
            snapshot_counter = 0
    logging.info('Model run complete...')

    #############################################################################

    print(Timer.timers) 

    logging.info('Writing iteration information to file...')
    timestr = time.strftime("%Y%m%d-%H%M%S")
    jsn = json.dumps(plot_snapshot)
    f = open("../plots/run-info-" + timestr + ".json", "w")
    f.write(jsn)
    f.close()

    logging.info('Plotting flux information...')
    plot.flux_info(particle_flux_list, to_file=True)
    logging.info('Model execution complete.')

    #############################################################################

if __name__ == '__main__':
    main()

    
    
