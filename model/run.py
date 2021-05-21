import numpy as np
import argparse
from tqdm import tqdm
from codetiming import Timer
import json
import logging
import time

import model_logic as ml 
import parameters as pm
import plotting as plot

ITERATION_HEADER = ("""\n
    --------------------- Iteration {iteration} ---------------------
    """)
   
ITERATION_TEMPLATE = ("""\n
    # of entrainment events: {e_events}\n
    Particles to be entrained: {particles}\n                          
                      """)


# Parse command line arguments and validate parameters
parser = argparse.ArgumentParser()
parser.add_argument('--no-plots', action='store_true')
args = parser.parse_args()
no_plots = args.no_plots

ml.validate_parameters()   
parameters = (
    f"Pack = {pm.Pack}, x_max = {pm.x_max}, set_diam = {pm.set_diam},"
    f"num_subregions = {pm.num_subregions}, level_limit = {pm.level_limit},"
    f"n_iterations = {pm.n_iterations}, lambda_1 = {pm.lambda_1},"
    f"normal_dist = {pm.normal_dist}, mu = {pm.mu}, sigma = {pm.sigma} \n" 
)  
    
# Run the model
d = np.divide(np.multiply(np.divide(pm.set_diam, 2), pm.set_diam), pm.set_diam)
h = np.sqrt(np.square(pm.set_diam) - np.square(d))

bed_particles, bed_length = ml.build_streambed()   
model_particles = ml.set_model_particles(bed_particles, h)

subregions = ml.define_subregions(bed_length, pm.num_subregions)

particle_flux_list = []
plot_snapshot = {}
plot_snapshot['param'] = parameters 
for iteration in tqdm(range(pm.n_iterations)):
    
    # print(ITERATION_HEADER.format(iteration=iteration))  
    # Calculate number of entrainment events per-subregion for iteration
    e_events = np.random.poisson(pm.lambda_1, None)
    # Retrieve the Event Particles
    event_particles = ml.get_event_particles(e_events, subregions,
                                             model_particles)
    # print(ITERATION_TEMPLATE.format(
    #                             e_events=len(event_particles), 
    #                             particles=event_particles))     
    
    model_particles, particle_flux = ml.run_entrainments(model_particles, 
                                                         bed_particles, 
                                                         event_particles, 
                                                         pm.normal_dist, 
                                                         h)
    particle_flux_list.append(particle_flux)
    
    # Re-calculate avail vertices for plotting.
    available_vertices = ml.compute_available_vertices(model_particles, 
                                                       bed_particles)

    plot_snapshot[iteration] = [np.ndarray.tolist(bed_particles), np.ndarray.tolist(model_particles), np.ndarray.tolist(available_vertices)]

print(Timer.timers) 

timestr = time.strftime("%Y%m%d-%H%M%S")

print("Writing plotting information to file...")
json = json.dumps(plot_snapshot)
f = open("../plots/run-info-" + timestr + ".json", "w")
f.write(json)
f.close()
print("Finished writing plotting information to file...")

plot.flux_info(particle_flux_list, to_file=True)
print("\n\nModel run complete\n\n")

    
    
