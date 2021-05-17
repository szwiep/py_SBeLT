import numpy as np
import argparse
from tqdm import tqdm
from codetiming import Timer


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
    
# Run the model
bed_particles, bed_length = ml.build_streambed()   
model_particles = ml.set_model_particles(bed_particles)

subregions = ml.define_subregions(bed_length, pm.num_subregions)

particle_flux_list = []

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
                                                         pm.normal_dist)
    particle_flux_list.append(particle_flux)
    
    # Re-calculate avail vertices for plotting.
    available_vertices = ml.compute_available_vertices(model_particles, 
                                                       bed_particles)
    if no_plots:
        continue
    # TODO: Store iteration data/state, allow call to plots only if desired
    else:
        plot.stream(iteration, bed_particles, model_particles, pm.x_max, 10, 
                        available_vertices, to_file=True)
print(Timer.timers)
    
plot.flux_info(particle_flux_list, to_file=True)
    
    
