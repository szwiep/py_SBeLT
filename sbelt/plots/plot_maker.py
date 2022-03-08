import numpy as np
import h5py
import argparse
import os

import plotting


def main(filename, save_location, iter_min, iter_max):
    # Path to run-info file
    fp_in = filename
    # Path to output directory
    fp_out = save_location + '/'
    # Range to plot
    iteration_range = [iter_min, iter_max] 
    # Plot height
    y_limit = 10

    # Make appropriate sub-directory if it doesn't exist already
    # Don't overwrite data without user's permission
    try: 
        os.mkdir(save_location)
        print(f'Creating output directory...')
    except FileExistsError:
        while True: 
            response = input(f'Directory {save_location}/ already exists. You could be overwriting existing data, continue (Y/N)? ')
            if response.lower() == 'y':
                print(f'Continuing with plotting script...')
                break 
            elif response.lower() == 'n':
                print('Exiting plotting script...')
                return 

    with h5py.File(fp_in, 'r') as f:
        subregion_max = f['params']['num_subregions'][()]
        ds_boundary_key = f'subregion-{subregion_max-1}-flux' 
        # Plot stream bed for iter_min to iter_max
        #   This method assumes plotting iteration 0 means plotting
        #   the model particles at the start of iteration 0. Since the
        #   model particles are stored (runpy) at the end of each iteration, 
        #   we need to plot 1 index back.
        print(f'Plotting stream bed...')
        for iter in range(iteration_range[0], iteration_range[1]+1):
            if iter == 0:
                model_particles = np.array(f['initial_values']['model'])
            else:
                model_particles = np.array(f[f'iteration_{iter-1}']['model'])
            plotting.stream(iter,   np.array(f['initial_values']['bed']), 
                                    model_particles, 
                                    f['params']['x_max'][()], 
                                    y_limit, 
                                    fp_out)
        print(f'Creating gif of stream bed...')
        # Create a gif of the stream bed for iter_min to iter_max
        plotting.stream_gif(iter_min, iter_max, fp_out)
       
        print(f'Plotting flux and age plot...')
        # Plot downstream boundary crossing histogram
        plotting.crossing_info(f['final_metrics']['subregions'][ds_boundary_key][()], 
                                f['params']['n_iterations'][()], 1, fp_out)
        # Plot timeseries of particle age and downstream crossings
        plotting.crossing_info2(f['final_metrics']['subregions'][ds_boundary_key][()], 
                                np.array(f['final_metrics']['avg_age']),
                                f['params']['n_iterations'][()], 1, fp_out)

def parse_arguments():
    parser = argparse.ArgumentParser(description='Plotting script using the Python shelf files created by a py_BeRCM model run')
    parser.add_argument('path_to_file', help='Path to hdf5 file being plotted')
    parser.add_argument('save_location', help='Path to plot output location')
    parser.add_argument('iter_min', help='First iteration for stream plot')
    parser.add_argument('iter_max', help='Final iteration for stream plot')
    args = parser.parse_args()
    return args.path_to_file, args.save_location, int(args.iter_min), int(args.iter_max)

if __name__ == '__main__':
    # TODO: Replace save location param with the filename parsed of path and format
    path_to_file, save_location, iter_min, iter_max = parse_arguments()
    main(path_to_file, save_location, iter_min, iter_max)