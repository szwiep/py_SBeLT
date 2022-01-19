import numpy as np
import h5py
import argparse
import os

from matplotlib import pyplot as plt
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
    try: 
        os.mkdir(save_location)
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
        # Plot streambed
        for iter in range(iteration_range[0], iteration_range[1]+1):
            plotting.stream(iter,   np.array(f['initial_values']['bed']), 
                                    np.array(f[f'iteration_{iter}']['model']), 
                                    f['params']['x_max'][()], 
                                    y_limit, 
                                    fp_out)
       
        # Plot downstream boundary crossing histogram
        plotting.crossing_info(f['final_metrics']['subregions'][ds_boundary_key][()], 
                                f['params']['n_iterations'][()], 1, fp_out)
        plotting.crossing_info2(f['final_metrics']['subregions'][ds_boundary_key][()], 
                                np.array(f['final_metrics']['avg_age']),
                                f['params']['n_iterations'][()], 1, fp_out)

def parse_arguments():
    parser = argparse.ArgumentParser(description='Plotting script using the Python shelf files created by a py_BeRCM model run')
    parser.add_argument('filename', help='Path to hdf5 file being plotted')
    parser.add_argument('save_location', help='Path to plot output location')
    parser.add_argument('iter_min', help='First iteration for streambed plot')
    parser.add_argument('iter_max', help='Final iteration for streambed plot')
    args = parser.parse_args()
    return args.filename, args.save_location, int(args.iter_min), int(args.iter_max)

if __name__ == '__main__':
    filename, save_location, iter_min, iter_max = parse_arguments()
    main(filename, save_location, iter_min, iter_max)