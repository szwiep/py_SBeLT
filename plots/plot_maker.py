import numpy as np
import shelve
import argparse
import os

import plotting

from matplotlib import pyplot as plt


def main(filename, save_location, iter_min, iter_max):
    # Path to run-info file
    fp_in = filename
    # Path to output directory
    fp_out = save_location + '/'
    # Range to plot
    iteration_range = [iter_min, iter_max] 
    # Plot height
    y_limit = 10

    # Make appropriate sub-directory
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

    with shelve.open(fp_in, 'r') as shelf:
        subregion_max = shelf['param']['num_subregions']
        ds_boundary_key = f'subregion-{subregion_max-1}-flux' 

        # # # Plot streambed
        # for iter in range(iteration_range[0], iteration_range[1]+1):
        #     plotting.stream(iter, np.array(shelf['bed']), np.array(shelf[str(iter)][0]), 
        #                     shelf['param']['x_max'], y_limit, np.array(shelf[str(iter)][1]), fp_out)
       
        # # Flux and age information at Downstream boundary
        # plotting.flux_info(shelf[ds_boundary_key], shelf['param']['n_iterations'], 1000, fp_out)
        # plotting.flux_info2(shelf[ds_boundary_key], np.array(shelf['avg_age']), shelf['param']['n_iterations'], 1000, fp_out)
        # plotting.flux_info3(shelf[ds_boundary_key], np.array(shelf['avg_age']), np.array(shelf['age_range']), 
        #                     shelf['param']['n_iterations'], 1000, fp_out)
        
        plotting.heat_map(shelf, shelf['param']['n_iterations'], 1000, fp_out)
       


def parse_arguments():
    parser = argparse.ArgumentParser(description='Plotting script using the Python shelf files created by a py_BeRCM model run')
    parser.add_argument('filename', help='Path to run-info file to be plotted')
    parser.add_argument('save_location', help='Where to save plots')
    parser.add_argument('iter_min', help='First iteration to plot stream')
    parser.add_argument('iter_max', help='Final iteration to plot stream')
    args = parser.parse_args()
    return args.filename, args.save_location, int(args.iter_min), int(args.iter_max)


if __name__ == '__main__':
    filename, save_location, iter_min, iter_max = parse_arguments()
    main(filename, save_location, iter_min, iter_max)