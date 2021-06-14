import numpy as np
import shelve
import argparse

import plotting


def main(filename, save_location, iter_min, iter_max):
    # Path to run-info file
    fp_in = filename
    # Path to output directory
    fp_out = save_location
    # Range to plot
    iteration_range = [iter_min, iter_max] 
    # Plot height
    y_limit = 10

    with shelve.open(fp_in, 'r') as shelf:
        for iter in range(iteration_range[0], iteration_range[1]+1):
            plotting.stream(iter, np.array(shelf["bed"]), np.array(shelf[str(iter)][0]), 
                            shelf["param"]["x_max"], y_limit, np.array(shelf[str(iter)][1]), fp_out)
        # Flux and age plots always plotted using information from ALL iterations
        plotting.flux_info(np.array(shelf["flux"]), shelf["param"]["n_iterations"], fp_out)
        plotting.flux_info2(np.array(shelf["flux"]), np.array(shelf["avg_age"]), shelf["param"]["n_iterations"], fp_out)
        plotting.flux_info3(np.array(shelf["flux"]), np.array(shelf["avg_age"]), np.array(shelf["age_range"]), 
                            shelf["param"]["n_iterations"], fp_out)

def parse_arguments():
    parser = argparse.ArgumentParser(description='Plotting script using the Python shelf files created by a py_BeRCM model run')
    parser.add_argument("filename", help="Path to run-info file to be plotted")
    parser.add_argument("save_location", help="Where to save plots")
    parser.add_argument("iter_min", help="First iteration to plot stream")
    parser.add_argument("iter_max", help="Final iteration to plot stream")
    args = parser.parse_args()
    return args.filename, args.save_location, int(args.iter_min), int(args.iter_max)


if __name__ == '__main__':
    filename, save_location, iter_min, iter_max = parse_arguments()
    main(filename, save_location, iter_min, iter_max)