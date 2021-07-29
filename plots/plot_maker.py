import numpy as np
import shelve
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

        # # Plot streambed
        for iter in range(iteration_range[0], iteration_range[1]+1):
            plotting.stream(iter, np.array(shelf['bed']), np.array(shelf[str(iter)][0]), 
                            shelf['param']['x_max'], y_limit, np.array(shelf[str(iter)][1]), fp_out)
       
        # Flux and age information at Downstream boundary
        plotting.flux_info(shelf[ds_boundary_key], shelf['param']['n_iterations'], 1000, fp_out)
        plotting.flux_info2(shelf[ds_boundary_key], np.array(shelf['avg_age']), shelf['param']['n_iterations'], 1000, fp_out)
        plotting.flux_info3(shelf[ds_boundary_key], np.array(shelf['avg_age']), np.array(shelf['age_range']), 
                            shelf['param']['n_iterations'], 1000, fp_out)
        
        plotting.heat_map(shelf, shelf['param']['n_iterations'], 1000, fp_out)

        ## Dynamic System Idea plotting below:
        # total_flux = np.zeros(1000000)
        # for subregion_id in range(0, subregion_max):
        #     subregion_key = f'subregion-{subregion_id}-flux' 
        #     total_flux = total_flux + shelf[subregion_key]
        #     print(f'Sub {subregion_id}: {shelf[subregion_key]}')
        # # Long term average
        # total_flux_avg = np.average(total_flux)
        # # Colvulve, subsample, get avg of subsample
        # total_flux_conv = np.convolve(total_flux, np.ones(1000)/1000, mode='valid')
        # ss_flux = total_flux_conv[0::1000]
        # ss_flux_avg = np.average(ss_flux)
        #  # Long term average - subsamples average
        # ss_total_flux_avg =  total_flux_avg - ss_flux_avg

        # # Convulve and subsample age as done to flux
        # avg_age = np.array(shelf['avg_age'])
        # avg_age_conv = np.convolve(avg_age, np.ones(1000)/1000, mode='valid')
        # ss_avg_age = avg_age_conv[0::1000]

        # # Subtract average from total flux
        # zero_mean_flux = ss_flux - ss_total_flux_avg

        # print(len(ss_avg_age), len(zero_mean_flux))
        # num_iters = np.arange(0, 1000)
        # plt.plot(zero_mean_flux, ss_avg_age,  alpha=0.7)
        # plt.scatter(zero_mean_flux, ss_avg_age, c = num_iters, cmap='inferno')
        # plt.colorbar(orientation='vertical',fraction=0.046,label='Iteration (Subsampled over window of 1000)')
        # plt.xlabel('Mean Particle Age')
        # plt.ylabel('Standardized Total Flux (Total Flux - Average Flux)')
        # plt.savefig('dynam_sysScatter6Sim35Line.png', format='png')
        # print(zero_mean_flux)


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