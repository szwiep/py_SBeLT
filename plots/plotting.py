"""
All things plotting related.
Most of this code created by Dr. Shawn Chartrand.
"""
import numpy as np
import matplotlib
from matplotlib import pyplot as plt
from matplotlib.collections import PatchCollection
from matplotlib.collections import LineCollection
from matplotlib.patches import Circle
from scipy.optimize import curve_fit
from scipy.special import factorial
import seaborn as sns

def stream(iteration, bed_particles, model_particles, x_lim, y_lim,
                available_vertices, fp_out):
    """ Plot the complete stream from 0,0 to x_lim and y_lim. Bed particles 
    are plotted as light grey and model particles are dark blue. Allows
    for closer look at state of a subregion of the stream during simulation """
    plt.clf()
    fig = plt.figure(1)
    fig.set_size_inches(20, 6.5)
    ax = fig.add_subplot(1, 1, 1, aspect='equal')
    # NOTE: xlim and ylim modified for aspec ratio -- WIP
    ax.set_xlim((-2, x_lim))
    ax.set_ylim((0, y_lim))
    
    radius_array = np.asarray((bed_particles[:,1] / 2.0), dtype=float)
    x_center = bed_particles[:,0]
    y_center_bed = np.zeros(np.size(x_center))
    plt.rcParams['image.cmap'] = 'gray'
    ## This method of plotting circles comes from Stack Overflow questions\32444037
    ## Note that the patches won't be added to the axes, instead a collection will.
    patches = []
    for x1, y1, r in zip(x_center, y_center_bed, radius_array):
        circle = Circle((x1, y1), r)
        patches.append(circle)
    p = PatchCollection(patches, color="#BDBDBD", alpha=0.9, linewidths=(0, ))
    ax.add_collection(p)
    
    x_center_m = model_particles[:,0]
    y_center_m = model_particles[:,2]
    patches1 = []
    for x2, y2, r in zip(x_center_m, y_center_m, model_particles[:,1]/2):
        circle = Circle((x2, y2), r)
        patches1.append(circle)
    p_m = PatchCollection(patches1, cmap=matplotlib.cm.RdGy, edgecolors='black')
    p_m.set_array(model_particles[:,5])
    ax.add_collection(p_m)
    ### FOR TESTING: Plot various vertex types 
    # for xc in vertex_idx:
    #     plt.axvline(x=xc, color='b', linestyle='-')
    # for xxc in available_vertices:
    #     plt.axvline(x=xxc, color='b', linestyle='-', linewidth=0.25)
    # for green in chosen_vertex:
    #     plt.axvline(x=green, color='g', linestyle='-')
    ### 
    plt.colorbar(p_m,orientation='horizontal',fraction=0.046, pad=0.1,label='Particle Age (iterations since last hop)')
    plt.title(f'Iteration {iteration}')

    filename = f'iter{iteration}.png'
    plots_path = fp_out + filename
    plt.savefig(plots_path, format='png',)
        
    return

def flux_info(particle_flux_list, iterations, subsample, fp_out):
    plt.clf()
    fig = plt.figure(figsize=(8,7))
    ax = fig.add_subplot(1, 1, 1)
    bins = np.arange(-0.5, 11.5, 1) # fixed bin size
    plt.title('Histogram of Particle Flux, I = %i iterations' % iterations, fontsize=10, style='italic')
    plt.xlabel('Downstream Flux (particle count)')
    plt.ylabel('Fraction')
    ax.set_xlim((-1, max(bins)+1))
    ax.set_xticks([-1, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11])
    hist, bin_edges = np.histogram(particle_flux_list, bins=bins, density=True)

    # calculate binmiddles
    bin_middles = 0.5*(bin_edges[1:] + bin_edges[:-1])
    plt.bar(bin_middles,hist,color='lightgray')
    #fig.savefig('./ScriptTest/DownstreamFluxHistogram.pdf', format='pdf', dpi=600)

    # poisson function, parameter lamb is the fit parameter
    def poisson(k, lamb):
        return (lamb**k/factorial(k)) * np.exp(-lamb)

    # fit with curve_fit
    parameters, cov_matrix = curve_fit(poisson, bin_middles, hist)
    plt.plot(bin_middles, poisson(bin_middles, *parameters), color='black', marker='o', fillstyle = 'none', markersize=4, lw=0, markeredgecolor='black', markeredgewidth=1, label='Poisson PMF Fit')

    plt.legend(loc='upper right',frameon=0)
    filename = 'FluxDownstreamBoundaryHist.png'
    fi_path = fp_out + filename
    fig.savefig(fi_path, format='png', dpi=600)

    #####
    flux_list_avg = np.convolve(particle_flux_list, np.ones(subsample)/subsample, mode='valid')

    flux_list = flux_list_avg[0::subsample]
    Time = np.arange(1,  iterations + 1, subsample)
    Flux_CS = np.cumsum(flux_list)
    
    plt.clf()
    fig = plt.figure(figsize=(8,7))
    ax1 = fig.add_subplot(1,1,1)
    ax1.plot(Time, flux_list, 'lightgray')
    plt.title('Timeseries of Particle Flux at Downstream Boundary')
    ax1.set_xlabel('Numerical Step')
    ax1.set_ylabel('Particle Flux')
    ax1.set_yticks([0, 1, 2, 3, 4, 5, 6, 7, 8])
    ax2 = ax1.twinx()

    ax2.plot(Time, Flux_CS, 'black')
    ax2.set_ylabel('Particle Flux Cumulative Sum', color='black', rotation=270, labelpad=15)
    ax2.tick_params('y', colors='black')
    
    fig.tight_layout()
    filenameCS = 'FluxDownstreamBoundary_2YAx.png'
    fiCS_path = fp_out + filenameCS
    fig.savefig(fiCS_path, format='png', dpi=600)
        
def flux_info2(particle_flux_list, particle_age_list, n_iterations, subsample, fp_out):
    plt.clf()
    fig = plt.figure(figsize=(8,7))
    ax3 = fig.add_subplot(1, 1, 1)
    #####
    flux_list_avg = np.convolve(particle_flux_list, np.ones(subsample)/subsample, mode='valid')
    age_list_avg = np.convolve(particle_age_list, np.ones(subsample)/subsample, mode='valid')

    flux_list = flux_list_avg[0::subsample]
    age_list = age_list_avg[0::subsample]
    Time = np.arange(1,  n_iterations + 1, subsample)

    fig = plt.figure(figsize=(8,7))
    ax3 = fig.add_subplot(1,1,1)
    ax3.plot(Time, flux_list, 'lightgray')
    plt.title('Timeseries of Particle Flux at Downstream Boundary')
    ax3.set_xlabel('Numerical Step')
    ax3.set_ylabel('Particle Flux')
    ax3.set_yticks([0, 1, 2, 3, 4, 5, 6, 7, 8])
    ax4 = ax3.twinx()

    ax4.plot(Time, age_list, 'black')
    ax4.set_ylabel('Particle Age (# of iterations)', color='black', rotation=270, labelpad=15)
    ax4.tick_params('y', colors='black')
    
    fig.tight_layout()
    filename = 'FluxDownstreamBoundary_Age.png'
    fi_path = fp_out + filename
    fig.savefig(fi_path, format='png', dpi=600)

        
def flux_info3(particle_flux_list, particle_age_list,particle_rage_list, n_iterations, subsample, fp_out):
    plt.clf()
    fig = plt.figure(figsize=(8,7))
    ####
    flux_list_avg = np.convolve(particle_flux_list, np.ones(subsample)/subsample, mode='valid')
    age_list_avg = np.convolve(particle_age_list, np.ones(subsample)/subsample, mode='valid')
    age_range_list_avg = np.convolve(particle_rage_list, np.ones(subsample)/subsample, mode='valid')

    flux_list = flux_list_avg[0::subsample]
    age_list = age_list_avg[0::subsample]
    age_range_list = age_range_list_avg[0::subsample]

    Time = np.arange(1,  n_iterations + 1, subsample)

    # fig = plt.figure(figsize=(9,8))
    # ax5 = fig.add_subplot(1,1,1)
    # ax5.plot(Time, flux_list, 'lightgray')
    # plt.title('Timeseries of Particle Flux at Downstream Boundary')
    # ax5.set_xlabel('Numerical Step')
    # ax5.set_ylabel('Particle Flux')
    # ax5.set_yticks([0, 1, 2, 3, 4, 5, 6, 7, 8])
    # ax6 = ax5.twinx()

    # plt.plot(Time, age_list, linewidth='0.5', color='lightgray')
    # plt.scatter(Time, age_list, s = 15, c = age_range_list, cmap = 'RdGy')
    # ax6.set_ylabel('Mean Particle Age (# of iterations)', color='black', rotation=270, labelpad=15)
    # ax6.tick_params('y', colors='black')
    # plt.colorbar(orientation='horizontal',fraction=0.046, pad=0.2,label='Particle Age Range (max age - min age)')
    
    # fig.tight_layout()
    
    x = np.linspace(1,  n_iterations + 1, subsample)
    y = age_list
    dydx = age_range_list 
    flux_cumsum = np.cumsum(flux_list)

    # Create a set of line segments so that we can color them individually
    # This creates the points as a N x 1 x 2 array so that we can stack points
    # together easily to get the segments. The segments array for line collection
    # needs to be (numlines) x (points per line) x 2 (for x and y)
    points = np.array([x, y]).T.reshape(-1, 1, 2)
    segments = np.concatenate([points[:-1], points[1:]], axis=1)

    fig, axs = plt.subplots(figsize=(8,5))

    # Create a continuous norm to map from data points to colors
    norm = plt.Normalize(dydx.min(), dydx.max())
    lc = LineCollection(segments, cmap='RdGy', norm=norm)
    # Set the values used for colormapping
    lc.set_array(dydx)
    lc.set_linewidth(1)
    line = axs.add_collection(lc)

    plt.title('Timeseries of Particle Flux at Downstream Boundary')
    axs.autoscale(enable=True, axis='both')
    axsTwin = axs.twinx()
    axs.set_xlabel('Numerical Step')
    axs.set_ylabel('Particle Flux')
    axsTwin.set_ylabel('Mean Particle Age (# of iterations)', color='black', rotation=270, labelpad=15)
    fig.colorbar(line, ax=axs, pad=0.13)

    plt.plot(x, flux_cumsum, 'black')

    # fig.tight_layout()
    

    filename = 'FluxDownstreamBoundary_Rage.png'
    fi_path = fp_out + filename
    fig.savefig(fi_path, format='png', dpi=600)


def heat_map(shelf, n_iterations, window_subsample, fp_out):
     # Plot heatmap of all boundary crossings
        plt.clf()
        flux_list = []
        for i in range(shelf['param']['num_subregions']):
            key = f'subregion-{i}-flux' 
            flux_full = shelf[key]
            flux_list_avg = np.convolve(flux_full, np.ones(window_subsample), mode='valid') / window_subsample
            flux_list_ss = flux_list_avg[0::window_subsample]

            flux_list.append(flux_list_ss)

        # the content of labels of these yticks
        # print(flux_list)
        flux_list = np.transpose(flux_list)
        flux_heatmap = sns.heatmap(flux_list, cmap="coolwarm")
        fig = flux_heatmap.get_figure()

        
        filename = 'FluxAllBoundaries_Heatmap.png'
        fi_path = fp_out + filename
        fig.savefig(fi_path, format='png', dpi=600)