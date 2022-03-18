"""
All things plotting related.
Logic created by Dr. Shawn Chartrand.
"""
import numpy as np
import matplotlib
from matplotlib import pyplot as plt
from matplotlib.collections import PatchCollection
from matplotlib.patches import Circle
from scipy.optimize import curve_fit
from scipy.special import factorial

def stream(iteration, bed_particles, model_particles, x_lim, y_lim, out_location=None, out_name=None):
    """ Plot the complete stream from 0,0 to x_lim and y_lim. Bed particles 
    are plotted as light grey and model particles are plotted in a colour
    range dependant on their age. Allows for closer look at state of a subregion of the 
    stream during simulation. Also makes for fun gifs!


    Args:
        iteration: the iteration of the stream being plotted 
        bed_particles: array of all bed particles
        model_particles: array of all model particles
        x_lim: length of stream to plot
        y_lim: height of stream to plot
        out_location: save location
        out_name: filename (ignored if out_location not set) 
    """
    if out_location is not None:
        if out_name is None:
            raise ValueError('The out_name argument must be set if saving file.')
    fig = plt.figure(1)
    # fig.set_size_inches(20, 6.5)
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

    plt.colorbar(p_m,orientation='horizontal',fraction=0.046, pad=0.1,label='Particle Age (iterations since last hop)')
    plt.title(f'Iteration {iteration}')

    if out_location is None:
        plt.show()
    else:
        filename = f'iter{iteration}.png'
        plots_path = out_location + filename
        plt.savefig(plots_path, format='png',)
    return

def downstream_boundary_hist(particle_crossing_list, iterations, out_location=None, out_name=None):
    """Histogram of downstream particle crossings per iteration
    
    Args:
        particle_crossing_list: array of # of crossings per iteration
        iteration: number of iterations
        subsample: value to subsample by (use if data too large)
        out_location: save location. Default=None will print plot to screen
            but will not save.
        out_name: filename (ignored if out_location not set) 
    """
    if out_location is not None:
        if out_name is None:
            raise ValueError('The out_name argument must be set if saving file.')

    fig = plt.figure(figsize=(8,7))
    ax = fig.add_subplot(1, 1, 1)
    bins = np.arange(-0.5, 11.5, 1) # fixed bin size
    plt.title(f'Histogram of Particle Crossing, I = {iterations} iterations', fontsize=10, style='italic')
    plt.xlabel('Downstream Crossing (particle count)')
    plt.ylabel('Fraction')
    ax.set_xlim((-1, max(bins)+1))
    ax.set_xticks([-1, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11])
    hist, bin_edges = np.histogram(particle_crossing_list, bins=bins, density=True)

    # calculate binmiddles
    bin_middles = 0.5*(bin_edges[1:] + bin_edges[:-1])
    plt.bar(bin_middles,hist,color='lightgray')
    #fig.savefig('./ScriptTest/DownstreamCrossingHistogram.pdf', format='pdf', dpi=600)

    # poisson function, parameter lamb is the fit parameter
    def poisson(k, lamb):
        return (lamb**k/factorial(k)) * np.exp(-lamb)

    # fit with curve_fit
    parameters, cov_matrix = curve_fit(poisson, bin_middles, hist)
    plt.plot(bin_middles, poisson(bin_middles, *parameters), color='black', marker='o', fillstyle = 'none', markersize=4, lw=0, markeredgecolor='black', markeredgewidth=1, label='Poisson PMF Fit')

    plt.legend(loc='upper right',frameon=0)
    if out_location is None:
        plt.show()
    else:
        fi_path = out_location + out_name  + '.png'
        fig.savefig(fi_path, format='png', dpi=600)
    return

def downstream_boundary_ts(particle_crossing_list, iterations, subsample, out_location=None, out_name=None):
    """Time series plot of downstream particle crossings per iteration
    
    Args:
        particle_crossing_list: array of # of crossings per iteration
        iteration: number of iterations
        subsample: value to subsample by (use if data too large)
        out_location: save location
        out_name: filename (ignored if out_location not set) 
    """
    if out_location is not None:
        if out_name is None:
            raise ValueError('The out_name argument must be set if saving file.')
    crossing_list_avg = np.convolve(particle_crossing_list, np.ones(subsample)/subsample, mode='valid')

    crossing_list = crossing_list_avg[0::subsample]
    Time = np.arange(1,  iterations + 1, subsample)
    Crossing_CS = np.cumsum(crossing_list)
    
    fig = plt.figure(figsize=(8,7))
    ax1 = fig.add_subplot(1,1,1)
    ax1.plot(Time, crossing_list, 'lightgray')
    plt.title('Timeseries of Particle Crossing at Downstream Boundary')
    ax1.set_xlabel('Numerical Step')
    ax1.set_ylabel('Particle Crossing')
    ax1.set_yticks([0, 1, 2, 3, 4, 5, 6, 7, 8])
    ax2 = ax1.twinx()

    ax2.plot(Time, Crossing_CS, 'black')
    ax2.set_ylabel('Particle Crossing Cumulative Sum', color='black', rotation=270, labelpad=15)
    ax2.tick_params('y', colors='black')
    
    fig.tight_layout()
    if out_location is None:
        plt.show()
    else:
        fiCS_path = out_location + out_name  + '.png'
        fig.savefig(fiCS_path, format='png', dpi=600)
    return



def crossing_info_age(particle_crossing_list, particle_age_list, n_iterations, subsample, out_location=None, out_name=None):
    """Plot of particle crossing/flux vs. particle age. For very large values the
    data can be subsampled using the subsample parameter
    
    Args:
        particle_crossing_list: array of # of crossings per iteration
        particle_crossing_list: array of average particle age per iteration
        n_iteration: number of iterations
        subsample: value to subsample by (use if data too large)
        out_location: save location
        out_name: filename (ignored if out_location not set) 
    
    """
    if out_location is not None:
        if out_name is None:
            raise ValueError('The out_name argument must be set if saving file.')
    fig = plt.figure(figsize=(8,7))
    ax3 = fig.add_subplot(1, 1, 1)
    #####
    crossing_list_avg = np.convolve(particle_crossing_list, np.ones(subsample)/subsample, mode='valid')
    age_list_avg = np.convolve(particle_age_list, np.ones(subsample)/subsample, mode='valid')

    crossing_list = crossing_list_avg[0::subsample]
    age_list = age_list_avg[0::subsample]
    Time = np.arange(1,  n_iterations + 1, subsample)

    fig = plt.figure(figsize=(8,7))
    ax3 = fig.add_subplot(1,1,1)
    ax3.plot(Time, crossing_list, 'lightgray')
    plt.title('Timeseries of Particle Crossing at Downstream Boundary')
    ax3.set_xlabel('Numerical Step')
    ax3.set_ylabel('Particle Crossing')
    ax3.set_yticks([0, 1, 2, 3, 4, 5, 6, 7, 8])
    ax4 = ax3.twinx()

    ax4.plot(Time, age_list, 'black')
    ax4.set_ylabel('Particle Age (# of iterations)', color='black', rotation=270, labelpad=15)
    ax4.tick_params('y', colors='black')
    
    fig.tight_layout()
    if out_location is None:
        plt.show()
    else:
        fi_path = out_location + out_name + '.png'
        fig.savefig(fi_path, format='png', dpi=600)