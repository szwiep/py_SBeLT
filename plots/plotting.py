"""
All things plotting related.
Most of this code created by Dr. Shawn Chartrand.
"""
import numpy as np
import matplotlib
from matplotlib import pyplot as plt
from matplotlib.collections import PatchCollection
from matplotlib.patches import Circle
from scipy.optimize import curve_fit
from scipy.special import factorial

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
    plt.savefig(plots_path)
        
    return

def flux_info(particle_flux_list, iterations, to_file):
    fig = plt.figure(10)
    ax = fig.add_subplot(1, 1, 1)
    bins = np.arange(-0.5, 11.5, 1) # fixed bin size
    plt.title('Histogram of Particle Flux, I = %i iterations' % n_iterations, fontsize=10, style='italic')
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
    if to_file:
        fig.savefig('../plots/FluxDownstreamBoundaryHist.png', format='png', dpi=600)
    else: 
        plt.show()

    #####
    Flux_CS = np.cumsum(particle_flux_list)
    Time = np.arange(1, iterations + 1, 1)
    fig = plt.figure(20)
    ax1 = fig.add_subplot(1,1,1)
    ax1.plot(Time, particle_flux_list, 'lightgray')
    plt.title('Timeseries of Particle Flux at Downstream Boundary')
    ax1.set_xlabel('Numerical Step')
    ax1.set_ylabel('Particle Flux')
    ax1.set_yticks([0, 1, 2, 3, 4, 5, 6, 7, 8])
    ax2 = ax1.twinx()

    ax2.plot(Time, Flux_CS, 'black')
    ax2.set_ylabel('Particle Flux Cumulative Sum', color='black', rotation=270, labelpad=15)
    ax2.tick_params('y', colors='black')
    
    fig.tight_layout()
    if to_file:
        fig.savefig('../plots/FluxDownstreamBoundary_2YAx.png', format='png', dpi=600)
    else:
        plt.show()
        
def flux_info2(particle_flux_list, particle_age_list, to_file):
    plt.clf()
    fig = plt.figure(20)
    ax3 = fig.add_subplot(1, 1, 1)
    #####
    Time = np.arange(1,  n_iterations + 1, 1)
    fig = plt.figure(20)
    ax3 = fig.add_subplot(1,1,1)
    ax3.plot(Time, particle_flux_list, 'lightgray')
    plt.title('Timeseries of Particle Flux at Downstream Boundary')
    ax3.set_xlabel('Numerical Step')
    ax3.set_ylabel('Particle Flux')
    ax3.set_yticks([0, 1, 2, 3, 4, 5, 6, 7, 8])
    ax4 = ax3.twinx()

    ax4.plot(Time, particle_age_list, 'black')
    ax4.set_ylabel('Particle Age (# of iterations)', color='black', rotation=270, labelpad=15)
    ax4.tick_params('y', colors='black')
    
    fig.tight_layout()
    if to_file:
        fig.savefig('../plots/FluxDownstreamBoundary_Age.png', format='png', dpi=600)
    else:
        plt.show()
        
def flux_info3(particle_flux_list, particle_age_list,particle_rage_list, to_file):
    plt.clf()
    fig = plt.figure(10)
    #####
    Time = np.arange(1,  n_iterations + 1, 1)
    fig = plt.figure(20)
    ax5 = fig.add_subplot(1,1,1)
    ax5.plot(Time, particle_flux_list, 'lightgray')
    plt.title('Timeseries of Particle Flux at Downstream Boundary')
    ax5.set_xlabel('Numerical Step')
    ax5.set_ylabel('Particle Flux')
    ax5.set_yticks([0, 1, 2, 3, 4, 5, 6, 7, 8])
    ax6 = ax5.twinx()

    plt.scatter(Time, particle_age_list, s = 15, c = particle_rage_list, cmap = 'RdGy' )
    ax6.set_ylabel('Mean Particle Age (# of iterations)', color='black', rotation=270, labelpad=15)
    ax6.tick_params('y', colors='black')
    plt.colorbar(orientation='horizontal',fraction=0.046, pad=0.2,label='Particle Age Range (max age - min age)')
    
    fig.tight_layout()
    if to_file:
        fig.savefig('../plots/FluxDownstreamBoundary_Rage.png', format='png', dpi=600)
    else:
        plt.show()