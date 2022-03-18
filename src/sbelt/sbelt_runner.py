""" 
This module is responsible for executing a run of the sbelt numerical model.
Default parameter/argument values are provided for the run function. See
the documentation in the project repository for more information on the 
arguments and how to change them. Have fun and entrain on!

Examples:
    The model execution can be initiated by executing the script with an 
    appropriate interpreter::
    
        $ python sbelt_runner.py

    or, by directly calling the run function once the module has been imported::

        sbelt_runner.run()

Attributes:
    ITERATION_HEADER: String used to delineate iterations in INFO-level logs
    ENTRAINMENT_HEADER: String used to idenitfy event particle being
                            entrained each iteration in INFO-level logs

Todo:
    * setuptools console_scripts and '$ python sbelt_runner.py' executions 
        can only use the default values currently. Need to add in sys.argv 
        parsing to let parameters be user-defined in these cases. For now
        users using sbelt from source code can pass desired arguments to
        the run call in the ``if __name_ == '__main__' function.``.
"""
import numpy as np
import h5py
import logging
from tqdm import tqdm

from sbelt import utils
from sbelt import logic

ITERATION_HEADER = ('Beginning iteration {iteration}...')
ENTRAINMENT_HEADER = ('Entraining particles {event_particles}')
logging.getLogger(__name__)

def run(iterations=1000, bed_length=100, particle_diam=0.5, particle_pack_dens = 0.78, \
                num_subregions=4, level_limit=3, poiss_lambda=5, gauss=False, gauss_mu=1, \
                gauss_sigma=0.25, data_save_interval=1, height_dependant_entr=False, \
                out_path='.', out_name='sbelt-out'): 
    """ Execute an sbelt run. 

    This function is responsible for calling appropriate logic
    and storing relevant information as outlined in project documentation. 
    It is the intended 
        
    An sbelt run consists of (generally) the following:
        (1) argument/parameter validation
        (2) building the stream data structures
        (3) running x iterations of entrainments over the stream
        (4) storing particle and flux data from each iteration
        

    Args:
        iterations: An int indicating number of iterations for model.
        bed_length: A float indicating the length of the stream to build 
        particle_diam: A float. 
        particle_pack_dens: A float indicating the packing fraction of 
            the model particles.
        num_subregions: An int representing the number of subregions to 
            define in the stream.
        level_limit: An int representing the maximum number of levels.
            permitted in-stream at any time (i.e how many particles high to stack)
        poiss_lambda: A float representing Lamba for poisson dist, 
            used to determine the number of entrainment events.
        gauss: A boolean flag indicating which distribution to sample from for 
            hop calculations. True=Normal, False=logNormal.
        gauss_mu: A float representing mean/expectation of the logNormal/Normal
            distribution for hop calculations.
        gauss_sigma: A float representing standard deviation of logNormal/Normal 
            distribution for hop calculations.
        data_save_interval: An int representing how often to record model 
            particle arrays (e.g 1=save every iteration, 2=save every other).
        height_dependant_entr: A boolean flag indicating whether model 
            automatically entrains particles that are on the level limit.
        out_path: A string representing the relative location to save the output.
        out_name: A string representing the name of the output file.
    """ 
    #############################################################################
    # validate parameters
    #############################################################################
    
    parameters = locals()
    utils.validate_arguments(parameters)

    #############################################################################
    #  Create model data and data structures
    # TODO: Better names for d, h variables. What would be more intuitive?
    #############################################################################

    print(f'Building Bed and Model particle arrays...')
    # Pre-compute d and h values for particle elevation placement
        # see d and h here: https://math.stackexchange.com/questions/2293201/
    d = np.divide(np.multiply(np.divide(particle_diam, 2), 
                                        particle_diam), 
                                        particle_diam)
    h = np.sqrt(np.square(particle_diam) - np.square(d))
    # Build the required structures for entrainment events
    bed_particles, model_particles, model_supp, subregions = build_stream(parameters, h)
    print(f'Bed and Model particles built.')

    #############################################################################
    #  Create entrainment data and data structures
    #############################################################################

    particle_age_array = np.ones(iterations)*(-1) # -1 represents an untouched element
    particle_range_array = np.ones(iterations)*(-1)
    snapshot_counter = 0

    #############################################################################

    h5py_filename = f'{out_name}.hdf5'
    hdf5_path = f'{out_path}/{h5py_filename}'

    with h5py.File(hdf5_path, "a") as f: 
        
        grp_p = f.create_group(f'params')
        for key, value in parameters.items():
            grp_p[key] = value

        grp_iv = f.create_group(f'initial_values')
        grp_iv.create_dataset('bed', data=bed_particles)
        grp_iv.create_dataset('model', data=model_particles)

        #############################################################################
        #  Entrainment iterations
        #############################################################################
        
        print(f'Model and event particle arrays will be written to {hdf5_path} every {data_save_interval} iteration(s).')
        print(f'Beginning entrainments...')
        for iteration in tqdm(range(iterations)):
            logging.info(ITERATION_HEADER.format(iteration=iteration))
            snapshot_counter += 1

            # Calculate number of entrainment events iteration
            e_events = np.random.poisson(parameters['poiss_lambda'], None)
            # Select n (= e_events) particles, per-subregion, to be entrained
            event_particle_ids = logic.get_event_particles(e_events, subregions,
                                                        model_particles, 
                                                        level_limit, 
                                                        height_dependant_entr)
            logging.info(ENTRAINMENT_HEADER.format(event_particles=event_particle_ids))
            # Determine hop distances of all event particles
            unverified_e = logic.compute_hops(event_particle_ids, model_particles, gauss_mu,
                                                    gauss_sigma, normal=gauss)
            # Compute available vertices based on current model_particles state
            avail_vertices = logic.compute_available_vertices(model_particles, 
                                                        bed_particles,
                                                        particle_diam,
                                                        level_limit,
                                                        lifted_particles=event_particle_ids)
            # Run entrainment event                    
            model_particles, model_supp, subregions = entrainment_event(model_particles, 
                                                                    model_supp,
                                                                    bed_particles, 
                                                                    event_particle_ids,
                                                                    avail_vertices, 
                                                                    unverified_e,
                                                                    subregions,
                                                                    iteration,  
                                                                    h)
            # Compute age range and average age, store in np arrays
            age_range = np.max(model_particles[:,5]) - np.min(model_particles[:,5])
            particle_range_array[iteration] = age_range

            avg_age = np.average(model_particles[:,5]) 
            particle_age_array[iteration] = avg_age

            # Record per-iteration information 
            if (snapshot_counter == data_save_interval):
                grp_i = f.create_group(f"iteration_{iteration}")
                grp_i.create_dataset("model", data=model_particles, compression="gzip")
                grp_i.create_dataset("event_ids", data=event_particle_ids, compression="gzip")
                snapshot_counter = 0

        #############################################################################
        # Store flux and age information
        #############################################################################
        
        print(f'Writting flux and age information to file...')
        grp_final = f.create_group(f'final_metrics')
        grp_sub = grp_final.create_group(f'subregions')
        for subregion in subregions:
            name = f'{subregion.getName()}-flux'
            flux_list = subregion.getFluxList()
            grp_sub.create_dataset(name, data=flux_list, compression="gzip")

        grp_final.create_dataset('avg_age', data=particle_age_array, compression="gzip")
        grp_final.create_dataset('age_range', data=particle_range_array, compression="gzip")
        print(f'Finished writing flux and age information.')

        print(f'Model run finished successfully.')
    return

#############################################################################
# Helper functions
#############################################################################

def build_stream(parameters, h):
    """ Build the data structures which define a stream.       

    Build array of m bed particles and array of n model particles. 
    Model particles will be assigned (x,y) positions which represent 
    resting on top of the bed particles. At the end of the build, each 
    model particle will have 2 support particles from the bed recorded 
    in an array. Finally, define an array of subregion objects. 

    Args:
        parameters: A dictionary of the 13 parameters required by the model.
            
        h: Geometric value used in calculations of particle placement. See
            in-line and project documentation for further explanation.

    Returns:
        bed_particles: An m-7 NumPy array representing the stream's m bed
            particles. For example::

            [[x, diam, y, uid, active, age, loops], ... ,[x, diam, y, uid, active, age, loops]]
        
        All bed particles should share the same diameter and y (elevation valu), all uids 
        must be negative, and to represent 'static' bed particles active = 0 and loops = 0. 

        model_particles: An n-7 NumPy array representing the stream's n model 
            particles in their initial placement. For example::
        
            [[x, diam, y, uid, active, age, loops], ... ,[x, diam, y, uid, active, age, loops]]

        model_supp: An n-2 NumPy array with the uids of the two particles supporting each 
            model particle. For example::
        
            model_supp = [[[-1,-2]], ... ,[-3,-4]]

        The above example states that model particle with uid 0 (model_supp[0]) is supported
        by bed particles with uids -1 and -2. Similarly, the model particle with uid n 
        (model_supp[n]) is supported by bed particles with uids -3 and -4.

        subregions: An array of Subregion objects
    
    """
    bed_particles = logic.build_streambed(parameters['bed_length'], parameters['particle_diam'])
    empty_model = np.empty((0, 7))      
    available_vertices = logic.compute_available_vertices(empty_model, bed_particles, parameters['particle_diam'],
                                                        parameters['level_limit'])    
    # Create model particle array and set on top of bed particles
    model_particles, model_supp = logic.set_model_particles(bed_particles, available_vertices, parameters['particle_diam'], 
                                                        parameters['particle_pack_dens'],  h)
    # Define stream's subregions
    subregions = logic.define_subregions(parameters['bed_length'], parameters['num_subregions'], parameters['iterations'])
    return bed_particles,model_particles, model_supp, subregions


def entrainment_event(model_particles, model_supp, bed_particles, event_particle_ids, avail_vertices, 
                                                                    unverified_e, subregions, iteration, h):
    """ This function mimics a single entrainment event through
    calls to the entrainment-related logic functions. 
    
    An entrainment event consists of:
        (1) Moving a set of the model particles 
            (event particles) downstream.
        (2) Recording crossings/flux across each downstream 
            boundary of the subregions.
        (3) Updating model particles states (i.e can it be
            selected for entrainment next iter?) and ages
    
    Args:
        model_particles: An n-7 NumPy array representing the stream's 
            n model particles. 
        model_supp: An n-2 NumPy array with the uids of the two particles supporting each 
            model particle (e.g model_supp[j] = supports for model particle j). 
        bed_particles: An m-7 NumPy array representing the stream's m 
            bed particles.
        event_particle_ids: A NumPy array of k uids representing the model particles
            that have been selected for entrainment.
        
    Returns:
        model_particles: Updated model_particles (Args) with updated age, location, 
            loops, and states, based on entrainment event placements.
        model_supp: An updated model_supports (Args) based on placements.
        subregions: Python array of Subregion objects with updated flux lists.
    """

    initial_x = model_particles[event_particle_ids][:,0]
    model_particles, model_supp = logic.move_model_particles(unverified_e, 
                                                                model_particles,
                                                                model_supp, 
                                                                bed_particles, 
                                                                avail_vertices,
                                                                h)
    final_x = model_particles[event_particle_ids][:,0]
    subregions = logic.update_flux(initial_x, final_x, iteration, subregions)
    model_particles = logic.update_particle_states(model_particles, model_supp)
    # Increment age at the end of each entrainment
    model_particles = logic.increment_age(model_particles, event_particle_ids)

    return model_particles, model_supp, subregions

if __name__ == '__main__':
    # assign argument values here if running from source code!
    run()
