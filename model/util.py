import logging
import argparse
# import parameters
import math

def validate_parameters():
    """ A quick and dirty check of 
    parameters for crash-able violations.

    Raises:
    ------
    TypeError -- Invalid param type
    ValueError -- Invalid param value

    """ 
    logging.info('Validating parameters...')
    # if parameters.Pack <= 0:
    #     raise ValueError('Pack <=0. Pack fraction must be greater than 0')
    
    # if parameters.x_max <= 0:
    #     raise ValueError('x_max <=0. The stream must be longer than 0mm long')

    # if parameters.set_diam <= 0:
    #     raise ValueError('set_diam <= 0. A grain must have a diameter greater than 0')
    
    # if (math.remainder(parameters.x_max, parameters.num_subregions) != 0):
    #     raise ValueError('X_max must be divisible by the number of subregions')
        
    # if type(parameters.num_subregions) != int:
    #     raise TypeError('num_subregions needs to be Integer')
    
    # if parameters.n_iterations < 0:
    #     raise ValueError('n_iterations is negative. Cannot run negative number of iterations')

    # if parameters.sigma <= 0:
    #     raise ValueError('Sigma cannot be negative')
        
    # if type(parameters.n_iterations) != int:
    #     raise TypeError('n_iterations needs to be Integer')
    
    # if type(parameters.normal_dist) != bool:
    #     raise TypeError('normal_dist needs to be Boolean')
