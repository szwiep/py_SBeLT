"""
A module for any utility-related functions for sbelt.
"""
import re 

def validate_arguments(parameters):
    """ Validate arguments for types/ranges.

    Args:
        parameters: A dictionary of the 13 parameters 
            required by the model. For example:

            {'particle_diam': 0.70,
            'bed_length': 100,
                ...
            'out_name': 'sbelt-out'} 
    
    Raises: 
        ValueError: if any argument is invalid. error 
            message provides further information.
    """
    # TODO: a lot of repeated code here - could be made prettier/simpler
    boolean_type_msg = "{failing_var} must be of type boolean (True/False)."
    boolean_type_vars = ['gauss', 'height_dependant_entr']
    for key in boolean_type_vars:
        if not isinstance(parameters[key], bool):
            raise ValueError(boolean_type_msg.format(failing_var=key))
    
    int_type_msg = "{failing_var} must be of type int."
    int_type_vars = ['bed_length', 'num_subregions', 'level_limit', 'iterations', 'data_save_interval']
    for key in int_type_vars:
        if not isinstance(parameters[key], int):
            raise ValueError(int_type_msg.format(failing_var=key))
    
    number_type_msg = "{failing_var} must be of type int or float."
    number_type_vars = ['particle_pack_dens', 'particle_diam', 'poiss_lambda', 'gauss_mu', 'gauss_sigma' ]
    for key in number_type_vars:
        if not isinstance(parameters[key], (int, float)):
            raise ValueError(number_type_msg.format(failing_var=key))

    string_type_msg = "{failing_var} must be of type string."
    string_type_vars = ['out_path', 'out_name']
    for key in string_type_vars:
        if not isinstance(parameters[key], str):
            raise ValueError(string_type_msg.format(failing_var=key))

    greater_than_0_msg = "{failing_var} must be > 0."
    greater_than_0_vars = ['bed_length','particle_pack_dens', 'particle_diam', 'num_subregions', 'level_limit', \
                                'iterations', 'gauss_sigma', 'data_save_interval']
    for key in greater_than_0_vars:
        if parameters[key] <= 0:
            raise ValueError(greater_than_0_msg.format(failing_var=key))
    
    geq_than_0_msg = "{failing_var} must be >= 0."
    geq_than_0_vars = ['poiss_lambda', 'gauss_mu']
    for key in geq_than_0_vars:
        if parameters[key] < 0:
            raise ValueError(geq_than_0_msg.format(failing_var=key))
    
    valid_filename_msg = "{failing_var} cannot contain spaces or invalid characters."
    valid_filename_vars = ['out_name']
    for key in valid_filename_vars:
        if not re.match(r"^[a-zA-Z\d\-\_]*$", parameters[key]):
            raise ValueError(valid_filename_msg.format(failing_var=key))

    # Check that particle_diam and num_subregions play well with bed_length!
    if parameters['bed_length'] % parameters['particle_diam'] != 0:
        bed_diam_msg = (
            "Invalid configuration of bed_length and particle_diam parameters: "
            " bed_length must be divisible by particle_diam."
        )
        print(bed_diam_msg)
        raise ValueError("bed_length must be divisible by set_diam")

    if parameters['bed_length'] % parameters['num_subregions'] != 0:
        bed_subr_msg = (
            "Invalid configuration of bed_length and num_subregions parameters: "
            " bed_length must be divisible by num_subregions."
        )
        print(bed_subr_msg)
        raise ValueError("bed_length must be divisible by num_subregions")

    return 