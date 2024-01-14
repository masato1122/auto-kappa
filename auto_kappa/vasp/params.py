# -*- coding: utf-8 -*-
#
# params.py
#
# Copyright (c) 2024 Masato Ohnishi
#
# This file is distributed under the terms of the MIT license.
# Please see the file 'LICENCE.txt' in the root directory
# or http://opensource.org/licenses/mit-license.php for information.
#
#
import sys
import math
import numpy as np
import glob

import logging
logger = logging.getLogger(__name__)

def get_amin_parameter(directory, lattice, **args):
    """ Get and return AMIN """
    ### get AMIN parameters
    from auto_kappa import default_amin_parameters
    amin_params = default_amin_parameters.copy()
    amin_params.update(args)
    
    ### get number of errors
    num_errors = get_number_of_errors(directory)
    if num_errors < amin_params['num_of_errors']:
        return None
    
    ### if error exists,
    amin = None
    for j in range(3):
        length = np.linalg.norm(lattice[j])
        if length > amin_params['tol_length']:
            return amin_params['value']
    return None

def get_number_of_errors(directory):
    """ Get and return the number of errors in the given directory. """
    num_errors = 0
    ### number of errors
    for suffix in ["tar", "tar.gz"]:
        line = directory + "/error.*." + suffix
        fns = glob.glob(line)
        num_errors += len(fns)

    ####
    #line = directory + "/INCAR"
    #fns = glob.glob(line)
    #num_errors += len(fns)
    return num_errors

