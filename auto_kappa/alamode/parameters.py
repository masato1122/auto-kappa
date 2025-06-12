# -*- coding: utf-8 -*-
import math
import numpy as np
import pandas as pd

import logging
logger = logging.getLogger(__name__)

def set_parameters_evec(inp, primitive_matrix, scell_matrix):
    """ """
    ###### supercell matrix wrt primitive cell
    mat_p2s_tmp = np.dot(
            np.linalg.inv(primitive_matrix),
            scell_matrix
            )
    
    mat_p2s = np.rint(mat_p2s_tmp).astype(int)
    diff_max = np.amax(abs(mat_p2s - mat_p2s_tmp))
    if diff_max > 1e-3:
        msg = "\n CAUTION: please check the cell size of primitive "\
                "and supercell\n"
        msg += str(diff_max)
        logger.warning(msg)
    
    ### commensurate points
    from auto_kappa.structure.crystal import get_commensurate_points
    comm_pts = get_commensurate_points(mat_p2s)
    inp.update({'printevec': 1})
    inp.set_kpoint(kpoints=comm_pts)

def set_parameters_kappa(
        inp, kpts=None, nac=None, 
        isotope=2, kappa_coherent=1, 
        tmin=50, tmax=1000, dt=50,
        **kwargs):
    """ """
    params = {
            "kpts": kpts,
            "nac": nac,
            "isotope": isotope,
            "kappa_coherent": kappa_coherent,
            "tmin": tmin,
            "tmax": tmax,
            "dt": dt,
            }
    params.update(kwargs)
    inp.update(params)
