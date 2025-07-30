# 
# compat.py
# 
# This script ...
# 
# Created on July 30, 2025
# Copyright (c) 2025 Masato Ohnishi
#
# This file is distributed under the terms of the MIT license.
# Please see the file 'LICENCE.txt' in the root directory
# or http://opensource.org/licenses/mit-license.php for information.
# 
import os
import numpy as np
import ase.io
import glob

from auto_kappa.structure import get_primitive_structure_spglib

import logging
logger = logging.getLogger(__name__)

def get_previously_used_structure(base_dir, tolerance=1e-5):
    """ Compare the optimized structure and the previously used structure
    
    Args
    ----
    base_dir (str): 
        Base directory containing the structure files.
    
    tolerance (float): 
        Tolerance for comparing cell parameters.
    
    Returns
    -------
    ase.Atoms or None:
        Returns the primitive cell of the previously used structure.
    
    """
    if base_dir.startswith('/'):
        base_dir = "./" + os.path.relpath(base_dir, os.getcwd())
    file_opt = base_dir + '/relax/structures/POSCAR.super'
    file_used = base_dir + '/harm/force/prist/POSCAR'
    
    if os.path.exists(file_used):
        structure_used = ase.io.read(file_used)
        prim_used = get_primitive_structure_spglib(structure_used)
        
        if os.path.exists(file_opt):
            structure_opt = ase.io.read(file_opt)
            prim_opt = get_primitive_structure_spglib(structure_opt)
        else:
            msg = "\n Read the previously used structure from %s" % file_used
            logger.info(msg)
            return prim_used
    else:
        if os.path.exists(file_opt):
            structure_opt = ase.io.read(file_opt)
            prim_opt = get_primitive_structure_spglib(structure_opt)
            msg = "\n Read the optimized structure from %s" % file_opt
            logger.info(msg)
            return prim_opt
        else:
            return None

    diff_cell = prim_opt.cell - prim_used.cell
    if np.any(np.abs(diff_cell) > tolerance):
        msg  = "\n Caution: the optimized structure and the previously used structure "
        msg += "\n          have different cell parameters."
        msg += "\n\n Lattice constant of the optimized structure:"
        msg += "\n (%s)" % file_opt
        msg += "\n\n %15.8f %15.8f %15.8f" % tuple(prim_opt.cell[0])
        msg += "\n %15.8f %15.8f %15.8f" % tuple(prim_opt.cell[1])
        msg += "\n %15.8f %15.8f %15.8f" % tuple(prim_opt.cell[2])
        msg += "\n\n Lattice constant of the previously used structure:"
        msg += "\n (%s)" % file_used
        msg += "\n\n %15.8f %15.8f %15.8f" % tuple(prim_used.cell[0])
        msg += "\n %15.8f %15.8f %15.8f" % tuple(prim_used.cell[1])
        msg += "\n %15.8f %15.8f %15.8f" % tuple(prim_used.cell[2])
        logger.warning(msg)
    
    msg = "\n Read the previously used structure from %s" % file_used
    logger.info(msg)    
    return prim_used
