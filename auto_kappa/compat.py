# 
# compat.py
# 
# This script helps to maintain compatibility with different versions of libraries or 
# repeated calculations.
# 
# Created on July 30, 2025
# Copyright (c) 2025 Masato Ohnishi
#
# This file is distributed under the terms of the MIT license.
# Please see the file 'LICENCE.txt' in the root directory
# or http://opensource.org/licenses/mit-license.php for information.
# 
import os
import ase.io
from pymatgen.core import Structure
from pymatgen.analysis.structure_matcher import StructureMatcher

from auto_kappa.structure import change_structure_format
from auto_kappa.structure.crystal import get_primitive_structure_spglib, convert_primitive_to_unitcell
from auto_kappa.structure.comparison import cells_equal, atoms_equal

import logging
logger = logging.getLogger(__name__)

def get_previously_used_structure(base_dir, prim_matrix, to_primitive=False, tolerance=1e-5):
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
    
    ## Optimized structures
    files_opt = {}
    structures_opt = {}
    for type in ['prim', 'unit', 'super']:
        files_opt[type] = base_dir + '/relax/structures/POSCAR.' + type
        if os.path.exists(files_opt[type]) == False:
            return None
        else:
            structures_opt[type] = ase.io.read(files_opt[type])
    
    ## Previously used structures
    file_used = base_dir + '/harm/force/prist/POSCAR'
    structures_used = {}
    if os.path.exists(file_used) == False:
        return None
    else:
        structures_used['super'] = ase.io.read(file_used)
    
    structures_used['prim'] = get_primitive_structure_spglib(structures_used['super'])
    structures_used['unit'] = convert_primitive_to_unitcell(structures_used['prim'], prim_matrix)
    
    ## Compare two structures using StructureMatcher
    type = 'unit'
    struct_opt = structures_opt[type]
    struct_used = structures_used[type]
    
    if isinstance(struct_opt, Structure) == False:
        struct_opt = change_structure_format(struct_opt, format='pmg-structure')
    if isinstance(struct_used, Structure) == False:
        struct_used = change_structure_format(struct_used, format='pmg-structure')

    matcher = StructureMatcher(
        ltol=1e-3, stol=1e-3, angle_tol=1.0,
        primitive_cell=to_primitive,
        scale=False)
    
    match = matcher.fit(struct_opt, struct_used)
    
    if match:
        ## Good! The optimized structure matches the previously used structure.
        msg = "\n The optimized structure is read from %s." % files_opt[type]
        logger.info(msg)
        return struct_opt
    else:
        if cells_equal(struct_opt.cell, struct_used.cell) == False:
            msg  = "\n Caution: the optimized structure and previously used structure "
            msg += "\n          have different cell parameters."
            msg += "\n\n Lattice constant of the optimized structure:"
            msg += "\n (%s)" % files_opt['super']
            msg += "\n\n %15.8f %15.8f %15.8f" % tuple(struct_opt.cell[0])
            msg += "\n %15.8f %15.8f %15.8f" % tuple(struct_opt.cell[1])
            msg += "\n %15.8f %15.8f %15.8f" % tuple(struct_opt.cell[2])
            msg += "\n\n Lattice constant of the previously used structure:"
            msg += "\n (%s)" % file_used
            msg += "\n\n %15.8f %15.8f %15.8f" % tuple(struct_used.cell[0])
            msg += "\n %15.8f %15.8f %15.8f" % tuple(struct_used.cell[1])
            msg += "\n %15.8f %15.8f %15.8f" % tuple(struct_used.cell[2])
        elif atoms_equal(struct_opt, struct_used, ignore_order=False):
            msg  = "\n Caution: the optimized structure and previously used structure do not match."
        
        msg += "\n\n The '%s' structure is obtained from %s." % (type, file_used)
        logger.warning(msg)
        return struct_used
