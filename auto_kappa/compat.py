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
import sys
import numpy as np
import ase.io

# from auto_kappa.structure import change_structure_format
from auto_kappa.structure.crystal import (
    get_primitive_structure_spglib, 
    convert_primitive_to_unitcell
)
from auto_kappa.structure.comparison import match_structures, cells_equal, atoms_equal

import logging
logger = logging.getLogger(__name__)

def get_previously_used_structure(base_dir, prim_matrix):
    """ Compare the optimized structure and the previously used structure
    
    Args
    ----
    base_dir (str): 
        Base directory containing the structure files.
    prim_matrix (np.ndarray):
        Primitive matrix
    
    Returns
    -------
    ase.Atoms or None:
        Returns the primitive cell of the previously used structure.
    """
    base_dir = ("./" + os.path.relpath(base_dir, os.getcwd())
                if os.path.isabs(base_dir) else base_dir)
    
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
    match_wo_order = match_structures(struct_opt, struct_used, 
                                      ignore_order=True, verbose=True,
                                      ltol=1e-7, stol=1e-7)
    # match_with_order = match_structures(struct_opt, struct_used, 
    #                                     ignore_order=False, verbose=True,
    #                                     ltol=1e-7, stol=1e-7)
    
    if match_wo_order:
        ## Good! The optimized structure matches the previously used structure.
        # msg = "\n The optimized structure is read from %s." % files_opt[type]
        # logger.info(msg)
        # return struct_opt
        return struct_used
    else:
        cell_opt = struct_opt.cell
        cell_used = struct_used.cell
        msg = ""
        if cells_equal(cell_opt, cell_used) == False:
            msg = "\n *** Warning ***"
            msg += "\n The optimized structure and previously used structure have different "
            msg += "\n cell parameters."
            msg += "\n\n Supercell of the optimized structure:"
            msg += "\n (%s)" % files_opt['super']
            msg += "\n\n %15.8f %15.8f %15.8f" % tuple(cell_opt[0])
            msg += "\n %15.8f %15.8f %15.8f" % tuple(cell_opt[1])
            msg += "\n %15.8f %15.8f %15.8f" % tuple(cell_opt[2])
            msg += "\n\n Supercell of the previously used structure:"
            msg += "\n (%s)" % file_used
            msg += "\n\n %15.8f %15.8f %15.8f" % tuple(cell_used[0])
            msg += "\n %15.8f %15.8f %15.8f" % tuple(cell_used[1])
            msg += "\n %15.8f %15.8f %15.8f" % tuple(cell_used[2])
            msg += "\n\n The discrepancy may occur becuase the previous calculations used "
            msg += "\n optimized structure in an old calculation."
            msg += "\n If you want to use the optimized structure,"
            msg += "\n i.e. %s, " % files_opt['unit']
            msg += "\n please submit the job again with a different output directory name."
            
            cell_diff = cell_opt - cell_used
            max_cell_diff = np.max(abs(cell_diff.flatten()))
            if max_cell_diff > 0.1:
                msg += "\n\n Error: The maximum difference in the cell parameters is %.3f." % max_cell_diff
                
                base_dir = ("./" + os.path.relpath(base_dir, os.getcwd())
                          if os.path.isabs(base_dir) else base_dir)
                msg += "\n\n You may need to run a new job as follows:"
                msg += f"\n > mv {base_dir} {base_dir}-old"
                msg += f"\n > mkdir {base_dir}"
                msg += f"\n > cp -r {base_dir}-old/relax {base_dir}/"
                
                logger.error(msg)
                sys.exit()
        
        elif atoms_equal(struct_opt, struct_used, ignore_order=False):
            msg  = "\n Caution: the optimized structure and previously used structure do not match."
        
        msg += "\n\n The '%s' structure is obtained from %s." % (type, file_used)
        logger.warning(msg)
        return struct_used
