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

# from auto_kappa.structure import change_structure_format
from auto_kappa.structure.crystal import (
    get_primitive_structure_spglib, 
    convert_primitive_to_unitcell
)
from auto_kappa.structure.comparison import match_structures, cells_equal, atoms_equal

import logging
logger = logging.getLogger(__name__)


from auto_kappa.structure.crystal import get_primitive_structure_spglib, convert_primitive_to_unitcell
from auto_kappa.structure import change_structure_format
from phonopy.structure.cells import get_supercell

def make_structures(unitcell, prim_matrix, scell_matrix, format='ase'):
    from auto_kappa.structure import change_structure_format
    from phonopy import Phonopy
    phonon = Phonopy(
            change_structure_format(unitcell, format='phonopy'),
            scell_matrix,
            primitive_matrix=prim_matrix
            )
    unit = change_structure_format(phonon.unitcell , format=format) 
    prim = change_structure_format(phonon.primitive , format=format) 
    sc   = change_structure_format(phonon.supercell , format=format)
    return {"unit": unit, "prim": prim, "super": sc}

def make_structures2(unitcell, prim_matrix, scell_matrix, format='ase'):
    unit = change_structure_format(unitcell , format=format)
    prim = get_primitive_structure_spglib(unitcell)
    prim = change_structure_format(prim, format=format)
    sc = get_supercell(
            change_structure_format(unitcell, format='phonopy'),
            scell_matrix)
    sc = change_structure_format(sc, format=format)
    return {"unit": unit, "prim": prim, "super": sc}

def get_previously_used_structure(base_dir, prim_matrix, scell_matrix):
    """ Compare the optimized structure and the previously used structure
    
    Args
    ----
    base_dir (str): 
        Base directory containing the structure files.
    
    prim_matrix (np.ndarray):
        Primitive matrix
    
    scell_matrix (np.ndarray):
        Supercell matrix
    
    Returns
    -------
    ase.Atoms or None:
        Returns the primitive cell of the previously used structure.
    
    """
    if base_dir.startswith('/'):
        base_dir = "./" + os.path.relpath(base_dir, os.getcwd())
    
    # ## Check "./relax/freeze-1/CONTCAR" and "./relax/structures/POSCAR.unit"
    # file1 = base_dir + '/relax/freeze-1/CONTCAR'
    # file2 = base_dir + '/relax/structures/POSCAR.unit'
    # if os.path.exists(file1) and os.path.exists(file2):
    #     atoms1 = ase.io.read(file1, format='vasp')
    #     atoms2 = ase.io.read(file2, format='vasp')
    #     structures1 = make_structures(atoms1, prim_matrix, scell_matrix, format='ase')
    #     structures2 = make_structures(atoms2, prim_matrix, scell_matrix, format='ase')
    #     if match_structures(atoms1, atoms2):
    #         file_used_sc = base_dir + '/harm/force/prist/POSCAR'
    #         atoms_sc = ase.io.read(file_used_sc, format='vasp')
            
    #         out1 = match_structures(atoms_sc, structures1['super'], ignore_order=False)
    #         out2 = match_structures(atoms_sc, structures2['super'], ignore_order=False)
    #         out12 = match_structures(structures1['super'], structures2['super'], ignore_order=False)
    #         print(out1, out2, out12)
    
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
    type = 'super'
    struct_opt = structures_opt[type]
    struct_used = structures_used[type]
    match = match_structures(struct_opt, struct_used, verbose=False)
    
    if match:
        ## Good! The optimized structure matches the previously used structure.
        msg = "\n The optimized structure is read from %s." % files_opt[type]
        logger.info(msg)
        return struct_opt
    else:
        cell_opt = struct_opt.cell
        cell_used = struct_used.cell
        msg = ""
        if cells_equal(cell_opt, cell_used, tol=1e-5) == False:
            msg = "\n ## Caution ##"
            msg += "\n The optimized structure and previously used structure have "
            msg += "\n different cell parameters."
            msg += "\n\n Lattice constant of the optimized structure:"
            msg += "\n (%s)" % files_opt['super']
            msg += "\n\n %15.8f %15.8f %15.8f" % tuple(cell_opt[0])
            msg += "\n %15.8f %15.8f %15.8f" % tuple(cell_opt[1])
            msg += "\n %15.8f %15.8f %15.8f" % tuple(cell_opt[2])
            msg += "\n\n Lattice constant of the previously used structure:"
            msg += "\n (%s)" % file_used
            msg += "\n\n %15.8f %15.8f %15.8f" % tuple(cell_used[0])
            msg += "\n %15.8f %15.8f %15.8f" % tuple(cell_used[1])
            msg += "\n %15.8f %15.8f %15.8f" % tuple(cell_used[2])
        elif atoms_equal(struct_opt, struct_used, ignore_order=False):
            msg  = "\n Caution: the optimized structure and previously used structure do not match."
        
        msg += "\n\n The '%s' structure is obtained from %s." % (type, file_used)
        logger.warning(msg)
        return struct_used
