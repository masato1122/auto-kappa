# 
# comparison.py
# 
# This script helps to compare crystal structures.
# 
# Created on July 30, 2025
# Copyright (c) 2025 Masato Ohnishi
#
# This file is distributed under the terms of the MIT license.
# Please see the file 'LICENCE.txt' in the root directory
# or http://opensource.org/licenses/mit-license.php for information.
# 
import numpy as np
import ase
from ase.geometry import get_distances
from pymatgen.core import Structure
from pymatgen.analysis.structure_matcher import StructureMatcher

from auto_kappa.structure import change_structure_format, get_primitive_structure_spglib

import logging
logger = logging.getLogger(__name__)

def match_structures(struct1, struct2, atol=1e-5, ltol=1e-5, stol=1e-5, angle_tol=0.1, 
                     ignore_order=False, primitive_cell=False, verbose=True):
    """ Check if two structures are the same using pymatgen's StructureMatcher.
    """
    if isinstance(struct1, Structure) == False:
        struct1 = change_structure_format(struct1, format='pmg-structure')
    if isinstance(struct2, Structure) == False:
        struct2 = change_structure_format(struct2, format='pmg-structure')
    
    ## Get primitive structures if requested
    if primitive_cell:
        struct1 = get_primitive_structure_spglib(struct1, format='pmg-structure')
        struct2 = get_primitive_structure_spglib(struct2, format='pmg-structure')
    
    ## Check number of atoms
    if len(struct1) != len(struct2):
        if verbose:
            logger.info("\n Number of atoms mismatch")
        return False
    
    if ignore_order:
        matcher = StructureMatcher(ltol=ltol, stol=stol, angle_tol=angle_tol, 
                                   scale=False, primitive_cell=primitive_cell)
        match = matcher.fit(struct1, struct2)
        return match
    else:
        ## Check cell size
        if not np.allclose(struct1.lattice.matrix, struct2.lattice.matrix, atol=atol):
            if verbose:
                logger.info("\n Lattice mismatch")
            return False
        
        ## Check atomic numbers
        if not np.allclose(struct1.atomic_numbers, struct2.atomic_numbers):
            if verbose:
                logger.info("\n Atomic number mismatch")
            return False
        
        ## Check atomic positions
        _, D_len = get_distances(struct1.cart_coords, struct2.cart_coords, 
                                 cell=struct1.lattice.matrix, pbc=True)
        disp_max = np.max(np.diag(D_len))
        if disp_max > atol:
            if verbose:
                logger.info("\n Cartesian coordinates mismatch")
            return False
        
        return True

def cells_equal(cell1, cell2, rtol=1e-5):
    """ Check if two cells are equal within a tolerance. """
    return np.allclose(cell1, cell2, rtol=rtol)

def atoms_equal(atoms1, atoms2, tol=1e-5, ignore_order=False):
    """ Check if two structures are equal within a tolerance. """
    
    if not isinstance(atoms1, ase.Atoms):
        atoms1 = change_structure_format(atoms1, format='ase')
    if not isinstance(atoms2, ase.Atoms):
        atoms2 = change_structure_format(atoms2, format='ase')
    
    if len(atoms1) != len(atoms2):
        return False
    
    if not np.allclose(atoms1.get_cell(), atoms2.get_cell(), atol=tol):
        return False

    if ignore_order:
        # sort atomic species and positions
        def get_sorted_data(atoms):
            symbols = atoms.get_chemical_symbols()
            positions = atoms.get_positions()
            data = list(zip(symbols, positions))
            data.sort(key=lambda x: (x[0], *x[1]))
            return data
        d1 = get_sorted_data(atoms1)
        d2 = get_sorted_data(atoms2)
    else:
        d1 = list(zip(atoms1.get_chemical_symbols(), atoms1.get_positions()))
        d2 = list(zip(atoms2.get_chemical_symbols(), atoms2.get_positions()))
    
    for (s1, p1), (s2, p2) in zip(d1, d2):
        if s1 != s2 or not np.allclose(p1, p2, atol=tol):
            return False
    
    return True
