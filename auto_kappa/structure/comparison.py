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
from ase.geometry import get_distances, find_mic
from pymatgen.core import Structure
from pymatgen.analysis.structure_matcher import StructureMatcher

from auto_kappa.structure import change_structure_format, get_primitive_structure_spglib

import logging
logger = logging.getLogger(__name__)

def get_structure_matcher(atol=1e-5, ltol=1e-5, stol=1e-5, angle_tol=0.1, 
                          scale=False, primitive_cell=False):
    """ Get a StructureMatcher object with the specified tolerances. """
    return StructureMatcher(ltol=ltol, stol=stol, angle_tol=angle_tol, 
                            scale=scale, primitive_cell=primitive_cell)

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
        matcher = get_structure_matcher(
            atol=atol, ltol=ltol, stol=stol, angle_tol=angle_tol, 
            scale=False, primitive_cell=primitive_cell)
        match = matcher.fit(struct1, struct2)
        if match == False and verbose:
            logger.info("\n Atomic order mismatch")
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


# def generate_mapping_s2p(supercell: ase.Atoms, primitive: ase.Atoms, tol=1e-3):
#     """ Generate mapping from supercell to primitive cell.
    
#     Parameters:
#         supercell (ase.Atoms): The supercell structure.
#         primitive (ase.Atoms): The primitive cell structure.
#         tol (float): Tolerance for atomic position matching (in fractional coordinates).
    
#     Returns:
#         map_s2p (np.ndarray): Index mapping from supercell atoms to primitive cell atoms.
#         shift (np.ndarray): Lattice translation vectors (integers) applied to match atoms.
#     """
#     frac_super = supercell.get_scaled_positions()
#     frac_prim = primitive.get_scaled_positions()
    
#     nat_sc = len(supercell)
#     nat_p = len(primitive)
    
#     map_s2p = np.full(nat_sc, -1, dtype=int)
#     shift = np.zeros((nat_sc, 3), dtype=int)
    
#     # Find mapping
#     for i, pos_s in enumerate(frac_super):
#         found = False
#         for j, pos_p in enumerate(frac_prim):
#             # Compute minimum image difference
#             diff, _ = find_mic(pos_s - pos_p, cell=np.eye(3), pbc=True)
#             dist = np.linalg.norm(diff)
#             if dist < tol:
#                 map_s2p[i] = j
#                 shift[i] = np.round(pos_s - pos_p - diff).astype(int)
#                 found = True
#                 break
#         if not found:
#             raise RuntimeError(f" No matching atom found for supercell atom {i}")
    
#     return map_s2p, shift

def generate_mapping_s2p(supercell, primitive, tol_zero=1e-3):
    """ Generate mapping from supercell to primitive cell 
    """
    scell = supercell.cell.array
    pcell = primitive.cell.array
    frac_coords_super = supercell.get_scaled_positions(wrap=True)
    frac_coords_prim = primitive.get_scaled_positions(wrap=True)
    
    convertor = np.dot(scell, np.linalg.inv(pcell))
    convertor = np.round(convertor).astype(float)
    # print(convertor)
    
    nat_sc = len(supercell)
    shift = np.zeros((nat_sc, 3))
    map_s2p = np.zeros(nat_sc, dtype=int)
    
    for iat in range(nat_sc):
        xtmp = frac_coords_super[iat]
        xnew = np.dot(xtmp, convertor)
        
        iloc = -1
        diff_min = 1000.
        
        for jat in range(len(primitive)):
            xp = frac_coords_prim[jat]
            xdiff = (xnew - xp) % 1.0
            xdiff[xdiff >= 0.5] -= 1.0
            diff = np.linalg.norm(xdiff)
            diff_min = min(diff_min, diff)
            
            if diff < tol_zero:
                iloc = jat
                break
            
        if iloc == -1:
            msg = f"\n min. diff. : {diff_min:.6f}"\
                    "\n Warning : Equivalent atom not found. Check the relaxed structure."
            logger.error(msg)
            raise RuntimeError("Equivalent atom not found")
        
        map_s2p[iat] = iloc
        shift[iat] = np.round(xnew - frac_coords_prim[iloc])
        
    return map_s2p, shift
