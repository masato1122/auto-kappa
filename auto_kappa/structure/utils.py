#
# utils.py
#
# Copyright (c) 2022 Masato Ohnishi
#
# This file is distributed under the terms of the MIT license.
# Please see the file 'LICENCE.txt' in the root directory
# or http://opensource.org/licenses/mit-license.php for information.
#
import numpy as np
# from scipy.spatial.distance import cdist
from scipy.optimize import linear_sum_assignment

from ase.data import atomic_numbers, atomic_masses
import ase
from ase.geometry import get_distances
import pymatgen.core.structure as str_pmg
from phonopy.structure.atoms import PhonopyAtoms
from scipy.spatial import cKDTree

import logging
logger = logging.getLogger(__name__)

def change_structure_format(structure, format='pymatgen-IStructure'):
    """ Convert from arbitrary crystal format to an arbitrary crystal format

    Args
    -----
    structure (Structure object):

    Return
    ------
    pymatgen's IStructure object
    """
    if (isinstance(structure, str_pmg.Structure) or
        isinstance(structure, str_pmg.IStructure)):
        
        ## from pymatgen's (I)Structure object
        lattice = structure.lattice.matrix
        all_symbols = []
        for specie in structure.species:
            all_symbols.append(specie.name)
        coords = structure.frac_coords
    
    elif (isinstance(structure, ase.Atoms) or
          isinstance(structure, PhonopyAtoms)):
        
        ## from ase's Atoms and phonopy's PhonopyAtoms
        lattice = structure.cell
        all_symbols = structure.get_chemical_symbols()
        coords = structure.get_scaled_positions()
        
    else:
        logger.error(" Structure type {} is not supported".format(
            type(structure)))
        exit()
    
    ## set atomic numbers
    numbers = []
    for name in all_symbols:
        numbers.append(atomic_numbers[name])
    
    ## return the structure
    form = format.lower()
    if 'pymatgen' in form or 'pmg' in form:
        
        if format == 'pymatgen-structure' or format == 'pmg-structure':
            return str_pmg.Structure(lattice, all_symbols, coords)
        else:
            return str_pmg.IStructure(lattice, all_symbols, coords)
        
    elif form == 'ase' or form == 'atoms':
        
        return ase.Atoms(
            cell=lattice,
            scaled_positions=coords,
            numbers=numbers,
            pbc=True
            )
    
    elif form == 'phonopyatoms' or form == 'phonopy':
        
        masses = []
        for name in all_symbols:
            masses.append(atomic_masses[atomic_numbers[name]])
        
        return PhonopyAtoms(
            cell=lattice,
            symbols=all_symbols,
            masses=masses,
            scaled_positions=coords,
            pbc=True
            )
    else:
        
        logger.warning(" Structure type '{}' is not supported. "\
                "The structure type did not changed.".format(format))
        return structure

### Use pymatgen.analysis.structure_matcher.StructureMatcher
# def get_mapping_indices(mat1: np.ndarray, mat2: np.ndarray, cell=None, pbc=True,
#                         tolerance: float = 0.001, verbose=False) -> np.ndarray:
#     """ Get mapping indices from mat2 to mat1 using the Hungarian algorithm.
#     The function finds the optimal assignment that minimizes the total distance
#     between points in mat1 and mat2. If any matched pair exceeds the specified
#     tolerance, an error is raised.
    
#     Args
#     ------
#     mat1, mat2: np.ndarray
    
#     Returns
#     -------
#     mapping: np.ndarray
#         mapping[i] gives the index in mat1 that corresponds to mat2[i].
#         Example: if mapping[0] = 3, then mat2[0] corresponds to mat1[3].
#         mat2[i] --> mat1[mapping[i]]
        
#     How to use
#     -----------
#     >>> map_2to1 = get_mapping_indices(mat1, mat2, tolerance=0.001)
#     >>> for i2, i1 in enumerate(map_2to1):
#     >>>     print(f" mat2[{i2}] --> mat1[{i1}]")
#     """
#     assert mat1.shape == mat2.shape, "Input matrices must have the same shape."
    
#     ## calculate distance matrix between two sets of points
#     # dist_matrix = cdist(mat2, mat1)  # (n, n) matrix
#     _, dist_matrix = get_distances(mat2, mat1, cell=cell, pbc=pbc)
    
#     ## Solve the assignment problem to minimize total distance
#     row_ind, col_ind = linear_sum_assignment(dist_matrix)
    
#     ## Get the distances for the matched pairs
#     matched_distances = dist_matrix[row_ind, col_ind]
    
#     ## Check if any distance exceeds the tolerance
#     if np.any(matched_distances > tolerance):
#         if verbose:
#             max_dist = np.max(matched_distances)
#             msg = f"\n Error: Atom mapping failed. "
#             msg += f"Maximum distance {max_dist:.6f} exceeds tolerance {tolerance:.6f}."
#             msg += "\n Please check the input structures."
#             logger.warning(msg)
#         return None
    
#     mapping = np.zeros_like(row_ind)
#     mapping[row_ind] = col_ind
#     return mapping

# def get_translation_vector(struct_a, struct_b, tol=1e-3):
#     """ Find the translation vector that maps struct_b to struct_a.
#     The function assumes that both structures have the same atomic species
#     and number of atoms, but the order of atoms may differ. """
    
#     struct_a = change_structure_format(struct_a, format='ase')
#     struct_b = change_structure_format(struct_b, format='ase')
    
#     if len(struct_a) != len(struct_b):
#         raise ValueError(" The number of atoms in the two structures must be the same.")
    
#     coords_a = struct_a.get_positions()
#     coords_b = struct_b.get_positions()
    
#     ## Translate b's coordinates to match a's atomic positions
#     for target in coords_a:
#         delta = target - coords_b[0]  # Translation vector (Cartesian)
#         shifted_b = coords_b + delta
        
#         # struct_b_shifted = ase.Atoms(
#         #     cell=struct_b.cell,
#         #     positions=shifted_b,
#         #     symbols=struct_b.get_chemical_symbols(),
#         #     pbc=True
#         # )
        
#         mapping_b2a = get_mapping_indices(coords_a, shifted_b, cell=struct_a.cell, pbc=True, tolerance=tol)
#         if mapping_b2a is not None:
#             return delta
#     else:
#         raise ValueError("No valid translation vector found.")
    

# def get_translation_vector(struct_a, struct_b, tol=1e-3):
#     """ Find the translation vector that maps struct_b to struct_a.
#     The function assumes that both structures have the same atomic species
#     and number of atoms, but the order of atoms may differ. """
    
#     struct_a = change_structure_format(struct_a, format='pymatgen')
#     struct_b = change_structure_format(struct_b, format='pymatgen')
    
#     if len(struct_a) != len(struct_b):
#         raise ValueError(" The number of atoms in the two structures must be the same.")
    
#     frac_a = np.array([site.frac_coords for site in struct_a])
#     frac_b = np.array([site.frac_coords for site in struct_b])
    
#     ## Translate b's coordinates to match a's atomic positions
#     translation_candidates = []
    
#     # Try all atoms in a for the first atom in b
#     for target in frac_a:
#         delta = (target - frac_b[0]) % 1.0  # Translation vector (fractional)
#         shifted_b = (frac_b + delta) % 1.0
        
#         # Evaluate consistency of this translation
#         translation_candidates.append(delta)
#         tree = cKDTree(frac_a)
#         dists, _ = tree.query(shifted_b, k=1)
        
#         print(np.max(dists))
        
#         if np.all(dists < tol):
#             translation_frac = delta
#             break
#     else:
#         raise ValueError("No valid translation vector found.")
    
#     # Convert to Cartesian coordinates
#     translation_cart = struct_a.lattice.get_cartesian_coords(translation_frac)
#     return translation_cart

