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
from auto_kappa.structure import change_structure_format

def cells_equal(cell1, cell2, tol=1e-5):
    """ Check if two cells are equal within a tolerance. """
    return np.allclose(cell1, cell2, atol=tol)

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


# def mirror_atoms(atoms, axis='x'):
#     """ Create a mirror-symmetric structure along a specified axis. """
#     mirrored = atoms.copy()
#     pos = mirrored.get_positions()
#     idx = {'x': 0, 'y': 1, 'z': 2}[axis]
#     pos[:, idx] *= -1
#     mirrored.set_positions(pos)
#     return mirrored

# def is_mirror(atoms1, atoms2, axis='x', tol=1e-3):
#     mirrored = mirror_atoms(atoms1, axis)
#     def sorted_data(atoms):
#         syms = atoms.get_chemical_symbols()
#         pos = atoms.get_positions()
#         data = list(zip(syms, pos))
#         data.sort(key=lambda x: (x[0], *x[1]))
#         return data
#     d1 = sorted_data(mirrored)
#     d2 = sorted_data(atoms2)
#     if len(d1) != len(d2):
#         return False
#     for (s1, p1), (s2, p2) in zip(d1, d2):
#         if s1 != s2 or not np.allclose(p1, p2, atol=tol):
#             return False
#     return True

# def check_mirror_symmetry(atoms1, atoms2, tol=1e-3):
#     results = {}
#     for axis in ('x', 'y', 'z'):
#         results[axis] = is_mirror(atoms1, atoms2, axis=axis, tol=tol)
#     return results
    
