#
# crystal.py
#
# This file helps to tread crystal structures and their reciprocal space.
#
# Copyright (c) 2022 Masato Ohnishi
#
# This file is distributed under the terms of the MIT license.
# Please see the file 'LICENCE.txt' in the root directory
# or http://opensource.org/licenses/mit-license.php for information.
#
import sys
import numpy as np
from scipy.optimize import linear_sum_assignment
import spglib 
import ase, ase.data

import pymatgen.core.structure as str_pmg
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from pymatgen.io.vasp import Kpoints

from auto_kappa.structure import change_structure_format
from auto_kappa.structure.two import get_normal_index

import logging
logger = logging.getLogger(__name__)

def get_automatic_kmesh(
    struct_init, 
    reciprocal_density=1500, 
    # grid_density=0.01,
    method='reciprocal_density', dim=3):
    
    structure = change_structure_format(struct_init, format='pmg')
    
    if method == 'reciprocal_density':    
        force_gamma = False
        vol = structure.lattice.reciprocal_lattice.volume           ### A^3
        kppa = reciprocal_density * vol * structure.num_sites       ### grid density
        kpts = Kpoints.automatic_density(structure, kppa, force_gamma=force_gamma).kpts[0]
        
        if dim == 2:            
            norm_idx_xyz = get_normal_index(structure, base='xyz')
            kpts = np.array(kpts)
            kpts[norm_idx_xyz] = 1
            kpts = kpts.tolist()
    
    #elif method == 'grid_density':
    #    Kpoints.automatic_density(structure, grid_density)
    
    else:
        msg = f"\n Method '{method}' is not supported. Please use 'reciprocal_density'."
        logger.error(msg)
        sys.exit()
    
    return list(kpts)

def get_commensurate_points(supercell_matrix):
    """ Get commensurate q-points.

    Args
    ------
    supercell_matrix : ndarray, int, shape=(3,3)
        supercell matrix with respect to the primitive cell

    Return
    -------
    q_pos : ndarray, float
        shape=(N, 3), N = det(supercell_matrix).
        Commensurate q-points

    """
    import ase
    from ase.build import make_supercell
    rec_prim = ase.Atoms(cell=np.identity(3), pbc=True)
    rec_prim.append(ase.Atom(position=np.zeros(3), symbol="H"))
    rec_scell = make_supercell(rec_prim, supercell_matrix.T)
    q_pos = rec_scell.get_scaled_positions()
    q_pos = np.where(q_pos > 1.-1e-15, q_pos-1., q_pos)
    return q_pos

#def get_primitive_structure(structure, primitive_matrix=None, format='ase'):
#    """ Return the primitive cell created by Phonopy. If primitive_matrix is not
#    given, it will be suggested by SpacegroupAnalyzer in Pymatgen.
#    """
#    if primitive_matrix is None:
#        primitive_matrix = get_primitive_matrix(structure)
#    
#    atoms_ph = get_primitive_phonopy(
#            change_structure_format(structure, format='phonopy'),
#            primitive_matrix)
#    
#    return change_structure_format(atoms_ph, format=format)

#def get_primitive_matrix(structure):
#    """ Return the primitive matrix of the given structure suggested by
#    SpacegroupAnalyzer in Pymatgen. Different formats of structure are
#    available. """
#    str_pmg = change_structure_format(structure, format="pymatgen")
#    spg_analyzer = SpacegroupAnalyzer(str_pmg)
#    pmat = spg_analyzer.get_conventional_to_primitive_transformation_matrix()
#    return pmat

def get_primitive_structure_spglib(structure, format='ase'):
    
    prim = get_standardized_structure(
            structure, to_primitive=True, format=format, version='spglib')
    
    if prim is None:
        logger.warning(" WARNING: the primitive could not be found.")

    return prim

def get_standardized_structure_spglib(struct_orig, to_primitive=False, format='ase'):
    """ Get a standardized cell shape with spglib and return its structure
    
    Args
    ------
    structure : ase, phonopy, pymatgen, ...
    
    version : string, 'new' or 'old'
        If 'old', spglib is used.
        If 'new', pymatgen is used.

    """
    structure = change_structure_format(struct_orig, format='phonopy')
    
    ## Both should provide the exactly same result.
    try:
        ## for new verison of spglib
        out = spglib.standardize_cell(structure, to_primitive=to_primitive)
    except Exception:
        ## PhonopyAtom => tuple
        cell = (
                structure.get_cell(),
                structure.get_scaled_positions(),
                structure.get_atomic_numbers()
                )
        out = spglib.standardize_cell(cell, to_primitive=to_primitive)
    
    ## make the structure with the given format
    atoms = ase.Atoms(cell=out[0], pbc=True, scaled_positions=out[1], numbers=out[2])
    return change_structure_format(atoms, format=format)

def get_standardized_structure(struct_orig, 
        to_primitive=False, format='ase', version='spglib'):
    """ Get a standardized cell shape and return its structure
    
    Args
    ------
    structure : ase, phonopy, pymatgen, ...
    
    version : string, 'new' or 'old'
        If 'old', spglib is used.
        If 'new', pymatgen is used.

    """
    if version == 'spglib':
        
        atoms = get_standardized_structure_spglib(
                struct_orig, 
                to_primitive=to_primitive,
                format=format,
                )
         
    elif version == 'pymatgen' or version == 'pmg':
        
        if (isinstance(struct_orig, str_pmg.Structure) == False and
                isinstance(struct_orig, str_pmg.IStructure) == False):
            structure = change_structure_format(struct_orig, format='pmg')
        else:
            structure = struct_orig
        
        spg_analyzer = SpacegroupAnalyzer(structure)
        
        if to_primitive:
            struct_stand = spg_analyzer.get_primitive_standard_structure()
        else:
            struct_stand = spg_analyzer.get_conventional_standard_structure()
        
        cell_stand = struct_stand.lattice.matrix    
        scaled_positions = struct_stand.frac_coords
        numbers = struct_stand.atomic_numbers
        
        ## 
        atoms = ase.Atoms(
                cell=cell_stand, pbc=True,
                scaled_positions=scaled_positions,
                numbers=numbers,
                )
    
    return change_structure_format(atoms, format=format)

def convert_primitive_to_unitcell(primitive, primitive_matrix, format='ase'):
    from phonopy.structure.cells import get_supercell
    mat_p2u = np.linalg.inv(primitive_matrix)
    mat_p2u = np.array(np.sign(mat_p2u) * 0.5 + mat_p2u, dtype="intc")
    unitcell = get_supercell(
        change_structure_format(primitive, format='phonopy'),
        mat_p2u)
    unitcell = change_structure_format(unitcell, format=format)
    return unitcell

def _make_new_atoms(cell, scaled_positions, numbers, pbc=True, center=False):
    """ Return an ase-Atoms object:

    Args
    -----
    cell : shape=(3,3)

    scaled_positions : shape=(natoms,3)
        scaled positions

    numbers : shape=(natoms)
        IDs of chemical symbol

    """
    # --- chemical symbols in the primitive cell
    symbols = []
    for ia, num in enumerate(numbers):
        symbols.append(ase.data.chemical_symbols[num])

    # --- make the primitive cell
    atoms_new = ase.Atoms(symbols, np.dot(scaled_positions, cell), cell=cell)
    if pbc:
        atoms_new.pbc = pbc
    if center:
        atoms_new.center()
    return atoms_new

def get_formula(str_orig):    
    structure = change_structure_format(str_orig, format='pmg-istructure')
    return structure.composition.reduced_formula
    
def get_spg_number(str_orig):
    """ Return the international space group number
    """
    atoms = change_structure_format(str_orig, format='ase')
    cell = (
            atoms.cell,
            atoms.get_scaled_positions(),
            atoms.numbers)
    dataset = spglib.get_symmetry_dataset(cell)
    try:
        return dataset.number
    except:    
        return dataset['number']

