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
import numpy as np
import warnings

from phonopy.structure.cells import get_primitive as get_primitive_phonopy
import spglib 
import ase, ase.data
from ase.data import atomic_numbers, chemical_symbols

import pymatgen.core.structure as str_pmg
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from pymatgen.io.vasp import Kpoints
from phonopy.structure.atoms import PhonopyAtoms

#def make_supercell(atoms0, P):
#    """ Make a supercell of the given structure
#    Args
#    ------
#    atoms0 : 
#        structure
#    
#    P : float, shape=(3,3)
#    """
#    return ase.build.make_supercell(atoms0, P)

def get_automatic_kmesh(struct_init, reciprocal_density=1500, 
        grid_density=0.01, method='reciprocal_density'):
    
    structure = change_structure_format(struct_init, format='pmg')
    
    if method == 'reciprocal_density':
        vol = structure.lattice.reciprocal_lattice.volume
        kppa = reciprocal_density * vol * structure.num_sites
        kpts = Kpoints.automatic_density(structure, kppa).kpts[0]
     
    #elif method == 'grid_density':
    #    Kpoints.automatic_density(structure, grid_density)
    
    else:
        warnings.warn(" Error: %s is not supported." % method)
        exit()
    
    return kpts

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
        warnings.warn(" WARRNING: the primitive could not be found.")

    return prim

def get_standardized_structure_spglib(
        struct_orig, to_primitive=False, format='ase'):
    """ Get a standardized cell shape with spglib and return its structure
    
    Args
    ------
    structure : ase, phonopy, pymatgen, ...
    
    version : string, 'new' or 'old'
        If 'old', spglib is used.
        If 'new', pymatgen is used.

    """
    structure = change_structure_format(struct_orig, format='phonopy')
    
    ### Get the standardized structure
    out = spglib.standardize_cell(structure, to_primitive=to_primitive)
    
    ### make the structure with the given format
    atoms = ase.Atoms(
            cell=out[0], pbc=True,
            scaled_positions=out[1],
            numbers=out[2],
            )
    
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
         
    elif version == 'pymatgen':
        
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


def change_structure_format(structure, format='pymatgen-IStructure'):
    """ Convert from arbitrary crystal format to an arbitrary crystal format

    Args
    -------
    structure (Structure object):

    Return
    -------
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
        warnings.warn(" Structure type {} is not supported".format(
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
        
        return PhonopyAtoms(
            cell=lattice,
            symbols=all_symbols,
            scaled_positions=coords,
            pbc=True
            )
    else:
        
        warnings.warn(" Structure type '{}' is not supported. "\
                "The structure type did not changed.".format(format))
        return structure

def get_formula(str_orig):    
    structure = change_structure_format(str_orig, format='pmg-istructure')
    return structure.composition.reduced_formula
    
def get_spg_number(str_orig):
    """ Regurn the international space group number
    """
    atoms = change_structure_format(str_orig, format='ase')
    cell = (
            atoms.cell,
            atoms.get_scaled_positions(),
            atoms.numbers)
    dataset = spglib.get_symmetry_dataset(cell)
    return dataset['number']

