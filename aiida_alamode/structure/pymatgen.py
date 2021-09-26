# Copyright (c) ppdb Team.
# Distributed under the terms of the MIT License.
# Written by M. Ohnishi

import numpy as np
import spglib
import warnings

import pymatgen.core.structure as str_pm
from ase import Atoms 
from phonopy.structure.atoms import PhonopyAtoms

from ase.data import chemical_symbols, atomic_numbers

def get_primitive_structure(structure, symprec=1e-5):
    """  Get primitive structure
    Args
    --------
    structure (pymatgen's (I)Structure, ase's Atoms object, or PhonopyAtoms):
        non-primitive structure
    
    Return
    --------
    primitive : pymatgen's (I)Structure
    """
    ## set the structure format to be IStructure
    istr = change_structure_format(structure, format='IStructure')
    
    ## extract structure info
    cell = istr.lattice.matrix
    scaled_positions = istr.frac_coords
    numbers = istr.atomic_numbers
    species = istr.species
    
    ## make a cell parameter
    cell = (cell, scaled_positions, numbers)
    
    ## find the primitive cell
    prim_cell = spglib.find_primitive(cell)
    
    ## make pymatgen's IStructure
    prim_species = []
    for num in prim_cell[2]:
        prim_species.append(chemical_symbols[num])
    primitive = str_pm.IStructure(prim_cell[0], prim_species, prim_cell[1])
    
    return primitive

def change_structure_format(structure, format='pymatgen-IStructure'):
    """ Convert to pytmagen's IStructure object 

    Args
    -------
    structure (Structure object):
    
    Return
    -------
    pymatgen's IStructure object
    """
    if (isinstance(structure, str_pm.Structure) or
            isinstance(structure, str_pm.IStructure)):
        ## from pymatgen's (I)Structure object
        lattice = structure.lattice.matrix
        species = structure.species
        coords = structure.frac_coords
    elif (isinstance(structure, Atoms) or
            isinstance(structure, PhonopyAtoms)):
        ## from ase's Atoms and phonopy's PhonopyAtoms
        lattice = structure.cell
        species = structure.get_chemical_symbols()
        coords = structure.get_scaled_positions()
    else:
        warnings.warn(" Structure type {} is not supported".format(
            type(structure)))

    ## set atomic numbers
    numbers = []
    for specie in species:
        numbers.append(atomic_numbers[specie.name])
    
    ## return the structure
    form = format.lower()
    if 'istructure' in form:
        return str_pm.IStructure(lattice, species, coords)
    elif form == 'pymatgen-structure':
        return str_pm.Structure(lattice, species, coords)
    elif form == 'ase' or form == 'atoms':
        return Atoms(
            cell=lattice,
            scaled_positions=coords,
            numbers=numbers,
            pbc=True
            )
    elif form == 'phonopyatoms' or form == 'phonopy':
        return PhonopyAtoms(
            cell=lattice,
            symbols=species,
            scaled_positions=coords,
            pbc=True
            )
    else:
        warnings.warm(" Structure type '{}' is not supported. "\
                "The structure type did not changed.".format(format))
        return structure

