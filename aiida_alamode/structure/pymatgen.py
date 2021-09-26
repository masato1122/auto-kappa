# Copyright (c) ppdb Team.
# Distributed under the terms of the MIT License.
# Written by M. Ohnishi

import numpy as np
import spglib
import warnings

import pymatgen.core.structure as str_pm
from ase import Atoms 
from phonopy.structure.atoms import PhonopyAtoms

from ase.data import chemical_symbols

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
    if (isinstance(structure, str_pm.Structure) or
            isinstance(structure, str_pm.IStructure)):
        cell = structure.lattice.matrix
        scaled_positions = structure.frac_coords
        numbers = structure.atomic_numbers
        species = structure.species
    elif (isinstance(structure, Atoms) or 
            isinstance(structure, PhonopyAtoms)):
        cell = structure.cell
        scaled_positions = structure.get_scaled_positions()
        numbers = structure.get_atomic_numbers()
        species = structure.get_chemical_symbols()
    else:
        raise ValueError(" Structure type %s is not supported"%(
            type(structure)))
    ##
    cell = (cell, scaled_positions, numbers)

    ## find the primitive cell
    prim_cell = spglib.find_primitive(cell)
    
    ## make pymatgen's IStructure
    prim_species = []
    for num in prim_cell[2]:
        prim_species.append(chemical_symbols[num])
    primitive = str_pm.IStructure(prim_cell[0], prim_species, prim_cell[1])
    
    return primitive

def convert_to_pymatgen_structure(structure):
    """ Convert to pytmagen's IStructure object 

    Args
    -------
    structure (Structure object):
    
    Return
    -------
    pymatgen's IStructure object
    """
    if isinstance(structure, str_pm.IStructure):
        return structure
    elif isinstance(structure, str_pm.Structure):
        lattice = structure.lattice.matrix
        species = structure.species
        coords = structure.frac_coords
    elif (isinstance(structure, Atoms) or
            isinstance(structure, PhonopyAtoms)):
        lattice = structure.cell
        species = structure.get_chemical_symbols()
        coords = structure.get_scaled_positions()
    else:
        warnings.warn(" Structure type {} is not supported".format(
            type(structure)))
    
    return str_pm.IStructure(lattice, species, coords)

