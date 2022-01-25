# Copyright (c) aiida-alamode Team.
# Distributed under the terms of the MIT License.
# Written by M. Ohnishi

import numpy as np
import spglib
import warnings

import pymatgen.core.structure as str_pmg
from ase.data import chemical_symbols
from .format import change_structure_format

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
    
    print("")
    print(" WARNING: get_primitive_structure may not work properly!!!")
    print("")

    ## make a cell parameter
    cell = (cell, scaled_positions, numbers)
    
    ## find the primitive cell
    prim_cell = spglib.find_primitive(cell)
    
    ### compare with the result of this
    #from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
    #sga = SpacegroupAnalyzer(structure)
    #prim = sga.find_primitive()
    
    ## make pymatgen's IStructure
    prim_species = []
    for num in prim_cell[2]:
        prim_species.append(chemical_symbols[num])
    primitive = str_pmg.IStructure(prim_cell[0], prim_species, prim_cell[1])
    
    return primitive

