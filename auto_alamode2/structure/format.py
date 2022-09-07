# Copyright (c) aiida-alamode Team.
# Distributed under the terms of the MIT License.
# Written by M. Ohnishi

import warnings

import pymatgen.core.structure as str_pmg
from ase import Atoms 
from phonopy.structure.atoms import PhonopyAtoms

from ase.data import atomic_numbers

def change_structure_format(structure, format='pymatgen-IStructure'):
    """ Convert to pytmagen's IStructure object 

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
    elif (isinstance(structure, Atoms) or
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
    if 'istructure' in form:
        return str_pmg.IStructure(lattice, all_symbols, coords)
    elif form == 'pymatgen-structure':
        return str_pmg.Structure(lattice, all_symbols, coords)
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
            symbols=all_symbols,
            scaled_positions=coords,
            pbc=True
            )
    else:
        warnings.warm(" Structure type '{}' is not supported. "\
                "The structure type did not changed.".format(format))
        return structure

