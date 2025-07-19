#
# utils.py
#
# Copyright (c) 2022 Masato Ohnishi
#
# This file is distributed under the terms of the MIT license.
# Please see the file 'LICENCE.txt' in the root directory
# or http://opensource.org/licenses/mit-license.php for information.
#
from ase.data import atomic_numbers
import ase
import pymatgen.core.structure as str_pmg
from phonopy.structure.atoms import PhonopyAtoms

import logging
logger = logging.getLogger(__name__)

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
        
        return PhonopyAtoms(
            cell=lattice,
            symbols=all_symbols,
            scaled_positions=coords,
            pbc=True
            )
    else:
        
        logger.warning(" Structure type '{}' is not supported. "\
                "The structure type did not changed.".format(format))
        return structure
