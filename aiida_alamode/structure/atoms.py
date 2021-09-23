# -*- coding: utf-8 -*-
import numpy as np
import spglib
from pymatgen.core.structure import Structure as str_pm
from ase import Atoms as str_ase

def guess_primitive_matrix(structure, symprec=1e-5):
    """  Guess primitive matrix
    Args
    --------
    structure (pymatgen Structure or ase Atoms object):
        non-primitive structure
    
    Reference
    -----------
    phonopy.structure.cells.guess_primitive_matrix
    """
    if isinstance(structure, str_pm):
        cell = structure.lattice.matrix
        scaled_positions = structure.frac_coords
        numbers = structure.atomic_numbers
    elif isinstance(structure, str_ase):
        cell = structure.cell
        scaled_positions = structure.get_scaled_positions()
        numbers = structure.get_atomic_numbers()
    else:
        raise ValueError(' Structure type (%s) is not supported.'%(
            type(structure)))
    ##
    cell = (cell, scaled_positions, numbers)
    dataset = spglib.get_symmetry_dataset(cell, symprec=1e-5)
    tmat = dataset['transformation_matrix']
    centring = dataset['international'][0]
    
    ##
    from phonopy.structure.cells import get_primitive_matrix_by_centring
    pmat = get_primitive_matrix_by_centring(centring)
    return np.array(
            np.dot(np.linalg.inv(tmat), pmat), 
            dtype='double', order='C'
            )


