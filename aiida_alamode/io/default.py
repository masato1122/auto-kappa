# -*- coding: utf-8 -*-
import numpy as np
from .variables import *

from pymatgen import loadfn
from ..structure import guess_primitive_matrix

def get_alm_defaults(structure=None, filename=None, guess_primitive=False):
    """ Generate default variables for alm
    
    Args
    -------
    structure (pymatgen.core.structure.Structure object): structure
    filename (str): structure file name
        'structure' or 'filename' must be given.
    guess_primitive (bool):
        If it's True, the primitive structure will be guessed with a module
        in Phonopy.
     
    Returns
    --------
    default_variables (dictionary):
    """
    if structure is None:
        if filename is None:
            raise ValueError(" Structure or filename must be given.")
        else:
            structure = loadfn(filename)
    
    default_variables = alm_variables.copy()
    
    ## guess the primitive structure
    
    fn = 'POSCAR-unitcell'
    import ase.io
    atoms = ase.io.read(fn, format='vasp')
    if guess_primitive:
        prim = guess_primitive_matrix(atoms)
    
    exit()

    print(prim)
    prim.to(filename='POSCAR.out')
    
    
    exit()
    
    return default_variables


#def _set_cell(cells):


