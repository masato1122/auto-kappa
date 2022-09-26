import numpy as np
import warnings

from phonopy.structure.cells import get_primitive as get_primitive_phonopy
import spglib
import ase, ase.data
from ase.data import atomic_numbers

import pymatgen.core.structure as str_pmg
from phonopy.structure.atoms import PhonopyAtoms

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

def get_primitive_structure(structure, primitive_matrix=None, format='ase'):
    """ Return the primitive cell created by Phonopy. If primitive_matrix is not
    given, it will be suggested by SpacegroupAnalyzer in Pymatgen.
    """
    if primitive_matrix is None:
        primitive_matrix = get_primitive_matrix(structure)
    
    atoms_ph = get_primitive_phonopy(
            change_structure_format(structure, format='phonopy'),
            primitive_matrix)
    
    return change_structure_format(atoms_ph, format=format)

def get_primitive_matrix(structure):
    """ Return the primitive matrix of the given structure suggested by
    SpacegroupAnalyzer in Pymatgen. Different formats of structure are
    available. """
    from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
    str_pmg = change_structure_format(structure, format="pymatgen")
    spg_analyzer = SpacegroupAnalyzer(str_pmg)
    prim = spg_analyzer.find_primitive()
    pmat = spg_analyzer.get_conventional_to_primitive_transformation_matrix()
    return pmat

def get_suggested_primitive(structure, symprec=1e-5, format='ase'):
    """ Find and return the primitive cell
    """
    atoms_orig = change_structure_format(structure, format='ase')
    
    cell, scaled_positions, numbers = spglib.find_primitive(
            atoms_orig, symprec=symprec)
    
    disp = atoms_orig.get_scaled_positions()[0] - scaled_positions[0]

    scaled_positions += disp
    
    return _make_new_atoms(cell, scaled_positions, numbers)


def get_standerdized_cell(prim):
    """ Return the conventional cell created by Spglib
    """
    cell, scaled_positions, numbers = spglib.refine_cell(prim)
    return _make_new_atoms(cell, scaled_positions, numbers)


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
    import re
    structure = change_structure_format(str_orig, format='pmg-istructure')
    words = structure.get_primitive_structure().formula.split()
    if len(words) == 1:
        return re.sub(r"[123456789]", "", words[0])
    else:
        prefix = ""
        for i in range(len(words)):
            prefix += words[i].replace("1", "")
        return prefix
    return None

