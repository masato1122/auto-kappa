import numpy as np

from phonopy.structure.cells import get_primitive as get_primitive_ph
from .format import change_structure_format
import spglib
import ase, ase.data

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

def get_primitive(structure, pmatrix, format='ase'):
    """ Return the primitive cell created by Phonopy
    """
    atoms_ph = get_primitive_ph(
            change_structure_format(structure, format='phonopy'),
            pmatrix)
    return change_structure_format(atoms_ph, format=format)

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

