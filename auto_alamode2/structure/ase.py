
from phonopy.structure.cells import get_primitive as get_primitive_ph
from .format import change_structure_format

def get_primitive(structure, pmatrix, format='ase'):
    atoms_ph = get_primitive_ph(
            change_structure_format(structure, format='phonopy'),
            pmatrix)
    return change_structure_format(atoms_ph, format=format)

