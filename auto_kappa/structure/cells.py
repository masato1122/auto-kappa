#
# cells.py
#
# Copyright (c) 2023 Masato Ohnishi
#
# This file is distributed under the terms of the MIT license.
# Please see the file 'LICENCE.txt' in the root directory
# or http://opensource.org/licenses/mit-license.php for information.
#
import numpy as np
import ase.io

def get_mat_u2p_spglib(atoms, format='ase'):
    """ Obtain the conventional to primitive transformation matrix of the given
    structure with Spglib and return the matrix.

    Args
    ------
    atoms : structure object

    format : string
        structure format of the return

    """
    import spglib
    from ase.data import chemical_symbols
    cell_orig = change_structure_format(atoms, format='phonopy')

    ### ver.1 w/o standardization
    #cell_prim = spglib.find_primitive(cell_orig)
    #cell_conv = spglib.refine_cell(cell_orig)
    #
    ### ver.2 with standardization
    cell_prim = spglib.standardize_cell(cell_orig, to_primitive=True)
    cell_conv = spglib.standardize_cell(cell_orig, to_primitive=False)

    if cell_prim is None or cell_conv is None:
        print("")
        print(" Error: Cannot obtain transformation matrix")
        print("")
        exit()

    ### transformation matrix: conventional to primitive
    mat_u2p = (cell_prim[0] @ np.linalg.inv(cell_conv[0])).T

    ### make an ASE Atoms obj.
    symbols = [chemical_symbols[num] for num in cell_conv[2]]
    unitcell = ase.Atoms(cell=cell_conv[0], pbc=True)
    isort = np.argsort(cell_conv[2])
    for ii in isort:
        unitcell.append(ase.Atom(
            position=cell_conv[1][ii],
            symbol=symbols[ii]
            ))
    return change_structure_format(unitcell, format=format), mat_u2p

