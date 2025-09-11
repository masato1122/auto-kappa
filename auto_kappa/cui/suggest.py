#
# suggest.py
#
# This file suggests the required parameters for the automation calculation.
#
# Copyright (c) 2023 Masato Ohnishi
#
# This file is distributed under the terms of the MIT license.
# Please see the file 'LICENCE.txt' in the root directory
# or http://opensource.org/licenses/mit-license.php for information.
#
import sys
import numpy as np

import ase
import spglib
from pymatgen.core.structure import Structure
from phonopy.structure.cells import get_supercell, get_primitive

from auto_kappa.structure import change_structure_format
from auto_kappa.structure.supercell import estimate_supercell_matrix

import logging
logger = logging.getLogger(__name__)

def suggest_structures_and_kmeshes(
    filename=None, structure=None, max_natoms=150, k_length=20, dim=3):
    """ Return the required parameters

    Parameters
    ----------
    filename : string
        structure filename

    structure : Structure obj
        structure object.
        ``filename`` or ``structure`` must be given.

    """
    ### check input parameters
    if filename is None and structure is None:
        msg = " Error: filename or structure must be given."
        logger.info(msg)
        sys.exit()
    
    ### prepare the structure obj.
    if filename is not None:
        struct_ase = change_structure_format(
                Structure.from_file(filename=filename),
                format='ase')
    else:
        struct_ase = change_structure_format(structure, format='ase')
    
    ### 2D: align the out-of-plane direction to the z-axis
    # if dim == 2:
    #     from auto_kappa.structure.two import set_out_of_plane_direction
    #     struct_ase = set_out_of_plane_direction(struct_ase)
        
    ### get the unitcell and the primitive matrix
    unitcell, prim_mat = get_unitcell_and_primitive_matrix(struct_ase)
    
    ### get the supercell matrix and supercell for FC2
    if dim == 3:
        sc_mat = estimate_supercell_matrix(unitcell, max_num_atoms=max_natoms)
    elif dim == 2:
        from auto_kappa.structure.two import estimate_supercell_matrix_2d
        sc_mat = estimate_supercell_matrix_2d(unitcell, max_num_atoms=max_natoms)
    else:
        logger.error("\n Error: dim must be 2 or 3.")
        sys.exit()
    
    supercell = change_structure_format(
            get_supercell(
                change_structure_format(unitcell, format='phonopy'), 
                sc_mat), 
            format='ase')
    
    ### get the primitive cell
    primitive = get_primitive(
                change_structure_format(unitcell, format='phonopy'),
                prim_mat,
                )
    try:
        primitive = get_primitive(
                change_structure_format(unitcell, format='phonopy'),
                prim_mat,
                )
        primitive = change_structure_format(primitive, format='ase')
    except Exception:
        msg = " Error: the primitive cell could not be obtained."
        logger.error(msg)
        sys.exit()

    ### collect obtained parameters
    structures = {"primitive": primitive,  "unitcell": unitcell,  "supercell": supercell}
    
    matrices = {"primitive": prim_mat, 
                "unitcell": np.identity(3).astype(int),
                "supercell": sc_mat}
    
    ### k-mesh
    kpts = {"primitive": klength2mesh(k_length, primitive.cell.array),
            "unitcell": klength2mesh(k_length, unitcell.cell.array),
            "supercell": klength2mesh(k_length, supercell.cell.array)}
    
    if dim == 2:
        for name in kpts:
            kpts[name][2] = 1
    
    return structures, matrices, kpts

def get_unitcell_and_primitive_matrix(structure):
    """ Get unitcell and primitive matrix of the structure written in the given
    file. The structure is standardized with Spglib.
    
    Parameters
    ----------
    structure : ASE Atoms obj.
        Crystal structure. Different cell shapes such as primitive, unit, and 
        supercell are accepted.
    
    Returns
    -------
    unitcell : Structure obj.
        conventional structure with ASE Atoms obj
    
    primitive_matrix : shape=(3,3), float
        Transformation matrix from the unitcell to the primitive cell.
        Row vectors
    
    """
    ### structure preparation
    cell = (structure.cell, structure.get_scaled_positions(), structure.numbers)
    
    ### get the unitcell
    ##
    ## Caution: spglib.standardize_cell() may change the out-of-plane 
    ## direction of 2D structures.
    ##
    cell_std = spglib.standardize_cell(cell, to_primitive=False)
    unitcell = ase.Atoms(
            cell=cell_std[0], pbc=True,
            scaled_positions=cell_std[1],
            numbers=cell_std[2])
    
    # ver.2
    # from pymatgen.symmetry.analyzer import SpacegroupAnalyzer as spg_analyzer
    # spg = spg_analyzer(change_structure_format(structure, 'pmg'))
    # conv = spg.get_conventional_standard_structure()
    # print(conv)
    # conv.to("POSCAR")
    # exit()
    
    ### primitive matrix
    cell_prim = spglib.standardize_cell(cell, to_primitive=True)
    primitive_matrix = np.dot(cell_prim[0], np.linalg.inv(cell_std[0])).T
    
    return unitcell, primitive_matrix

### Copy from phonopy/structure/grid_points.py in Phonopy
def klength2mesh(length, lattice, rotations=None):
    """Convert length to mesh for q-point sampling.
    
    This conversion for each reciprocal axis follows VASP convention by
    
    >>> N = max(1, int(l * |a|^* + 0.5))
    
    int means rounding down, not rounding to nearest integer.
    
    Parameters
    ----------
    length : float
        Length having the unit of direct space length.
    
    lattice : array_like
        Basis vectors of primitive cell in row vectors.
        dtype='double', shape=(3, 3)
    
    rotations : array of int, shape=(3,), optional
        Rotation matrices in real space. When given, mesh numbers that are
        symmetrically reasonable are returned. Default is None.

    Returns
    -------
    array_like : int, shape=(3,)
    
    Note
    -----
    This function is copied from Phonopy library.
    
    """
    rec_lattice = np.linalg.inv(lattice)
    rec_lat_lengths = np.sqrt(np.diagonal(np.dot(rec_lattice.T, rec_lattice)))
    mesh_numbers = np.rint(rec_lat_lengths * length).astype(int)

    if rotations is not None:
        from phonopy.structure.symmetry import get_lattice_vector_equivalence
        reclat_equiv = get_lattice_vector_equivalence(
            [r.T for r in np.array(rotations)])
        m = mesh_numbers
        mesh_equiv = [m[1] == m[2], m[2] == m[0], m[0] == m[1]]
        for i, pair in enumerate(([1, 2], [2, 0], [0, 1])):
            if reclat_equiv[i] and not mesh_equiv:
                m[pair] = max(m[pair])

    return np.maximum(mesh_numbers, [1, 1, 1])

