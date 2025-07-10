#
# two.py
#
# This script provides functions to treat 2D structures.
#
# Copyright (c) 2025 Masato Ohnishi
#
# This file is distributed under the terms of the MIT license.
# Please see the file 'LICENCE.txt' in the root directory
# or http://opensource.org/licenses/mit-license.php for information.
#
#
import numpy as np
from ase import Atoms
from ase.build import rotate

from auto_kappa.structure.crystal import change_structure_format

import logging
logger = logging.getLogger(__name__)

def get_out_of_plane_direction(structure) -> np.ndarray:
    """ Infer the out-of-plane direction from a structure.
    
    Parameters:
    - structure (Structure): pymatgen.Structure object
    - min_atoms (int): Minimum number of atoms to consider for out-of-plane direction inference.

    Returns:
    - array (float): Normal vector of the out-of-plane direction.
    """
    from sklearn.decomposition import PCA
    struct = structure.copy()
    struct = change_structure_format(struct, format='pmg-structure')
    
    # Make a supercell to ensure enough atoms for PCA
    struct.make_supercell([5, 5, 1])
    
    # Get atomic Cartesian coordinates
    coords = struct.cart_coords
    
    # Use PCA to find the out-of-plane direction (thinnest axis)
    pca = PCA(n_components=3)
    pca.fit(coords)
    
    normal_vec = pca.components_[-1]  # Get the last component (thinnest axis)
    
    if np.argmax(np.abs(normal_vec)) != 2:
        msg = "\n Warning: The out-of-plane direction is not the z-axis.\n"
        logger.warning(msg)
    
    return normal_vec

def align_to_z_axis(atoms: Atoms) -> Atoms:
    """ Rotate the atoms to align the out-of-plane direction with the z-axis.
    
    Parameters:
        atoms (ase.Atoms): Atoms object
    
    Returns:
        ase.Atoms: rotated Atoms object
    """
    # Get the out-of-plane direction
    normal = get_out_of_plane_direction(atoms)
    normal = normal / np.linalg.norm(normal)
    norm_idx = np.argmax(np.abs(normal))
    
    # Rotate the atoms to align the normal vector with the z-axis    
    atoms_rotated = atoms.copy()
    if norm_idx == 0:
        logger.info("\n Modify the out-of-plane direction from x-axis to the z-axis.")
        rotate(atoms_rotated, 'x', 'z', rotate_cell=True)
    elif norm_idx == 1:
        logger.info("\n Modify the out-of-plane direction from y-axis to the z-axis.")
        rotate(atoms_rotated, 'y', 'z', rotate_cell=True)
    elif norm_idx == 2:
        pass
    else:
        raise ValueError("Invalid normal vector index: {}".format(norm_idx))    
    return atoms_rotated

def estimate_supercell_matrix_2d(struct_orig, max_num_atoms=120):
    """ Estimate the supercell matrix for 2D structures.
    
    Parameters:
        structure : unitcell structure
        max_num_atoms (int): Maximum number of atoms in the supercell
        max_iter (int): Maximum number of iterations to find a suitable supercell
    
    Returns:
        np.ndarray: Supercell matrix
    """
    structure = change_structure_format(struct_orig, format='pmg-structure')
    
    # Estimate the supercell size based on maximum number of atoms
    n_atoms = len(structure)
    if n_atoms > max_num_atoms:
        sc_mat = np.diag([1, 1, 1])
    else:
        scale_factor = int(np.floor(np.sqrt(max_num_atoms / n_atoms)))
        sc_mat = np.diag([scale_factor, scale_factor, 1])
    return sc_mat

def get_thickness(structure, norm_idx=2):
    """ Calculate the thickness of a 2D structure.
    
    Parameters:
        structure : pymatgen.Structure object or ase.Atoms object
        out_idx (int): Index of the out-of-plane direction (default is 2 for z-axis)
    
    Returns:
        float: Thickness of the structure
    """
    from mendeleev import element
    struct = change_structure_format(structure, format='pmg-structure')
    chemical_symbols = [str(s) for s in struct.species]
    z_positions = struct.cart_coords[:, norm_idx]
    
    imax = np.argmax(z_positions)
    imin = np.argmin(z_positions)
    
    vdw_rad_max = element(chemical_symbols[imax]).vdw_radius * 0.01 # Convert to Angstrom
    vdw_rad_min = element(chemical_symbols[imin]).vdw_radius * 0.01 # Convert to Angstrom
    
    thickness = z_positions[imax] - z_positions[imin] + vdw_rad_max + vdw_rad_min
    
    return thickness
