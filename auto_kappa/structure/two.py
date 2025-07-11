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
from pathlib import Path
import json

from ase import Atoms
from ase.build import rotate

from auto_kappa.structure.crystal import change_structure_format

import logging
logger = logging.getLogger(__name__)

def get_vacuum_direction_by_gap(struct_orig) -> int:
    """ Estimate the vacuum direction by analyzing the real space gaps in a structure.
    
    Returns:
        int: Index of the direction with the largest vacuum gap (0 for x, 1 for y, 2 for z).
    """
    structure = change_structure_format(struct_orig, format='pmg-structure')
    coords = np.array([site.coords for site in structure.sites])
    lattice = structure.lattice
    vacuum_gaps = []
    
    for axis, vec in enumerate(lattice.matrix):
        unit_vec = vec / np.linalg.norm(vec)
        projections = coords @ unit_vec
        min_p, max_p = np.min(projections), np.max(projections)
        used = max_p - min_p
        total = np.linalg.norm(vec)
        gap = total - used
        vacuum_gaps.append(float(gap))
    
    return int(np.argmax(vacuum_gaps))

def get_out_of_plane_direction(structure) -> np.ndarray:
    norm_idx = get_vacuum_direction_by_gap(structure)
    normal_vec = np.zeros(3)
    normal_vec[norm_idx] = 1.0
    return normal_vec
    
# def get_out_of_plane_direction_old(structure) -> np.ndarray:
#     """ Infer the out-of-plane direction from a structure.
    
#     Parameters:
#     - structure (Structure): pymatgen.Structure object
#     - min_atoms (int): Minimum number of atoms to consider for out-of-plane direction inference.

#     Returns:
#     - array (float): Normal vector of the out-of-plane direction.
#     """
#     from sklearn.decomposition import PCA
#     struct = structure.copy()
#     struct = change_structure_format(struct, format='pmg-structure')
    
#     # Make a supercell to ensure enough atoms for PCA
#     struct.make_supercell([5, 5, 1])
    
#     # Get atomic Cartesian coordinates
#     coords = struct.cart_coords
    
#     # Use PCA to find the out-of-plane direction (thinnest axis)
#     pca = PCA(n_components=3)
#     pca.fit(coords)
    
#     normal_vec = pca.components_[-1]  # Get the last component (thinnest axis)
    
#     if np.argmax(np.abs(normal_vec)) != 2:
#         msg = "\n Warning: The out-of-plane direction is not the z-axis.\n"
#         logger.warning(msg)
    
#     return normal_vec

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
    struct = change_structure_format(structure, format='pmg-structure')
    chemical_symbols = [str(s) for s in struct.species]
    
    ## Prepare van der Waals radii as a dictionary
    try:
        current_path = Path(__file__).parent
        file_vdw_radius = current_path / "elements_vdw_radii.json"
        with open(file_vdw_radius, 'r', encoding='utf-8') as f:
            data = json.load(f)
            vdw_radii_dict = {
                symbol: data[symbol]['vdw_radius_ang'] 
                for symbol in chemical_symbols if symbol in data}
    except FileNotFoundError:
        from mendeleev import element
        vdw_radii_dict = {}
        for symbol in chemical_symbols:
            if symbol not in vdw_radii_dict:
                try:
                    vdw_radii_dict[symbol] = element(symbol).vdw_radius * 0.01  # Convert to Angstrom
                except Exception as e:
                    logger.warning(f"Could not get vdw radius for {symbol}: {e}")
    
    z_positions = struct.cart_coords[:, norm_idx]
    imax = np.argmax(z_positions)
    imin = np.argmin(z_positions)
    vdw_rad_max = vdw_radii_dict.get(chemical_symbols[imax], 0)  # Default to 0 if not found
    vdw_rad_min = vdw_radii_dict.get(chemical_symbols[imin], 0)  # Default to 0 if not found
    
    thickness = z_positions[imax] - z_positions[imin] + vdw_rad_max + vdw_rad_min
    return thickness
