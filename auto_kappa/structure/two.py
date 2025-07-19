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
import sys
import numpy as np
from pathlib import Path
import json

from ase import Atoms
from pymatgen.core import Structure, Lattice

from auto_kappa.structure import change_structure_format

import logging
logger = logging.getLogger(__name__)

def print_2d_system_notation():
    """ Print the notation for 2D systems. """
    msg =  "\n Note for 2D systems"
    msg += "\n ==================="
    msg += "\n 1. The out-of-plane direction is not always aligned with the z-axis "
    msg += "\n    due to the way VASP handles crystal symmetry."
    msg += "\n 2. Vacuum spacing is added along the out-of-plane direction."
    msg += "\n 3. The size of the vacuum spacing (t_vac) differs between VASP and ALAMODE calculations."
    msg += "\n    For ALAMODE, t_vac is set to be larger than the maximum atomic distance in the in-plane direction."
    msg += "\n    For VASP, t_vac is set to 20-30 Angstrom."
    logger.info(msg)

def print_length_info(structure):
    cutoff_harm = suggest_fc2_cutoff(structure)
    thickness = get_thickness(structure)
    diag_len_long = get_diagonal_length(structure, which='long')
    diag_len_short = get_diagonal_length(structure, which='short')
    msg  = f"\n - 2D material thickness   : {thickness:5.2f} Angstrom"
    msg += f"\n - Diagonal length (long)  : {diag_len_long:5.2f} Angstrom"
    msg += f"\n - Diagonal length (short) : {diag_len_short:5.2f} Angstrom"
    msg += f"\n => FC2 cutoff : {cutoff_harm:.2f} Angstrom"
    logger.info(msg)

def set_center(struct_orig) -> Structure:
    """ Set the position center of the structure to the center of the cell.
    
    Parameters:
        structure (Structure): pymatgen Structure object
    
    Returns:
        Structure: Centered structure
    """
    struct = change_structure_format(struct_orig, format='pmg-structure')
    coords = struct.cart_coords.copy()
    pos_center = np.mean(coords, axis=0)
    cell_center = np.mean(struct.lattice.matrix, axis=0)
    disp = cell_center - pos_center
    coords -= disp
    new_structure = Structure(
        lattice=struct.lattice,
        species=struct.species,
        coords=coords,
        coords_are_cartesian=True
    )
    return new_structure

def get_normal_index(struct_orig, base='abc') -> int:
    """ Estimate the vacuum direction by analyzing the real space gaps in a structure.
    
    Parameters:
        struct_orig: Structure object
        base (str): 'abc' for lattice vectors or 'xyz' for Cartesian axes
    
    Returns:
        int: Index of the direction with the largest vacuum gap
             For base='abc': (0 for a-axis, 1 for b-axis, 2 for c-axis)
             For base='xyz': (0 for x-axis, 1 for y-axis, 2 for z-axis)
    """
    structure = change_structure_format(struct_orig, format='pmg-structure')
    coords = np.array([site.coords for site in structure.sites])
    vacuum_gaps = []
    
    if base == 'abc':
        # Use lattice vectors
        cell = structure.lattice.matrix
        for axis, vec in enumerate(cell):
            unit_vec = vec / np.linalg.norm(vec)
            projections = coords @ unit_vec
            min_p, max_p = np.min(projections), np.max(projections)
            used = max_p - min_p
            total = np.linalg.norm(vec)
            gap = total - used
            vacuum_gaps.append(float(gap))
    
    elif base == 'xyz':
        # Use Cartesian axes
        for axis in range(3):
            coord_axis = coords[:, axis]
            min_p, max_p = np.min(coord_axis), np.max(coord_axis)
            used = max_p - min_p
            
            # Calculate total range in this Cartesian direction
            # Project lattice vectors onto this axis and find the span
            lattice_projections = []
            for i in range(2):
                for j in range(2):
                    for k in range(2):
                        point = (i * structure.lattice.matrix[0] + 
                                j * structure.lattice.matrix[1] + 
                                k * structure.lattice.matrix[2])
                        lattice_projections.append(point[axis])
            
            total = max(lattice_projections) - min(lattice_projections)
            gap = total - used
            vacuum_gaps.append(float(gap))
    
    else:
        msg = f" Error: Invalid base '{base}'. Use 'abc' or 'xyz'."
        logger.error(msg)
        sys.exit()
    
    norm_idx = int(np.argmax(vacuum_gaps))
    if base == 'abc':
        is_perpendicular(structure, norm_idx_abc=norm_idx)
    
    return norm_idx

# def get_normal_vector(structure, base='abc') -> np.ndarray:
#     norm_idx = get_normal_index(structure, base=base)
#     normal_vec = np.zeros(3)
#     normal_vec[norm_idx] = 1.0
#     return normal_vec

def is_perpendicular(struct_orig, norm_idx_abc, tol=1e-5) -> bool:
    """ Check if the angles of the structure are 90 degrees.
    
    Parameters:
        structure (Structure): pymatgen Structure object
        norm_idx (int): Index of the out-of-plane direction (default is None)
    
    Returns:
        bool: True if the angles are 90 degrees, False otherwise
    """
    structure = change_structure_format(struct_orig, format='pmg-structure')
    angles = structure.lattice.angles
    angles = np.delete(angles, norm_idx_abc)
    return all(abs(angle - 90.0) < tol for angle in angles)

def check_2d_structure(struct_orig: Atoms) -> bool:
    """ Check if the structure is a 2D structure.
    
    Parameters:
        structure (Atoms): Atoms object
    
    Returns:
        bool: True if the structure is 2D, False otherwise
    """
    if isinstance(struct_orig, Atoms) == False:
        atoms = change_structure_format(struct_orig, format='ase')
    else:
        atoms = struct_orig.copy()
    
    # if get_normal_index(atoms, base='xyz') != 2:
    #     msg = "\n Error: The out-of-plane direction is not aligned with the z-axis."
    #     logger.error(msg)
    #     sys.exit()
    
    norm_idx = get_normal_index(atoms, base='abc')
    if is_perpendicular(atoms, norm_idx):
        msg = "\n Error: The angles of the cell are not 90 degrees after switching axes."
        logger.error(msg)
        sys.exit()
    
    return 0

# def set_out_of_plane_direction(struct_orig) -> Atoms:
#     """ Rotate or/and switch axes of the structure to align the out-of-plane direction with the z-axis.
#     """
#     if isinstance(struct_orig, Atoms) == False:
#         atoms = change_structure_format(struct_orig, format='ase')
#     else:
#         atoms = struct_orig.copy()
    
#     norm_idx_abc = get_normal_index(atoms, base='abc')
#     norm_idx_xyz = get_normal_index(atoms, base='xyz')
    
#     if norm_idx_xyz == 2:
#         if norm_idx_abc == 2:
#             new_atoms = atoms.copy()
#         else:
#             new_atoms = switch_abc_axis(atoms, norm_idx_abc, norm_idx_xyz)
#             names = ['a', 'b', 'c']
#             msg = f"\n Switch {names[norm_idx_abc]}-axis to {names[norm_idx_xyz]}-axis."
#             logger.info(msg)
#     else:
#         x0 = [1, 0, 0]
#         y0 = [0, 1, 0]
#         z0 = [0, 0, 1]
#         new_atoms = atoms.copy()
#         if norm_idx_xyz == 0:
#             rotate(new_atoms, x0, z0, y0, y0) # x -> z
#             msg = "\n Rotate axis from x-axis to the z-axis."
#             logger.info(msg)
#         elif norm_idx_xyz == 1:
#             rotate(new_atoms, y0, z0, x0, x0) # y -> z
#             msg = "\n Rotate axis from y-axis to the z-axis."
#             logger.info(msg)
#         else:
#             msg = "\n Error: Invalid normal vector index."
#             logger.error(msg)
#             sys.exit()
        
#         norm_idx_abc = get_normal_index(new_atoms, base='abc')
#         if norm_idx_abc != 2:
#             new_atoms = switch_abc_axis(new_atoms, norm_idx_abc, 2)
    
#     ## Check if the out-of-plane direction is aligned with the z-axis
#     check_2d_structure(new_atoms)
    
#     return new_atoms
    
def switch_abc_axis(struct_orig: Atoms, idx1: int, idx2: int) -> Atoms:
    """ Switch the abc axes of the structure.
    """
    if isinstance(struct_orig, Atoms) == False:
        atoms = change_structure_format(struct_orig, format='ase')
    else:
        atoms = struct_orig
    
    new_cell = atoms.cell.copy()
    new_cell[idx1, :] = atoms.cell[idx2, :]
    new_cell[idx2, :] = atoms.cell[idx1, :]
    new_frac_coords = atoms.get_scaled_positions()
    new_frac_coords[:, idx1] = atoms.get_scaled_positions()[:, idx2]
    new_frac_coords[:, idx2] = atoms.get_scaled_positions()[:, idx1]
    
    new_atoms = Atoms(
        symbols=atoms.symbols, 
        scaled_positions=new_frac_coords,
        cell=new_cell, pbc=atoms.pbc)
    
    return new_atoms

# def align_to_z_axis(atoms: Atoms) -> Atoms:
#     """ Rotate the atoms to align the out-of-plane direction with the z-axis.
    
#     Parameters:
#         atoms (ase.Atoms): Atoms object
    
#     Returns:
#         ase.Atoms: rotated Atoms object
#     """
#     # Get the out-of-plane direction
#     norm_idx_abc = get_normal_index(atoms, base='abc')
#     norm_idx_xyz = get_normal_index(atoms, base='xyz')
#     if norm_idx_abc != norm_idx_xyz:
#         ## switch abc axis: norm_idx_abc -> norm_idx_xyz
#         atoms = switch_abc_axis(atoms, norm_idx_abc, norm_idx_xyz)
#         out = atoms.get_cell_lengths_and_angles()
#         if abs(out[3] - 90.) > 1e-5 or abs(out[4] - 90.) > 1e-5:
#             msg = "\n Error: The angles of the cell are not 90 degrees after switching axes."
#             logger.error(msg)
#             sys.exit()
        
#     # Rotate the atoms to align the normal vector with the z-axis
#     norm_idx_abc = get_normal_index(atoms, base='abc')
#     atoms_rotated = atoms.copy()
#     x0 = [1, 0, 0]
#     y0 = [0, 1, 0]
#     z0 = [0, 0, 1]
#     if norm_idx_xyz == 0:
#         logger.info("\n Modify the out-of-plane direction from x-axis to the z-axis.")
#         rotate(atoms_rotated, x0, z0, y0, y0, rotate_cell=True) # x -> z
#     elif norm_idx_xyz == 1:
#         logger.info("\n Modify the out-of-plane direction from y-axis to the z-axis.")
#         rotate(atoms_rotated, y0, z0, x0, x0, rotate_cell=True) # y -> z
#     elif norm_idx_xyz == 2:
#         pass
#     else:
#         raise ValueError("Invalid normal vector index: {}".format(norm_idx_xyz))
#     return atoms_rotated

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
    norm_idx = get_normal_index(structure, base='abc')
    
    # Estimate the supercell size based on maximum number of atoms
    n_atoms = len(structure)
    if n_atoms > max_num_atoms:
        sc_mat = np.diag([1, 1, 1])
    else:
        scale_factor = int(np.floor(np.sqrt(max_num_atoms / n_atoms)))
        sc_mat = np.eye(3) * scale_factor
        sc_mat[norm_idx, norm_idx] = 1
    
    return sc_mat

def get_thickness(structure, norm_idx=None) -> float:
    """ Calculate the thickness of a 2D structure.
    
    Parameters:
        structure : pymatgen.Structure object or ase.Atoms object
        out_idx (int): Index of the out-of-plane direction (default is 2 for z-axis)
    
    Returns:
        float: Thickness of the structure
    """
    struct = change_structure_format(structure, format='pmg-structure')
    chemical_symbols = [str(s) for s in struct.species]
    if norm_idx is None:
        norm_idx = get_normal_index(structure, base='xyz')
    
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

def set_vacuum_to_2d_structure(struct_orig, vacuum_thickness: float = 30.0) -> Structure:
    """ Add vacuum to a 2D structure by modifying the c-axis length. """
    
    structure = change_structure_format(struct_orig, format='pmg-structure')
    cell = structure.lattice.matrix
    norm_idx_abc = get_normal_index(structure, base='abc')
    norm_idx_xyz = get_normal_index(structure, base='xyz')
    
    ## Check the lattice vectors
    for j in range(3):
        if j != norm_idx_xyz and abs(cell[norm_idx_abc, j]) > 1e-5:
            msg = "\n Error: The c-axis of the structure is not aligned with the z-axis."
            logger.error(msg)
            sys.exit()
    
    ## Get the z-coordinates of the atoms
    out_coords = structure.cart_coords[:, norm_idx_xyz]
    out_center = np.mean(out_coords)
    material_thickness = get_thickness(structure, norm_idx=norm_idx_xyz)
    new_out_length = material_thickness + vacuum_thickness
    
    ## Calculate offset to center atoms in the z-direction
    offset = new_out_length / 2. - out_center
    
    ## Adjust atomic coordinates in the z-direction (cartesian coordinates)
    new_cart_coords = structure.cart_coords.copy()
    new_cart_coords[:, norm_idx_xyz] += offset
    
    ## Create new Structure
    new_cell = np.copy(cell)
    new_cell[norm_idx_abc, norm_idx_xyz] = new_out_length
    new_lattice = Lattice(new_cell)
    new_structure = Structure(
        lattice=new_lattice,
        species=structure.species,
        coords=new_cart_coords,
        coords_are_cartesian=True
    )
    return new_structure

def get_diagonal_length(struct_orig, norm_idx_abc=None, which='long') -> float:
    """ Calculate the longer diagonal length of a 2D structure. """
    
    if isinstance(struct_orig, Structure) == False:
        structure = change_structure_format(struct_orig, format='pmg-structure')
    else:
        structure = struct_orig.copy()
    
    if norm_idx_abc is None:
        norm_idx_abc = get_normal_index(structure, base='abc')
    
    cell = structure.lattice.matrix.copy()
    vecs = np.delete(cell, norm_idx_abc, axis=0)
    diag_vec1 = vecs[0] + vecs[1]
    diag_vec2 = vecs[0] - vecs[1]
    l1 = np.linalg.norm(diag_vec1)
    l2 = np.linalg.norm(diag_vec2)
    if which.startswith('long'):
        return max(l1, l2)
    elif which.startswith('short'):
        return min(l1, l2)
    else:
        msg = f" Error: Invalid 'which' parameter '{which}'. Use 'long' or 'short'."
        logger.error(msg)
        sys.exit()

def adjust_vacuum_size(orig_structure, scell_matrix=[[1,0,0], [0,1,0], [0,0,1]]) -> Structure:
    """ Modify the vacuum size of a 2D structure so that the cell size in the out-of-plane direction
    is larger than the in-plane dimensions.
    """
    structure = change_structure_format(orig_structure, format='pmg-structure')
    structure.make_supercell(scell_matrix)
    
    ## Get the current lattice parameters
    norm_idx_abc = get_normal_index(structure, base='abc')
    orig_length = np.linalg.norm(structure.lattice.matrix[norm_idx_abc, :])
    diag_length = get_diagonal_length(structure, norm_idx_abc=norm_idx_abc)
    mater_thickness = get_thickness(structure)
    max_length = diag_length + mater_thickness
    
    ## Calculate the new lattice parameters based on the supercell matrix
    new_structure = set_vacuum_to_2d_structure(orig_structure, vacuum_thickness=max_length * 1.2)
    
    msg = "\n Adjust the vacuum size to be larger than the maximum atomic distance (>= diagonal length + thickness):"
    msg += f"\n from {orig_length:.2f} to {max_length * 1.2:.2f} Angstrom."
    logger.info(msg)
    return new_structure

def suggest_fc2_cutoff(orig_sc, buffer=5.0) -> float:
    """ Suggest a cutoff for harmonic force constants based on the structure's diagonal length.
    
    Returns:
        float: Suggested cutoff value in Angstrom.
    """
    supercell = change_structure_format(orig_sc, format='pmg-structure')
    diag_length = get_diagonal_length(supercell)
    thickness = get_thickness(supercell)
    ## Suggest a cutoff based on the diagonal length
    cutoff_harm = diag_length + thickness + buffer
    return cutoff_harm
