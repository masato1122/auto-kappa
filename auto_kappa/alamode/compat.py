# 
# compat.py
# 
# Ensure compatibility across updated versions
# 
# Author      : M. Ohnishi
# Created on  : April 23, 2025
# 
# Copyright (c) 2025 Masato Ohnishi
#
# This file is distributed under the terms of the MIT license.
# Please see the file 'LICENCE.txt' in the root directory
# or http://opensource.org/licenses/mit-license.php for information.
#
import os
import sys
import glob
import subprocess
import shutil

import ase.io
# from ase.geometry import get_distances

from auto_kappa.io import AlmInput
from auto_kappa.structure import change_structure_format
# from auto_kappa.structure.two import get_normal_index
from auto_kappa.structure.comparison import match_structures

import logging
logger = logging.getLogger(__name__)

def _get_previously_suggested_structures(outdir):
    """ Read the previously suggested structures and their displacement patterns 
    to align with those from the previous version.
    """
    line = f"{outdir}/*/POSCAR"
    structures = {}
    for fn in glob.glob(line):
        try:
            key = fn.split("/")[-2]
            structures[key] = ase.io.read(fn)
        except:
            continue
    
    def _custom_sort_key(item):
        key = item[0]
        if key == 'prist':
            return (0, 0)
        try:
            numeric_key = int(key)
            return (1, numeric_key)
        except (ValueError, TypeError):
            return (2, str(key))
    
    return dict(sorted(structures.items(), key=_custom_sort_key))

def adjust_keys_of_suggested_structures(new_structures, outdir, tolerance=1e-5, dim=3):
    """ Sort the suggested structures and their displacement patterns 
    to align with those from the previous version.
    
    Args
    ------
    new_structures : dict
        The new structures
    
    outdir : str
        The output directory where the previous structures are stored.
    
    tolerance : float
        The tolerance for the distance between the new and previous structures.
    
    mag : float
        The magnitude of the atom displacement.
    """
    prev_structures = _get_previously_suggested_structures(outdir)
    
    ## Get index mapping from new to previous
    map_new2prev = {}
    assigned_keys = []
    for new_key, new_structure in new_structures.items():
        for prev_key, prev_structure in prev_structures.items():
            
            if prev_key in assigned_keys:
                continue
            
            # same = same_structures(
            #     new_structure, prev_structure,
            #     tolerance=tolerance, dim=dim)
            #
            same = match_structures(new_structure, prev_structure)
            
            if same:
                map_new2prev[new_key] = prev_key
                assigned_keys.append(prev_key)
                break
            
    ### Make a dict of structures with adjusted keys
    
    ## structures contained in prev_structures
    # avail_keys = [str(key) for key in list(new_structures.keys())]
    adjusted_key_structures = {}
    for new_key, prev_key in map_new2prev.items():
        adjusted_key_structures[prev_key] = new_structures[new_key]
        
    ## maximum key in prev_structures
    prev_key_max = 0
    for prev_key in prev_structures.keys():
        try:
            prev_key_max = max(prev_key_max, int(prev_key))
        except:
            pass
    
    ## new structures
    key_cur = prev_key_max + 1
    for new_key, new_structure in new_structures.items():
        if new_key not in map_new2prev:
            
            if new_key == 'prist':
                key = 'prist'
            else:
                key = str(key_cur)
                key_cur += 1
            
            adjusted_key_structures[key] = new_structure
    
    return adjusted_key_structures

# def min_cost_assignment(matrix):
#     from scipy.optimize import linear_sum_assignment
#     cost_matrix = np.array(matrix)
#     row_ind, col_ind = linear_sum_assignment(cost_matrix)
#     min_total_cost = cost_matrix[row_ind, col_ind].sum()
#     selected_elements = list(zip(row_ind, col_ind))
#     return selected_elements, min_total_cost

# def same_structures(struct1, struct2, tolerance=1e-7, dim=3):
#     """ Check whether the two structures are the same.
    
#     Args
#     ------
#     struct1 : ase.Atoms
#         The first structure.
    
#     struct2 : ase.Atoms
#         The second structure.
    
#     # cell : np.ndarray
#     #     The cell of the structure.
    
#     pristine : ase.Atoms
#         The pristine structure.
    
#     tolerance : float
#         The tolerance for the distance between the two structures.
    
#     mag : float
#         The magnitude of the atom displacement.
#     """
#     ## Check cell size
#     cell1 = struct1.cell.array
#     cell2 = struct2.cell.array
#     cell_diff = np.abs(cell1 - cell2)
    
#     pbc = [True, True, True]
#     if dim == 2:
#         norm_idx = get_normal_index(struct1)
#         cell_diff = np.delete(cell_diff, norm_idx, axis=1)
#         pbc[norm_idx] = False
#     else:
#         norm_idx = None
    
#     if np.amax(abs(cell_diff)) > tolerance:
#         # print("Cell size is different.")
#         return False
    
#     ## Check positions
#     pos1 = struct1.get_positions()
#     pos2 = struct2.get_positions()
#     D, D_len = get_distances(pos1, pos2, cell=struct1.cell, pbc=pbc)
    
#     ## Check the distance between atoms in the two structures
#     selected_elements, min_total_cost = min_cost_assignment(D_len)
    
#     n_atoms = len(struct1)
#     displacements = np.zeros((n_atoms, 3))
    
#     for i, j in selected_elements:
#         ## Check symbols of the atoms
#         if struct1[i].symbol != struct2[j].symbol:
#             print("symbol is different. (%s != %s)" % 
#                   (struct1[i].symbol, struct2[j].symbol))
#             return False
        
#         ## Check the displacement
#         disp = D[i, j]
#         displacements[i, :] = disp[:]
    
#     ## Check the displacement
#     ## Translational displacement is allowed
#     for j in range(3):
#         vmin = np.min(displacements[:, j])
#         vmax = np.max(displacements[:, j])
#         if vmax - vmin > tolerance:
#             # print("displacement is different.", j, vmax -vmin)
#             return False
    
#     return True
    
def check_directory_name_for_pristine(path_force, pristine):
    """ Check the directory name for the pristine structure.
    
    Args
    ------
    dir_force : str
        The directory name for the force calculation.
    
    pristine : ase.Atoms
        The pristine structure.
    """
    ## If "prist" directory already exists, return
    dir_prist = os.path.join(path_force, 'prist')
    if os.path.exists(dir_prist):
        return
    
    if os.path.exists(path_force) == False:
        return
    
    ## Get the directory names
    dirs_tmp = [entry.name for entry in os.scandir(path_force) if entry.is_dir()]
    
    labels = []
    for li in dirs_tmp:
        dir_i = os.path.join(path_force, li)
        fn = dir_i + "/POSCAR"
        if os.path.exists(fn):
            labels.append(li)
    
    ## Check the directory name for the pristine structure
    for lab in labels:
        fn = os.path.join(path_force, lab, "POSCAR")
        structure = ase.io.read(fn)
        # if same_structures(structure, pristine):
        if match_structures(structure, pristine):
            if lab != 'prist':
                dir1 = os.path.join(path_force, lab)
                msg = (
                    f"\n The directory name for the pristine structure "
                    f"was changed from \"{lab}\" to \"prist\".")
                logger.info(msg)
                os.rename(dir1, dir_prist)
                return 1
    
    return 0
    
    
def was_primitive_changed(struct_tmp, tol_prev, tol_new):
    """ Check whether the primitive cell was changed.
    
    Args
    ------
    structure : ase.Atoms
        The structure to be checked.
    
    tol_prev : float
        The tolerance for the previous version.
    
    tol_new : float
        The tolerance for the new version.
    """
    structure = change_structure_format(struct_tmp, format='pmg')
    
    prim_prev = structure.get_primitive_structure(tolerance=tol_prev)
    prim_new = structure.get_primitive_structure(tolerance=tol_new)
    
    if len(prim_prev) != len(prim_new):
        return True
    else:
        return False

def was_tolerance_changed(file_prev, new_params):
    """ Check whether the new parameters are different from the previous ones.
    
    Args
    ------
    file_prev : str
        The name of the previous ALAMODE input file.
    
    new_params : dict
        The new parameters to be compared with the previous ones.
    """
    if os.path.exists(file_prev) == False:
        return False
    
    ## read previous parameters
    prev_params = AlmInput().from_file(file_prev).as_dict()
    prev_tol = prev_params.get('tolerance')
    new_tol = new_params.get('tolerance')
    if prev_tol != new_tol:
        msg = f"\n Tolerance was changed from {prev_tol} to {new_tol}."
        logger.info(msg)
        return True
    else:
        return False

def backup_previous_results(directory, propt, prefix=None):
    """ Backup the previous results for ALAMODE in directory.
    """
    if os.path.exists(directory) == False:
        return
    
    ##
    if propt == 'suggest':
        cmd = f"rm {directory}/*"
        subprocess.run(cmd, shell=True)
    elif propt in ['fc2', 'fc3']:
        cmd = f"rm {directory}/{prefix}.* {directory}/{propt}.* {directory}/std_err.txt"
        subprocess.run(cmd, shell=True)
    elif propt in ['band', 'evec_commensurate', 'cv', 'kappa']:
        ##
        ## Make a backup directory and 
        ## all existing files are moved to the backup directory
        ##
        count = 1
        while True:
            out_backup = f"{directory}/backup{count}.tar.gz"
            if os.path.exists(out_backup) == False:
                break
            count += 1
        
        dir_backup = f"{directory}/backup{count}"
        
        ## Make a directory for backup
        os.makedirs(dir_backup, exist_ok=True)
        
        ## Move all files to the backup directory
        files = [f for f in os.listdir(directory) if os.path.isfile(os.path.join(directory, f))] 
        for f in files:
            if 'BORNINFO' in f:
                continue
            cmd = f"mv {directory}/{f} {directory}/backup{count}/"
            subprocess.run(cmd, shell=True)
        
        ## Compress the backup directory
        shutil.make_archive(
            base_name=dir_backup, 
            format='gztar', 
            root_dir=directory,
            base_dir=f"backup{count}")
        shutil.rmtree(dir_backup)
        
    else:
        print(directory, propt, prefix)
        cmd = f"\n Backup process for {propt} was not implemented yet."
        logger.info(cmd)
        sys.exit()
    

    # backup_dir = directory + "_backup"
    # if os.path.exists(backup_dir):
    #     os.system(f"rm -rf {backup_dir}")
    
    # os.system(f"mv {directory} {backup_dir}")

