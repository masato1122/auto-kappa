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
import numpy as np
import pandas as pd
import glob
import subprocess
import shutil

import ase.io
from ase.geometry import get_distances

from auto_kappa.io import AlmInput
from auto_kappa.structure import change_structure_format

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
    return structures

def adjust_keys_of_suggested_structures(new_structures, outdir, tolerance=1e-3, mag=None):
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
    
    if len(new_structures) == len(prev_structures):
        return prev_structures
    
    if 'prist' in prev_structures.keys():
        struct_prist = prev_structures['prist']
    else:
        struct_prist = None
    
    ## Get index mapping from new to previous
    map_new2prev = {}
    for new_key, new_structure in new_structures.items():
        p1 = new_structure.get_positions()
        for prev_key, prev_structure in prev_structures.items():
            p2 = prev_structure.get_positions()
            
            same = same_structures(
                p1, p2,
                cell=new_structure.cell,
                pristine=struct_prist.get_positions() if struct_prist else None,
                mag=mag,
                tolerance=tolerance,
                pbc=new_structure.pbc)
            
            if same:
                map_new2prev[new_key] = prev_key
                break
            
            
    ### Make a dict of structures with adjusted keys
    
    ## structures contained in prev_structures
    avail_keys = [str(key) for key in list(new_structures.keys())]
    adjusted_key_structures = {}
    for new_key, prev_key in map_new2prev.items():
        # print(f" new {new_key} -> prev {prev_key}")
        adjusted_key_structures[prev_key] = new_structures[new_key]
        avail_keys.remove(str(prev_key))
    
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
            
    # if 'cube' in outdir:
    #     ## remove the structures that are not in the previous version
    #     print(prev_structures.keys())
    #     print(adjusted_key_structures.keys())
    #     sys.exit()
    
    return adjusted_key_structures

def same_structures(
    p1, p2, cell=None, pristine=None, pbc=[True, True, True],
    tolerance=1e-3, mag=None):
    """ Check whether the two structures are the same.
    
    Args
    ------
    p1 : np.ndarray
        The positions of the first structure.
    
    p2 : np.ndarray
        The positions of the second structure.
    
    cell : np.ndarray
        The cell of the structure.
    
    pristine : ase.Atoms
        The pristine structure.
    
    tolerance : float
        The tolerance for the distance between the two structures.
    
    mag : float
        The magnitude of the atom displacement.
    """
    D, D_len = get_distances(p1, p2, cell=cell, pbc=pbc)
    dists = np.diagonal(D_len)
    if np.all(dists < tolerance):
        return True
    
    return False
    
    # if pristine is None:
    #     return False
    
    # ### Compare with the pristine structure
    # D1, D_len1 = get_distances(p1, pristine, cell=cell, pbc=pbc)
    # dists1 = np.diagonal(D_len1)
    # idx1 = np.where(abs(dists1) > tolerance)[0]
    
    # D2, D_len2 = get_distances(p2, pristine, cell=cell, pbc=pbc)
    # dists2 = np.diagonal(D_len2)
    # idx2 = np.where(abs(dists2) > tolerance)[0]
    
    # if len(idx1) == len(idx2) and len(idx1) == 1:
    #     iat1 = idx1[0]
    #     iat2 = idx2[0]
    #     if iat1 == iat2:
    #         iat_disp = iat1
    #         disp_iat = D[iat_disp, iat_disp]
    #         if abs(disp_iat - mag*2.) < tolerance:
    #             return True
    # return False

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

