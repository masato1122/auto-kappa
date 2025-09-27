# 
# kpoint.py
#
# This script parses "kpoint" section of ALAMODE log file.
#
# Created on September 16, 2025
# Copyright (c) 2025 Masato Ohnishi
#
# This file is distributed under the terms of the MIT license.
# Please see the file 'LICENCE.txt' in the root directory
# or http://opensource.org/licenses/mit-license.php for information.
# 
from auto_kappa.almlog.utils import (
    replace_symbols_to_blank, extract_afterward_lines,
    extract_lines_with_index, 
    read_section)

def _read_kpaths(lines):
    kpath_lines = extract_lines_with_index(lines, title="List of k paths")
    if kpath_lines is None:
        return None

    kpaths = []
    for line in kpath_lines:
        data = line.strip().split()
        
        ## index, name0, k0, name1, k1, nk
        index = int(data[0])
        name0 = data[1]
        k0 = [float(v) for v in data[2:5]]
        name1 = data[5]
        k1 = [float(v) for v in data[6:9]]
        nk = int(data[9])
        kpaths.append([{name0: k0}, {name1: k1}, nk])
    
    return kpaths

def _get_nkpts(lines, title):
    extracted_lines = extract_afterward_lines(lines, title=title, num_lines=3)
    nkpts = []
    for line in extracted_lines:
        nk = int(line.strip().split()[-1])
        nkpts.append(nk)
    return nkpts  

def _get_irred_kpts(lines, title):
    extract_lines = extract_lines_with_index(lines, title=title)
    reciprocal_kpoints = []
    weights = []
    for i, line in enumerate(extract_lines):
        data = line.split()
        index = int(data[0])
        if index != i+1:
            print(f" Warning: Irreducible k-point index mismatch: {index} != {i+1}")
        kpt = [float(v) for v in data[1:4]]
        weight = float(data[4])
        reciprocal_kpoints.append(kpt)
        weights.append(weight)
    return {
        'reciprocal_kpoints': reciprocal_kpoints,
        'weights': weights
    }

def read_kpoints(lines, search_lines_direct={}, search_lines_regex={}):
    """ Read k-point information from log file lines (kpoint section)

    Args:
        lines (list): list of lines in the log file
        search_lines_direct (dict): additional search patterns (direct string search)
        search_lines_regex (dict): additional search patterns (regex)
            The value must be a tuple of (pattern, dtype) where pattern is a regex pattern
            and dtype is a type or a tuple of types for multiple values.
    
    Returns:
        info (dict): dictionary containing symmetry information
    
    """
    info = {}
    
    for line in lines:
        
        line = line.strip().lower()
        data = line.split()
        if len(data) == 0:
            continue
        
        ## kpmode
        if line.lower().startswith("kpmode"):
            try:
                data = replace_symbols_to_blank(line).strip().split()
                info['kpmode'] = int(data[1])
            except (IndexError, ValueError):
                print(f" Error reading kpmode: {line}")
                pass
        
        ## nkpts
        title = "Gamma-centered uniform grid"
        if line.lower().startswith(title.lower()):
            info['nkpts'] = _get_nkpts(lines, title)

        ## List of irreducible k points
        title = "list of irreducible k points"
        if line.lower().startswith(title.lower()):
            info['irreducible_kpoints'] = _get_irred_kpts(lines, title)

    ## k-path
    info['kpath'] = _read_kpaths(lines)
    
    ## Search lines for type 1 (direct string search)
    search_lines1 = {
        "num_kpoints": ("Number of k points", int),
        "num_kpoints_irred": ("Number of irreducible k points", int),
    }
    search_lines1.update(search_lines_direct)

    ## Search lines for type 2 (regex)
    search_lines2 = {}
    search_lines2.update(search_lines_regex)
    
    info.update(read_section(lines, search_lines1, search_lines2))
    
    return info
