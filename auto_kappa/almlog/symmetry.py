# 
# symmetry.py
#
# This script parses "symmetry" section of ALAMODE log file.
#
# Created on September 14, 2025
# Copyright (c) 2025 Masato Ohnishi
#
# This file is distributed under the terms of the MIT license.
# Please see the file 'LICENCE.txt' in the root directory
# or http://opensource.org/licenses/mit-license.php for information.
# 
import re
from collections import defaultdict
from auto_kappa.almlog.utils import read_section

def _read_cell_atom_mapping(lines):
    """ Read Cell-Atom mapping from lines in the log file
    
    Args:
        lines (list): list of lines in the log file
    """
    mapping = defaultdict(list)
    current_cell = None
    
    title_format = r"CELL\s*|\s*ATOM"
    
    read_data = False
    for line in lines:
        line = line.strip()
        
        if not line:
            continue
        if re.match(title_format, line):
            read_data = True
            continue
        if not read_data:
            continue

        ## line: "cell index | integer ..."
        m = re.match(r"(\d+)\s*\|\s*(.*)", line)
        if m:
            current_cell = int(m.group(1))
            atoms = [int(x) for x in m.group(2).split()]
            mapping[current_cell].extend(atoms)
        # continue to the next line if the line starts with '|'
        else:
            m2 = re.match(r"\|\s*(.*)", line)
            if m2 and current_cell is not None:
                atoms = [int(x) for x in m2.group(1).split()]
                mapping[current_cell].extend(atoms)
            else:
                break
    
    return dict(mapping)

def read_symmetry(lines, search_lines_direct={}, search_lines_regex={}):
    """ Read symmetry information from log file lines (Symmetry section)
    
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
    
    ## Search lines for type 1 (direct string search)
    search_lines1 = {
        "num_sym_ops": ("Number of symmetry operations", int),
    }
    search_lines1.update(search_lines_direct)
    
    ## Search lines for type 2 (regex)
    v_int = r"\s*(\d+)\s*"
    v_str = r"\s*(\S+)\s*"
    search_lines2 = {
        "space_group"     : (fr"Space group:{v_str}\({v_int}\)", (str, int)),
        "num_trans_ops": (fr"There are {v_int} translation operations", int),
        "num_atoms_prim"     : (fr"Primitive cell contains {v_int} atoms", int),
    }
    search_lines2.update(search_lines_regex)
    
    info.update(read_section(lines, search_lines1, search_lines2))
    
    for line in lines:
        title = "**Cell-Atom Correspondens Below**"
        if line.lower().strip().startswith(title.lower()):
            info['cell_atom_mapping'] = _read_cell_atom_mapping(lines)
    
    return info
