# 
# dmatrix.py
# 
# This script parses "dynamical matrix" section of ALAMODE log file.
# 
# Created on September 16, 2025
# Copyright (c) 2025 Masato Ohnishi
#
# This file is distributed under the terms of the MIT license.
# Please see the file 'LICENCE.txt' in the root directory
# or http://opensource.org/licenses/mit-license.php for information.
# 
import numpy as np
import re
from auto_kappa.almlog.utils import (
    extract_afterward_lines
    )

def _get_dielectric(lines, title):
    extract_lines = extract_afterward_lines(lines, title=title, num_lines=3)
    dielectric = []
    for line in extract_lines:
        data = line.strip().split()
        if len(data) == 0:
            continue
        row = [float(v) for v in data]
        dielectric.append(row)
    return np.array(dielectric)

def _get_born_charges(lines, title):
    
    v_int = r"\s*(\d+)\s*"
    v_str = r"\s*(\S+)\s*"
    line0_style = fr"Atom{v_int}\({v_str}\)\s*:"

    elements = []
    born_charges = []
    for il, line in enumerate(lines):
        if not line:
            continue
        line_low = line.strip().lower()
        if line_low.startswith(title.lower()):
            count = 0
            while True:
                
                ## Break if the next atom line is not found
                if il + 1 + 4*count >= len(lines):
                    break
                
                ## Check the line style for atom index
                il0 = il + 1 + 4*count
                line0 = lines[il0].strip()
                match = re.match(line0_style, line0)
                if not re.match(line0_style, line0):
                    break
                
                index = int(match.group(1))
                element = match.group(2)
                elements.append(element)
                
                ## Read born charge for each atom
                for i in range(3):
                    line_i = lines[il0 + 1 + i].strip()
                    data = line_i.split()
                    if len(data) == 0:
                        continue
                    row = [float(v) for v in data]
                    born_charges.append(row)
                
                count += 1
    
    return (elements, np.array(born_charges))

def read_dmatrix(lines):
    """ Read dynamical matrix information from log file lines (Dynamical matrix section)
    
    Args:
        lines (list): list of lines in the log file
    Returns:
        info (dict): dictionary containing dynamical matrix information 
    """
    info = {}
    
    for line in lines:
        line = line.strip()
        data = line.split()
        if len(data) == 0:
            continue
        
        ## Non-analytic correction
        if line.lower().startswith("nonanalytic ="):
            info['nonanalytic'] = int(data[2])
        
        ## Dielectric constant tensor
        title = "Dielectric constant tensor in Cartesian coordinate"
        if line.lower().startswith(title.lower()):
            info['dielectric_tensor'] = _get_dielectric(lines, title)
        
        ## Born effective charge tensor
        title = "Born effective charge tensor in Cartesian coordinate"
        if line.lower().startswith(title.lower()):
            info['born_charges'] = _get_born_charges(lines, title)
    
    return info
