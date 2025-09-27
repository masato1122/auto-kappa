# 
# frequencies.py
# 
# This script obtains phonon frequencies from ALAMODE log files.
# 
# Created on September 16, 2025
# Copyright (c) 2025 Masato Ohnishi
#
# This file is distributed under the terms of the MIT license.
# Please see the file 'LICENCE.txt' in the root directory
# or http://opensource.org/licenses/mit-license.php for information.
# 
import re
import numpy as np
from auto_kappa.almlog.utils import extract_afterward_lines

def read_frequencies(lines):
    """ Read phonon frequencies from ALAMODE log file lines.
    """
    info = {}
    
    in_section = False
    
    blank = r"\s*"
    v_int = fr"{blank}(\d+){blank}"
    # v_str = fr"{blank}(\S+){blank}"
    v_number = fr"{blank}([-+]?(?:\d*\.\d+|\d+\.?)(?:[eE][-+]?\d+)?){blank}"
    
    # k point   422 : (  0.0000,  0.5000,  0.2959)
    title_patterns = [
        fr"# k point {v_int}:{blank}\({v_number},{v_number},{v_number}\)",
        fr"# Irred. k point {v_int}:{blank}\({v_number},{v_number},{v_number}\)",
    ]
    # Example: 1      0.0000 cm^-1  (      0.0000 THz )
    data_pattern = fr"{v_int}{v_number}cm\^-1{blank}\({blank}{v_number}THz{blank}\)"
    
    search_line = "Phonon frequencies below"
    kpoints = []
    frequencies = []
    for il, line in enumerate(lines):
        
        line = line.strip()
        
        if not line:
            continue
        
        if line.lower().startswith(search_line.lower()):
            in_section = True
            continue
        
        if in_section:
            if line.startswith('---') or line.startswith('==='):
                in_section = False
        
        if in_section:
            
            ## Get kpoints or irreducible kpoints
            for i in range(2):
                match_title = re.match(title_patterns[i], line)
                if match_title:
                    kpoints.append([])
                    frequencies.append([])
                    idx_k = match_title.group(1)
                    for j in range(3):
                        kpoints[-1].append(float(match_title.group(2+j)))
                    
            ## Get frequency data lines
            match_data = re.match(data_pattern, line)
            if match_data:
                try:
                    frequencies[-1].append(float(match_data.group(2)))
                except ValueError:
                    continue        
        
    info['kpoints'] = np.array(kpoints)
    info['frequencies'] = np.array(frequencies)
    return info
