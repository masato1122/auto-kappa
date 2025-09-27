# 
# fcs.py
# 
# This script parses "force constant" section of ALAMODE log file.
# 
# Created on September 16, 2025
# Copyright (c) 2025 Masato Ohnishi
#
# This file is distributed under the terms of the MIT license.
# Please see the file 'LICENCE.txt' in the root directory
# or http://opensource.org/licenses/mit-license.php for information.
# 
from auto_kappa.almlog.utils import (
    replace_symbols_to_blank)

def _read_max_deviation(lines):
    
    title = 'Maximum deviation from the translational invariance'
    
    max_deviations = {}
    for il, line in enumerate(lines):
        if not line:
            continue
        line_low = line.strip().lower()
        if line_low.startswith(title.lower()):
            count = 0
            while True:
                
                if il + 1 + count >= len(lines):
                    break
                
                line_i = lines[il + 1 + count].strip()
                data = line_i.split()
                if len(data) == 0:
                    continue
                
                if data[0].lower() == 'order':
                    order = int(data[1])
                    max_dev = float(data[-1])
                    max_deviations[order] = max_dev
                    count += 1
                else:
                    break
    return max_deviations

def _read_num_nonzero_fcs(lines):
    search_line = "Number of non-zero IFCs for"
    num_nonzero_fcs = {}
    for line in lines:
        if not line:
            continue
        line = replace_symbols_to_blank(line)
        line_low = line.strip().lower()
        if line_low.startswith(search_line.lower()):
            data = line.split()
            if len(data) == 0:
                continue
            order = int(data[-3])
            num_nonzero = int(data[-1])
            num_nonzero_fcs[order] = num_nonzero
    return num_nonzero_fcs

def read_fcs(lines):
    info = {}
    info['max_deviations'] = _read_max_deviation(lines)
    info['num_nonzero_fcs'] = _read_num_nonzero_fcs(lines)        
    return info
