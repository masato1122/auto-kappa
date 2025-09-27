# 
# utils.py
# 
# Created on September 13, 2025
# Copyright (c) 2025 Masato Ohnishi
#
# This file is distributed under the terms of the MIT license.
# Please see the file 'LICENCE.txt' in the root directory
# or http://opensource.org/licenses/mit-license.php for information.
# 
import numpy as np
import re

def parse_data_line(line, key=None, idx_dist=None, index=None, dtype=float):
    """ Extract data from a line of log file 
    
    Algorithms:
        1. Split the line into words
        2. Search for the keyword in the words
        3. If found, get the value at the specified index distance from the keyword
        4. If last is True, get the last value in the line
    
    Args:
        line (str): a line of log file
        key (str): keyword to search
        idx_dist (int): index distance from the keyword to the target value
        last (bool): if True, get the last value in the line
        dtype (type): data type of the target value (default: float)
    
    """
    if idx_dist is None and index is None:
        raise ValueError("Either idx_dist or last must be specified.")
    
    data = line.strip().split()
    
    if index is not None:
        try:
            value = dtype(data[index])
        except:
            value = None
    else:
        if key in data:
            idx_key = data.index(key)
            try:
                value = dtype(data[idx_key+idx_dist])
            except:
                value = None
        else:
            value = None
    
    return value

def replace_symbols_to_blank(line, symbols=["=", ",", ":", ";", "(", ")", "[", "]"]):
    for sym in symbols:
        line = line.replace(sym, " ")
    return line

def extract_afterward_lines(lines, title, num_lines=None):
    """ Extract lines after a specific title line for a given number of lines.
    
    Example:
        Mass of atomic species (u):
        S:   32.065000
        Ba:  137.327000
    
    Args:
        lines (list): list of lines in the log file
        title (str): title line to search
        num_lines (int): number of lines to read after the title line
    
    Returns:
        list: list of extracted lines
    """
    extracted_lines = []
    for il, line in enumerate(lines):
        if line.strip().lower().startswith(title.lower()):
            for i in range(num_lines):
                try:
                    extracted_lines.append(lines[il + 1 + i])
                except IndexError:
                    print(f" Error reading lines after '{title}'.")
                    break
            break
    return extracted_lines

def extract_lines_with_index(lines, title, initial_index=1, replace_symbols=True):
    """ Extract lines after a specific title line that start with an integer index (1, 2, 3, ...)
    The extraction stops when a line does not start with a sequential integer index.
    Symbols such as '=', ':', etc. are replaced with blanks if replace_symbols is True.
    
    Example:
    
    
    Args:
        lines (list): list of lines in the log file
        title (str): title line to search
        initial_index (int): initial index for the ordered list (default: 1)
        replace_symbols (bool): if True, replace symbols with blanks
    """
    extracted_lines = None
    for il, line in enumerate(lines):
        if line.strip().lower().startswith(title.lower()):
            extracted_lines = []
            count = 0
            
            try:
                while True:
                    
                    ## Check index range
                    if il + count + 1 >= len(lines):
                        break
                    
                    if replace_symbols:
                        l = replace_symbols_to_blank(lines[il + count + 1])
                    else:
                        l = lines[il + count + 1]
                    parts = l.strip().split()
                    
                    ## Check if the first part is an integer
                    if not parts[0].isdigit():
                        break
                    
                    ## Check if the index is sequential
                    idx_now = int(parts[0])
                    if idx_now != count + initial_index:
                        break
                    
                    extracted_lines.append(l)
                    count += 1
            
            except IndexError:
                print(f" Error reading lines after '{title}'.")
            
            break
    return extracted_lines

def get_position_list(lines, title, initial_index=1, replace_symbols=True):
    """ Parse ordered list (e.g. atomic positions) after a specific title line
    Each line must start with an integer index (1, 2, 3, ...). 
    If the index is not sequential, the parsing stops.
    Currently, each line must contain five or more columns: index, x, y, z, (...,) tag
    Symbols such as '=', ':', etc. are replaced with blanks if replace_symbols is True.
    
    Example:
        Atomic positions in fractional basis and atomic species
        1   3.403998e-01   4.096002e-01   1.208532e-01    1
        2   8.403998e-01   4.096002e-01   1.208532e-01    1
        3   3.403998e-01   9.096002e-01   1.208532e-01    1
        4   8.403998e-01   9.096002e-01   1.208532e-01    1
        5   3.403998e-01   4.096002e-01   3.708532e-01    1
        ...
        
    Args:
        lines (list): list of lines in the log file
        title (str): title line to search
        initial_index (int): initial index for the ordered list (default: 1)
        replace_symbols (bool): if True, replace symbols with blanks
    
    Returns:
        list: list of parsed values
    """
    extracted_lines = extract_lines_with_index(lines, title, initial_index=initial_index, 
                                               replace_symbols=replace_symbols)
    
    if extracted_lines is None or len(extracted_lines) == 0:
        return None, None, None
    
    indices = []
    coords = []
    tag_list = []
    
    for line in extracted_lines:
        parts = line.strip().split()
        if len(parts) < 5:
            continue
        try:
            indices.append(int(parts[0]))
            coords.append([float(v) for v in parts[1:4]])
            tag_list.append(parts[4])
        except ValueError:
            print(f" Error parsing line: {line}")
            continue
    
    return (np.asarray(indices), np.asarray(coords), tag_list)

def read_section(lines, search_lines1={}, search_lines2={}):
    """ Read section information from log file lines (e.g. Symmetry section)
    
    Args:
        lines (list): list of lines in the log file
        search_lines1 (dict): additional search patterns for type 1 (direct string search)
        search_lines2 (dict): additional search patterns for type 2 (regex)
            The value must be a tuple of (pattern, dtype) where pattern is a regex pattern
            and dtype is a type or a tuple of types for multiple values.
    
    Returns:
        info (dict): dictionary containing section information
        
    """
    info = {}
    for il, line in enumerate(lines):
        
        line = line.strip()
        if not line:
            continue
        
        ### Type 1: direct string search
        for key in search_lines1:
            if line.lower().startswith(search_lines1[key][0].lower()):
                val = parse_data_line(line, index=-1, dtype=search_lines1[key][1])
                if val is not None:
                    info[key] = val
                continue
        
        ### Type 2: regex search
        for key, (pattern, dtype) in search_lines2.items():
            match = re.search(pattern, line)
            if match:
                if isinstance(dtype, tuple):
                    val = tuple(d(match.group(i+1)) for i, d in enumerate(dtype))
                else:
                    val = dtype(match.group(1))
                if val is not None:
                    info[key] = val
                continue
    
    return info
