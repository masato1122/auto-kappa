# 
# variables.py
# 
# This script parses variables from ALAMODE log files.
# 
# Created on September 13, 2025
# Copyright (c) 2025 Masato Ohnishi
#
# This file is distributed under the terms of the MIT license.
# Please see the file 'LICENCE.txt' in the root directory
# or http://opensource.org/licenses/mit-license.php for information.
# 
from collections import defaultdict
from auto_kappa.io.alm import proc_val

def read_variables(lines, step=1):
    """ Read variables from log file lines
    """
    result = defaultdict(dict)
    current_section = None
    
    flag_read = False
    count = 0
    for line in lines:
        line = line.strip()
        
        ## 1. Skip empty lines
        if not line: # Skip empty lines
            continue
        
        ## 2. Start reading after "Input variables" line is found
        ##    otherwise, ignore lines before that
        if not flag_read:
            if "Input variables" in line:
                flag_read = True
            continue  # Ignore lines before the variables part
        
        ## 3. Stop reading when the second "------------------" line is found
        if "------------------" in line:
            count += 1
            if count >= 2:
                break
        
        ## Section name ("General", "Kpoint", "Analysis", etc.)
        if line.endswith(':'):
            current_section = line[:-1].strip().split()[0]
            continue
        
        ## 4. Parse key-value pairs
        for entry in line.split(';'):
            if '=' in entry:
                key, value = entry.split('=', 1)
                if value.strip() != "":
                    var_name = key.strip().split()[0].lower()
                    result[current_section.lower()][var_name] = value.strip()
    
    ## Process values (convert to appropriate types)
    params = {}
    for key1 in result.keys():
        params[key1] = {}
        for key2 in result[key1].keys():
            val = result[key1][key2]
            try:
                params[key1][key2] = proc_val(key2, val)
            except:
                params[key1][key2] = val
    
    if step == 2:
        return params
    elif step == 1:
        # Flatten the dictionary
        flat_params = {}
        for key1 in params.keys():
            for key2 in params[key1].keys():
                flat_params[key2] = params[key1][key2]
        return flat_params
    
