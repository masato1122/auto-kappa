# -*- coding: utf-8 -*-
import sys
import os.path
import numpy as np
import pandas as pd
import yaml

def write_output_yaml(filename, name, info, overwrite=True):
    """ Write an output directory name in a yaml file.

    Args
    =====
    filename : string
        yaml file name
    name : string
        name of a dictionary
    info : dict
        dictionary
    """
    if os.path.exists(filename) == False or overwrite == False:
        with open(filename, 'w') as f:
            data = {}
            data[name] = info
            yaml.dump(data, f, sort_keys=False)
    else:
        with open(filename, 'r+') as f:
            data = yaml.safe_load(f)
            name_add = name

            ### overwritten
            note_orig = info["note"]
            count = 2
            while name_add in data.keys():
                name_add = name + "(%d)" % count
                info["note"] = note_orig + ", overwritten"
                count += 1
            
            data[name_add] = info
            f.seek(0)
            yaml.dump(data, f, sort_keys=False)
    f.close()

def extract_data_from_file(filename, word, index=-1):
    """ Return data for ``word`` in ``filename``
    """
    lines = open(filename, 'r').readlines()
    values_get = []
    for line in lines:
        if word.lower() in line.lower():
            data = line.strip().split()[index]
            values_get.append(float(data))
    if len(values_get) == 0:
        return None
    else:
        return values_get

def convert_numpy(obj):
    """ Convert numpy data types to native Python types.
    
    How to use
    ----------
    >>> with open("output.json", "w", encoding="utf-8") as f:
    >>>     json.dump(data, f, indent=4, ensure_ascii=False, default=convert_numpy)
    """
    if isinstance(obj, (np.integer,)):
        return int(obj)
    elif isinstance(obj, (np.floating,)):
        return float(obj)
    elif isinstance(obj, (np.ndarray,)):
        return obj.tolist()
    return obj
