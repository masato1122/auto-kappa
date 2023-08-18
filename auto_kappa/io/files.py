# -*- coding: utf-8 -*-
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
            data[name] = info
            f.seek(0)
            yaml.dump(data, f, sort_keys=False)

