# -*- coding: utf-8 -*-
import numpy as np
from optparse import OptionParser
import pandas as pd

def main(options):
    
    from auto_alamode2 import default_vasp_parameters as params
    
    #######################
    keys_skip = ['md']
    #######################
    
    keys = []
    for key in params:
        if key in keys_skip:
            continue
        keys.append(key)
    
    names = []
    for key in keys:
        if key in keys_skip:
            continue
        for name in params[key]:
            if name not in names:
                names.append(name)
    ##
    df = pd.DataFrame(index=[name.upper() for name in names])
    for key in keys:
        if key == 'shared':
            continue
        each = []
        for name in names:
            if name in params['shared']:
                each.append(params['shared'][name])
            else:
                each.append('_')

            if name in params[key]:
                each[-1] = params[key][name]
        
        df[key] = each
    
    #print(df)
    df.to_csv('default_params.csv')
    
if __name__ == '__main__':
    parser = OptionParser()
    parser.add_option("-f", "--filename", dest="filename", type="string",
        help="input file name")
    (options, args) = parser.parse_args()
    #file_check(options.filename)
    main(options)

