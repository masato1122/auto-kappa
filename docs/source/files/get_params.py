# -*- coding: utf-8 -*-
import numpy as np
from optparse import OptionParser
import pandas as pd

def main(options):
    
    from auto_kappa import default_vasp_parameters as params
    
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
        if key in ['shared', 'relax']:
            continue
        each = []
        for name in names:
            
            if name in params['shared']:
                each.append(params['shared'][name])
            else:
                each.append('_')
            
            if 'relax' in name:
                if name in params['relax']:
                    each.append(params['relax'][name])
                else:
                    each.append('_')

            if name in params[key]:
                each[-1] = params[key][name]
        
        df[key] = each
    
    #print(df)
    outfile = "default_params.csv"
    df.to_csv(outfile)
    print(" Output", outfile)
    
if __name__ == '__main__':
    parser = OptionParser()
    parser.add_option("-f", "--filename", dest="filename", type="string",
        help="input file name")
    (options, args) = parser.parse_args()
    #file_check(options.filename)
    main(options)

