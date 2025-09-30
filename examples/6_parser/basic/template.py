#!/usr/bin/env python3
import os
import argparse
import glob

def _get_formula(base_dir):
    line = f"{base_dir}/harm/bandos/*.bands"
    fns = glob.glob(line)
    if len(fns) == 0:
        raise ValueError(" No band file found")
    fn = fns[0]
    return os.path.basename(fn).split(".")[0]

def main(options):
    
    formula = _get_formula(options.base_dir)
    

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Input parameters')
    parser.add_argument('--base_dir', dest='base_dir', type=str,
                        default="./mp-22862", help="directory of auto-kappa output")
    args = parser.parse_args()    
    main(args)
    