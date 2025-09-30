#!/usr/bin/env python3
import os
import argparse
import glob

from auto_kappa.io import BORNINFO

def main(options):
    
    from plot_bandos import get_formula
    formula = get_formula(options.base_dir)
    file_vasprun = f"{options.base_dir}/nac/vasprun.xml"
    file_fc2 = f"{options.base_dir}/harm/force/{formula}.xml"
    
    # BORNINFO object creats BORNINFO file from vasprun.xml 
    # adjusting the order of atoms for ALAMODE force constants file.
    born = BORNINFO(file_vasprun, file_fcs=file_fc2)
    born.write(outfile="BORNINFO")
    
    print("Primitive cell from vasprun.xml")
    print("===============================")
    print(born.prim_vasp)
    print()
    print("Primitive cell from ALAMODE force constants file")
    print("================================================")
    print(born.prim_fcs)
    
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Input parameters')
    parser.add_argument('--base_dir', dest='base_dir', type=str,
                        default="./mp-22862", help="directory of auto-kappa output")
    args = parser.parse_args()    
    main(args)
    