#!/usr/bin/env python3
# Filename: make_dfset.py
import argparse

from auto_kappa.io.vasp import get_dfset

def main(options):
    
    from plot_bandos import get_formula
    formula = get_formula(options.base_dir)    
    dir_force = f"{options.base_dir}/cube/force_fd"
    offset_xml = dir_force + '/prist/vasprun.xml'
    
    outfile = 'DFSET.cube'
    disps, force = get_dfset(dir_force, offset_xml=offset_xml, outfile=outfile)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Input parameters')
    parser.add_argument('--base_dir', dest='base_dir', type=str,
                        default="./mp-22862", help="directory of auto-kappa output")
    args = parser.parse_args()    
    main(args)
    