## 
## check_old_nac.py
##
## This script helps to compare the k-mesh for NAC calculation in old
## calculations and the suggested value by the new version of auto-kappa.
##
## If you switch on the option "print_only" (please see below), directories in
## the old calculation, "nac", "harm/bandos", "harm/evec", "cube/kappa_*", will
## be moved to "{output directory name such as mp-***}/small_kpts4nac".
## 
"""
How to Use
-----------
>>> cd (directory in which many calculations have done with auto-kappa)
>>> ls 
... mp-*** mp-*** mp-*** ...
>>> python check_old_nac.py
(This will only print Material IDs whose NAC caluclation was conducted with an inproper k-mesh.)
mp-***
mp-***
mp-***
...

To show the result only for a certain number of materials:

>>> python check_old_nac.py --nmater 5 
mp-***
mp-***
mp-***
...

To move directories related to NAC calculation:

>>> python check_old_nac.py --nmater 5 --print_only 0
mp-***
mp-***
mp-***
...

"""
# -*- coding: utf-8 -*-
import os.path
import numpy as np
from optparse import OptionParser

from pymatgen.core.structure import Structure
from pymatgen.io.vasp.inputs import Kpoints
from auto_kappa.cui.suggest import klength2mesh

def get_oldnac_info(dir_apdb):
    """ Read info of the old calculation from the following files
    - ./mp-***/nac/POSCAR
    - ./mp-***/nac/KPOINTS
    
    Args
    -----
    dir_apdb : string
        directory name for an old calculation such as ./mp-***
    
    Return
    -------
    structure :
        structure with Pymatgen format

    kpoints : list, shape=(3)
        k-mesh

    None : if information related to NAC cannot be extracted.
    
    """
    dir_nac = dir_apdb + "/nac"
    fn_pos = dir_nac + "/POSCAR"
    fn_kpoints = dir_nac + "/KPOINTS"
    flag_dir = os.path.exists(dir_nac)
    flag_pos = os.path.exists(fn_pos)
    flag_kpts = os.path.exists(fn_kpoints)

    if flag_dir:
        if flag_pos:
            structure = Structure.from_file(fn_pos)
            if flag_kpts:
                kpoints = Kpoints.from_file(fn_kpoints).kpts[0]
                return structure, kpoints
            else:
                print(" Error: KPOINTS file cannot be found.")
                return None
        else:
            print(" Error: POSCAR file cannot be found.")
            return None
    else:
        #print(" No NAC calculation")
        return None

def check_each_nac(dir_apdb):
    """ Compare the previously-used k-mesh (kpts_pre) and the suggested k-mesh 
    by the new version of auto-kappa (kpts_sug). If kpts_pre < kpts_sug, return 
    1.
    
    Args
    -----
    dir_apdb : string
        Output directory name of auto-kappa, which is given by the option,
        "--material_name"
    
    Return
    -------
    0 : previously used k-mesh is not too small.
    1 : previously used k-mesh is too small compared with the suggested value.

    """
    ### Old calculation
    out = get_oldnac_info(dir_apdb)
    
    if out is None:
        return -1
    
    structure = out[0]
    kpts_pre = out[1]
    
    ### Suggested value by auto-kappa
    k_length = 20.     ## default value in auto-kappa
    kpts_sug = klength2mesh(k_length, structure.lattice.matrix)
    
    ### Compare k-meshes
    flag_smaller = False
    for j in range(3):
        if kpts_pre[j] < kpts_sug[j]:
            flag_smaller = True
    
    if flag_smaller == False:
        #print(" No problem.")
        return 0
    else:
        #print(" suggested: ", kpts_sug, ", previous: ", kpts_pre)
        return 1

def move_directories(dir_apdb, prefix_new="small_kpts4nac"):
    """
    Under ``dir_apdb``, the following directories will be moved to
    ``small_kpts4nac``.
    ./nac
    ./harm/bandos
    ./harm/evec
    ./cube/kappa_*
    """
    import subprocess
    import glob

    dir_new = dir_apdb + "/" + prefix_new
    
    ##
    dir_moved = [
            ["nac"],
            ["harm", "bandos"],
            ["harm", "evec"],
            ["cube", "kappa_*"]
            ]
    
    for dd in dir_moved:
        
        ###
        if dd[0] == "nac":
            dir_pre = dir_apdb + "/" + dd[0]
            dir_move = dir_new
        else:
            dir_pre = dir_apdb + "/" + dd[0] + "/" + dd[1]
            dir_move = dir_new + "/" + dd[0]
        
        dirs = glob.glob(dir_pre)
        for dir_each in dirs:
            cmd = "mv %s %s/" % (dir_each, dir_move)
            os.makedirs(dir_move, exist_ok=True)
            try:
                #print(cmd)
                subprocess.run(cmd.split())
            except Exception:
                pass

def main(options):
    
    import glob
    line = options.prefix + "*"
    dirs = glob.glob(line)
    count = 0
    for dir_each in dirs:
        
        flag = check_each_nac(dir_each)
        
        if flag == 1:
            print(dir_each)
            if options.print_only == 0:
                move_directories(dir_each)
            
            count += 1
        
        ## stop at the middle of the process
        if options.nmater is not None:
            if count == options.nmater:
                exit()

if __name__ == '__main__':
    
    parser = OptionParser()
    
    parser.add_option("--prefix", dest="prefix", type="string",
            default="./mp-", help="prefix for old directories [./mp-]")
    
    parser.add_option("--nmater", dest="nmater", type="int",
            default=None, help="Number of materials to be treated [None]")
    
    #parser.add_option("--k_length", dest="k_length", type="float",
    #        default=20, help="k length [20]")
    
    parser.add_option("--print_only", dest="print_only", type="int",
            default=1, 
            help="only print message if old k-mesh is too small (0 or 1) [1]")
    
    (options, args) = parser.parse_args()
    main(options)

