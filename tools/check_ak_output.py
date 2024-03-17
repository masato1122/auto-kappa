#
# check_apdb_log_ver2.py
#
# This script helps to check if the automation calculation was finished or not.
#
# Copyright (c) 2024 Masato Ohnishi
#
# This file is distributed under the terms of the MIT license.
# Please see the file 'LICENCE.txt' in the root directory
# or http://opensource.org/licenses/mit-license.php for information.
#
"""
Hot to use
==========
>>> python check_ak_output.py --dir_apdb ${directory name}

``directory name`` is the directory name of the automation calculation to check.
If the calculation has not yet been finished, the return value is ``NotYet``.

"""
import sys
import os.path
import numpy as np
from optparse import OptionParser
import yaml
import glob
import ase.io

from auto_kappa import output_directories
from auto_kappa.io.vasp import wasfinished as wasfinished_vasp
from auto_kappa.structure.crystal import get_spg_number
from auto_kappa.alamode.io import wasfinished_alamode

POSSIBLE_STATUSES = {
        0: "NotYet",
        1: "Finished_harmonic",
        2: "Finished_cubic",
        3: "Symmetry_error",
        4: "Stop_with_error",
        }

def too_many_symmetry_errors(directory, tol_number=5):
    """ check symmetry error during energy minimization """
    for ext in ["tar", "tar.gz"]:
        dir_error = "%s/relax_error%d.%s" % (directory, tol_number, ext)
        if os.path.exists(dir_error):
            return True
    return False

def change_symmetry(directory):
    """ Check if the log file contains "symmetry change" or not. """
    logfile = directory + "/ak.log"
    if os.path.exists(logfile) == False:
        return False
    lines = open(logfile, 'r').readlines()
    for line in lines:
        if "Error: crystal symmetry was changed" in line:
            return True 
    return False

#def stop_with_error(directory):
#    """ Check if the log file contains "STOP THE CALCULATION" or not. """
#    logfile = directory + "/ak.log"
#    if os.path.exists(logfile) == False:
#        return False
#    
#    lines = open(logfile, 'r').readlines()
#    for line in lines:
#        if "STOP THE CALCULATION" in line:
#            return True
#    return False

def _get_directory_name_for_largesc(dir_orig):
    """ get a directory name for larger SC """
    
    line = dir_orig + "/sc-*"
    dirs = glob.glob(line)
    if len(dirs) == 0:
        return None
    
    ### get the largest SC
    nsc_max = 0
    dir_tmp = None
    for dd in dirs:
        data = dd.split("/")[-1].replace("sc-", "").split("x")
        nsc_i = 0
        for j in range(3):
            nsc_i += int(data[j])
        if nsc_i > nsc_max:
            nsc_max = nsc_i
            dir_tmp = dd
    return dir_tmp

def check_results_with_larger_sc(dir_orig):
    """ check results for larger supercell """
    dir_sc = _get_directory_name_for_largesc(dir_orig)
    if dir_sc is None:
        return None
    else:
        num = check_log_yaml(dir_sc)
        
        statuses = POSSIBLE_STATUSES
        if num in list(statuses.keys()):
            return statuses[num]
        else:
            return "Error"

#def check_log_yaml(dir_name, tol_zero=-1e-3):
#    """
#    Return
#    --------
#    0 : not yet finished
#    1 : finished at harmonic
#    2 : finished at cubic
#    """
#    emin = get_minimum_energy(dir_name)
#    
#    ### Band and DOS were not calculated.
#    if emin is None:
#        if change_symmetry(dir_name):
#            return 3
#        if too_many_symmetry_errors(dir_name):
#            return 3
#        if stop_with_error(dir_name):
#            return 4
#        #if too_many_relaxation_errors(dir_name, max_errors=200):
#        #    return 4
#        return 0
#    
#    ###
#    if emin < tol_zero:
#        """ with negative frequencies """
#        return 1
#    else:
#        """ w/o negative frequencies """
#        figname = dir_name + "/result/fig_kappa.png"
#        if os.path.exists(figname):
#            return 2
#        else:
#            return 0

## >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
## >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

def finish_relaxation(dir_base):
    """ """
    outdir = dir_base + "/" + output_directories["relax"]
    
    #### ver.1
    #xmlfile = outdir + "/vasprun.xml"
    #if os.path.exists(xmlfile):
    #    if wasfinished_vasp(outdir, filename='vasprun.xml'):
    #        return True
    
    ### get initial structure
    filename = outdir + "/full-1/POSCAR.orig"
    struct_orig = ase.io.read(filename, format='vasp')
    num_orig = get_spg_number(struct_orig)

    ### ver.2
    for dd in ["full-1", "full-2", "freeze-1"]:
        dir1 = outdir + "/" + dd
        if os.path.exists(dir1):
            if wasfinished_vasp(dir1) == False:
                print(" %s : %s" % (POSSIBLE_STATUSES[0], dd))
                return False
        else:
            print(" %s : %s" % (POSSIBLE_STATUSES[0], dd))
            return False
        
        ### symmetry check
        filename = dir1 + "/CONTCAR"
        struct = ase.io.read(filename, format='vasp')
        num_i = get_spg_number(struct)
        if num_orig != num_i:
            print(" %s : %s" % (POSSIBLE_STATUSES[3], dir1))
            return False
    ###
    return True

def get_minimum_energy(dir_name):
    """ Check minimum frequency """ 
    from auto_kappa.alamode.log_parser import get_minimum_frequency_from_logfile
    
    fmin = 1000
    for kind in ["band", "dos"]:
        logfile = dir_name + "/harm/bandos/%s.log" % kind
        if os.path.exists(logfile) == False:
            return None
        try:
            fmin_each = get_minimum_frequency_from_logfile(logfile)
            fmin = min(fmin, fmin_each["minimum_frequency"])
        except Exception:
            return None
    return fmin
    
def finish_harmonic(directory, neg_freq=-0.001):
    """ Check calculations for harmonic properties """
    
    ### check log files
    dir_bandos = directory + "/" + output_directories['harm']['bandos']
    for propt in ["band", "dos"]:
        logfile = dir_bandos + "/" + propt + ".log"
        if os.path.exists(logfile) == False:
            print(" %s : %s" % (POSSIBLE_STATUSES[0], dir_bandos))
            return False
        ##
        if wasfinished_alamode(logfile) == False:
            print(" %s : %s" % (POSSIBLE_STATUSES[0], dir_bandos))
            return False
    
    ###
    fmin = get_minimum_energy(directory)
    if fmin < neg_freq:
        print(" %s : %.3f" % (POSSIBLE_STATUSES[1], fmin))
        return False
    
    return True

def finish_cubic(directory):
    """ Check calculations for cubic properties """
    
    dir_cube = directory + "/cube"
    
    ### check
    if os.path.exists(dir_cube) == False:
        print(" %s : %s" % (POSSIBLE_STATUSES[0], dir_cube))
        return False
    
    ### kappa
    line = dir_cube + "/kappa_*"
    dirs = glob.glob(line)
    for dd in dirs:
        
        logfile = dd + "/kappa.log"
        
        ### check presence
        if os.path.exists(logfile) == False:
            print(" %s : %s" % (POSSIBLE_STATUSES[0], dd))
            return False
        
        ### check if finished or not
        if wasfinished_alamode(logfile) == False:
            print(" %s : %s" % (POSSIBLE_STATUSES[0], dd))
            return False
    
    return True

def check_result(directory):
    """ Check output figures"""
    
    dir_result = directory + "/" + output_directories["result"]
    
    names = {
            "harm": ["bandos"],
            "cube": ["cumu_frequency", "cumu_mfp", "kappa", "lifetime", "scat_rates"],
            "time": ["times"]}
    
    for order in ["harm", "cube", "time"]:
        for name in names[order]:
            figname = dir_result + "/fig_%s.png" % name
            if os.path.exists(figname) == False:
                print(" %s : %s" % (POSSIBLE_STATUSES[4], dir_result))
                sys.exit()
    return True

def finish_ak_calculation(directory, neg_freq=-0.001):
    """ Check output directory generated by auto-kappa """
    
    ### log file
    logfile = directory + "/ak.log"
    if os.path.exists(logfile) == False:
        print(" %s : %s" % (POSSIBLE_STATUSES[0], directory))
        sys.exit()
    
    ### check structure optimization
    flag_relax = finish_relaxation(directory)
    if flag_relax == False:
        return False
    
    ### check NAC
    dir_nac = directory + "/" + output_directories["nac"]
    if os.path.exists(dir_nac):
        if wasfinished_vasp(dir_nac) == False:
            print(" %s : %s" % (POSSIBLE_STATUSES[0], dir_nac))
            return False
    
    ### check harmonic
    flag_harm = finish_harmonic(directory, neg_freq=neg_freq)
    if flag_harm == False:
        return False
    
    ### check cubic
    flag_cube = finish_cubic(directory)
    if flag_cube == False:
        return False
    
    ### calculation with supercell
    #finish_sc(directory)

    ### check result
    flag = check_result(directory)
    if flag == False:
        return False

    return True

def main(options):
    
    flag = finish_ak_calculation(
            options.directory, 
            neg_freq=options.neg_freq)

    if flag:
        print(" DONE")

if __name__ == '__main__':
    parser = OptionParser()
    
    parser.add_option("-d", "--directory", dest="directory", type="string",
                help="directory name for APDB")
    
    parser.add_option("--neg_freq", dest="neg_freq", type="float",
                default=-0.001, help="negtive frequency [-0.001]")
    
    (options, args) = parser.parse_args()
    main(options)

