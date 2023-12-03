#
# check_apdb_log.py
#
# This script helps to check if the automation calculation was finished or not.
#
# Copyright (c) 2023 Masato Ohnishi
#
# This file is distributed under the terms of the MIT license.
# Please see the file 'LICENCE.txt' in the root directory
# or http://opensource.org/licenses/mit-license.php for information.
#
"""
Hot to use
==========
>>> python check_apdb_log.py --dir_apdb ${directory name}

``directory name`` is the directory name of the automation calculation to check.
If the calculation has not yet been finished, the return value is ``NotYet``.

"""
import os.path
import numpy as np
from optparse import OptionParser
import yaml
import glob

POSSIBLE_STATUSES = {
        0: "NotYet",
        1: "Finished_harmonic",
        2: "Finished_cubic",
        3: "Symmetry_error",
        4: "Stop_with_error",
        }

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
    
def too_many_symmetry_errors(directory, tol_number=5):
    """ check symmetry error during energy minimization """
    dir_error = "%s/relax_error%d.tar.gz" % (directory, tol_number)
    if os.path.exists(dir_error):
        return True
    else:
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

def symmetry_error(directory):
    """ Check if the log file contains "STOP THE CALCULATION" or not. """
    logfile = directory + "/ak.log"
    if os.path.exists(logfile) == False:
        return False
    
    lines = open(logfile, 'r').readlines()
    for line in lines:
        if "STOP THE CALCULATION" in line:
            return True
    return False

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

#def too_many_relaxation_errors(dir_name, max_errors=None):
#    """ Check number of errors for VASP calculation """
#    for dd in ["freeze-1", "full-2", "full-1"]:
#        line = dir_name + "/relax/" + dd + "/error.*.tar*"
#        dirs = glob.glob(line)
#        for deach in dirs:
#            num = int(deach.replace(dir_name, "").split("error.")[-1].split(".")[0])
#            if num >= max_errors:
#                return True
#        if len(dirs) > 0:
#            return False
#    return False
    
def check_log_yaml(dir_name, tol_zero=-1e-3):
    """
    Return
    --------
    0 : not yet finished
    1 : finished at harmonic
    2 : finished at cubic
    """
    emin = get_minimum_energy(dir_name)
    
    ### Band and DOS were not calculated.
    if emin is None:
        if change_symmetry(dir_name):
            return 3
        if too_many_symmetry_errors(dir_name):
            return 3
        if stop_with_error(dir_name):
            return 4
        #if too_many_relaxation_errors(dir_name, max_errors=200):
        #    return 4
        return 0
    
    ###
    if emin < tol_zero:
        """ with negative frequencies """
        return 1
    else:
        """ w/o negative frequencies """
        figname = dir_name + "/result/fig_kappa.png"
        if os.path.exists(figname):
            return 2
        else:
            return 0

def main(options):

    #flag_sc = check_results_with_larger_sc(options.dir_apdb)
    
    flag = check_log_yaml(options.dir_apdb)
    
    if flag in list(POSSIBLE_STATUSES.keys()):
        print(POSSIBLE_STATUSES[flag])
    else:
        print("Error")

if __name__ == '__main__':
    parser = OptionParser()
    
    parser.add_option("-d", "--dir_apdb", dest="dir_apdb", type="string",
            help="directory name for APDB")
    
    (options, args) = parser.parse_args()
    main(options)

