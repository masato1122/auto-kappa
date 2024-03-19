#
# check_ak_output.py
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
>>> dir_apdb = "mp-149"
>>> dir_apdb = "mp-149.tar.gz"
>>> python check_ak_output.py \
        --dir_apdb dir_apdb \
        --vol_relax 0 \
        --larger_sc 0

Note
=====

``dir_apdb`` is the directory name of the automation calculation to check.
If the calculation has not yet been finished, the return value is ``NotYet``.

"""
import sys
import os.path
import numpy as np
from optparse import OptionParser
import yaml
import glob
import ase.io
import tarfile

from auto_kappa import output_directories
from auto_kappa.io.vasp import wasfinished as wasfinished_vasp
from auto_kappa.structure.crystal import get_spg_number
from auto_kappa.alamode.io import wasfinished_alamode

### log definition
import logging
logger = logging.getLogger(__name__)
sh = logging.StreamHandler()
sh.setFormatter(logging.Formatter("%(message)s"))
logging.basicConfig(level=logging.DEBUG, handlers=[sh])

###
POSSIBLE_STATUSES = {
        0: "NotYet",
        1: "Finished_harmonic",
        2: "Finished_cubic",
        3: "Symmetry_error",
        4: "Stop_with_error",
        "sc": "Finished_sc",
        }

#def too_many_symmetry_errors(directory, tol_number=5):
#    """ check symmetry error during energy minimization """
#    for ext in ["tar", "tar.gz"]:
#        dir_error = "%s/relax_error%d.%s" % (directory, tol_number, ext)
#        if os.path.exists(dir_error):
#            return True
#    return False
#
#def change_symmetry(directory):
#    """ Check if the log file contains "symmetry change" or not. """
#    logfile = directory + "/ak.log"
#    if os.path.exists(logfile) == False:
#        return False
#    lines = open(logfile, 'r').readlines()
#    for line in lines:
#        if "Error: crystal symmetry was changed" in line:
#            return True 
#    return False
#
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
#
#def _get_directory_name_for_largesc(dir_orig):
#    """ get a directory name for larger SC """
#    
#    line = dir_orig + "/sc-*"
#    dirs = glob.glob(line)
#    if len(dirs) == 0:
#        return None
#    
#    ### get the largest SC
#    nsc_max = 0
#    dir_tmp = None
#    for dd in dirs:
#        data = dd.split("/")[-1].replace("sc-", "").split("x")
#        nsc_i = 0
#        for j in range(3):
#            nsc_i += int(data[j])
#        if nsc_i > nsc_max:
#            nsc_max = nsc_i
#            dir_tmp = dd
#    return dir_tmp
#
#def check_results_with_larger_sc(dir_orig):
#    """ check results for larger supercell """
#    dir_sc = _get_directory_name_for_largesc(dir_orig)
#    if dir_sc is None:
#        return None
#    else:
#        num = check_log_yaml(dir_sc)
#        
#        statuses = POSSIBLE_STATUSES
#        if num in list(statuses.keys()):
#            return statuses[num]
#        else:
#            return "Error"
#
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
def _read_poscar(filename, tar=None, file_tmp="_tmp.txt"):
    """ Read POSCAR file to get a structure object
    
    Args
    -----
    filename : string
    
    tar : tarfile.TarFile
    
    """
    if tar is None:
        file_pos = filename
    else:
        try:
            content = tar.extractfile(filename).read().decode('utf-8')
            with open(file_tmp, "w") as f:
                f.write(content)
            file_pos = file_tmp
        except Exception:
            print(" Error: %s" % filename)
            return None
    try:
        return ase.io.read(file_pos, format='vasp')
    except Exception:
        return None

def finish_relaxation(dir_base, tar=None):
    """ 
    Args
    -----
    dir_base : string
        Name of a directory or .tar.gz file
    
    tar : tarfile.TarFile
    """
    dir_harm = dir_base + "/harm"
    if _exists(dir_harm, tar=tar):
        return True
    else:
        msg = " %s : relax" % (POSSIBLE_STATUSES[0])
        logger.info(msg)
        return False
    
    ###########################################################
    outdir = dir_base + "/" + output_directories["relax"]

    ### get initial structure
    filename = outdir + "/full-1/POSCAR.orig"
    
    ### for tar.gz
    struct_orig = _read_poscar(filename, tar=tar)
    if struct_orig is None:
        msg = " %s : full-1" % (POSSIBLE_STATUSES[0])
        logger.info(msg)
        return False

    ### get structure
    num_orig = get_spg_number(struct_orig)
    
    ### list of filenames
    if tar:
        all_files = tar.getnames()
    else:
        all_files = []
    
    ### ver.2
    for dd in ["full-1", "full-2", "freeze-1"]:
        
        dir1 = outdir + "/" + dd

        if os.path.exists(dir1) or dir1 in all_files:
            
            if wasfinished_vasp(dir1, tar=tar) == False:
                msg = " %s : %s" % (POSSIBLE_STATUSES[0], dd)
                logger.info(msg)
                return False
        else:
            msg = " %s : %s" % (POSSIBLE_STATUSES[0], dd)
            logger.info(msg)
            return False
        
        ### symmetry check
        filename = dir1 + "/CONTCAR"
        struct = _read_poscar(filename, tar=tar)
        num_i = get_spg_number(struct)
        if num_orig != num_i:
            msg = " %s : %s" % (POSSIBLE_STATUSES[3], dir1)
            logger.info(msg)
            return False
    
    return True

def get_minimum_energy(dir_name, tar=None):
    """ Check minimum frequency """ 
    from auto_kappa.alamode.log_parser import get_minimum_frequency_from_logfile
    fmin = 1000
    for kind in ["band", "dos"]:
        logfile = dir_name + "/harm/bandos/%s.log" % kind
        try:
            fmin_each = get_minimum_frequency_from_logfile(logfile, tar=tar)
            fmin = min(fmin, fmin_each["minimum_frequency"])
        except Exception:
            return None
    
    return fmin
    
def finish_harmonic(dir_base, tar=None):
    """ Check calculations for harmonic properties """
    
    ### check log files
    dir_bandos = dir_base + "/" + output_directories['harm']['bandos']
    for propt in ["band", "dos"]:
        logfile = dir_bandos + "/" + propt + ".log"
        if _exists(logfile, tar=tar) == False:
            msg = " %s : %s" % (POSSIBLE_STATUSES[0], dir_bandos)
            logger.info(msg)
            return False
        ##
        if wasfinished_alamode(logfile, tar=tar) == False:
            msg = " %s : %s" % (POSSIBLE_STATUSES[0], dir_bandos)
            logger.info(msg)
            return False
    
    return True

def finish_cubic(dir_base, tar=None):
    """ Check calculations for cubic properties """
    
    dir_cube = dir_base + "/cube"
    
    ### check
    if _exists(dir_cube, tar=tar) == False:
        msg = " %s : %s" % (POSSIBLE_STATUSES[0], dir_cube)
        logger.info(msg)
        return False
    
    ### kappa
    if tar is None:
        line = dir_cube + "/kappa_*"
        dirs = glob.glob(line)
    else:
        line = dir_cube + "/kappa_"
        dirs = []
        for name in tar.getnames():
            if (name.startswith(line) and 
                    len(line.split("/")) == len(name.split("/"))):
                dirs.append(name)
    
    ###
    for dd in dirs:
        
        logfile = dd + "/kappa.log"
        
        ### check presence
        if _exists(logfile, tar=tar) == False:
            msg = " %s : %s" % (POSSIBLE_STATUSES[0], dd)
            logger.info(msg)
            return False
        
        ### check if finished or not
        if wasfinished_alamode(logfile, tar=tar) == False:
            msg = " %s : %s" % (POSSIBLE_STATUSES[0], dd)
            logger.info(msg)
            return False
    
    return True

def check_result(dir_base, tar=None):
    """ Check output figures"""
    
    dir_result = dir_base + "/" + output_directories["result"]
    
    names = {
            "harm": ["bandos"],
            "cube": ["cumu_frequency", "cumu_mfp", "kappa", "lifetime", "scat_rates"],
            }
    
    for order in ["harm", "cube"]:
        for name in names[order]:
            figname = dir_result + "/fig_%s.png" % name
            if _exists(figname, tar=tar) == False:
                msg = " %s : %s" % (POSSIBLE_STATUSES[4], dir_result)
                logger.info(msg)
                return False
    return True

def finish_strict_optimization(dir_base, tar=None):
    """ """
    dir_vol = dir_base + "/" + output_directories["relax"] + "/volume"
    for ff in ["POSCAR.init", "POSCAR.opt", "fig_bm.png"]:
        filename = dir_vol + "/" + ff
        if _exists(filename, tar=tar) == False:
            msg = " %s : %s" % (POSSIBLE_STATUSES[0], dir_vol)
            logger.info(msg)
            return False
    return True

def finish_larger_sc(dir_base, tar=None, neg_freq=None):
    """ """
    ### harmonic for the supercell
    flag_harm = finish_harmonic(dir_base, tar=tar)
    if flag_harm == False:
        return [False, "sc_harm"]
    else:
        
        ### check minimum frequency
        fmin = get_minimum_energy(dir_base, tar=tar)
        
        if fmin < neg_freq:
            msg = " %s : %.3f" % (POSSIBLE_STATUSES[1], fmin)
            logger.info(msg)
            return [True, "sc_harm"]
    
    ### cubic for the supercell
    flag_cube = finish_cubic(dir_base, tar=tar)
    if flag_cube == False:
        return [False, "sc_cube"]
    
    ### check result
    flag = check_result(dir_base, tar=tar)
    if flag == False:
        return [False, "sc_result"]
    
    return [True, "sc_cube"]

def finish_larger_supercells(dir_base, tar=None, neg_freq=None):
    """ check analysis with larger supercell(s) """    
    
    ### get names of supercells
    if tar is None:
        line = dir_base + "/sc-*"
        dirs = glob.glob(line)
    else:
        line = dir_base + "/sc-"
        dirs = []
        for name in tar.getnames():
            if (name.startswith(line) and 
                    len(line.split("/")) == len(name.split("/"))):
                dirs.append(name)
    
    ### check for each supercell
    for dir_sc in dirs:
        out = finish_larger_sc(dir_sc, tar=tar, neg_freq=neg_freq)
        if out[0] == False:
            msg = " %s : %s" % (POSSIBLE_STATUSES[0], dir_sc)
            logger.info(msg)
            return out
        else:
            if out[1] == "sc_harm":
                return out
    
    return [True, "sc"]

def _get_tar_prefix(tar):
    """ """
    data = tar.name.replace(os.getcwd() + "/", "").split("/")
    prefix = ""
    for i in range(len(data)-1):
        prefix += data[i] + "/"
    return prefix

def _get_tar_relative_names(tar):
    """ """
    prefix = _get_tar_prefix(tar)
    all_files = [prefix + name for name in tar.getnames()]
    return all_files

def _exists(filename, tar=None):
    """ """
    if tar is None:
        if os.path.exists(filename) == False:
            return False
    else:
        if filename not in tar.getnames():
            return False
    return True 

def finish_ak_calculation(dir_tmp, neg_freq=-0.001, vol_relax=None, larger_sc=None):
    """ Check output directory generated by auto-kappa """
    
    suffix = ".tar.gz"
    if dir_tmp.endswith(suffix):
        
        ### for .tar.gz file
        tar = tarfile.open(dir_tmp, 'r')
        dir_type = "tar.gz"
        dir_base = dir_tmp.split("/")[-1].replace(suffix, "")
    else:
        tar = None
        dir_type = "directory"
        dir_base = dir_tmp
    
    ### log file
    logfile = dir_base + "/ak.log"
    if _exists(logfile, tar=tar) == False:
        msg = " %s : %s" % (POSSIBLE_STATUSES[0], dir_base)
        logger.info(msg)
        return [False, "ak.log"]
    
    ### check structure optimization
    flag_relax = finish_relaxation(dir_base, tar=tar)
    if flag_relax == False:
        return [False, "relax"]
    
    ### check strict optimization
    if vol_relax == 1:
        flag_vol = finish_strict_optimization(dir_base, tar=tar)
        if flag_vol == False:
            return [False, "strict"]
    
    ### check NAC
    dir_nac = dir_base + "/" + output_directories["nac"]
    if _exists(dir_nac, tar=tar):
        if wasfinished_vasp(dir_nac, tar=tar) == False:
            msg = " %s : %s" % (POSSIBLE_STATUSES[0], dir_nac)
            logger.info(msg)
            return [False, "dir_nac"]
    
    ### check harmonic
    flag_harm = finish_harmonic(dir_base, tar=tar)
    
    if flag_harm == False:
        return [False, "harm"]
    else:
        
        ### check minimum frequency
        fmin = get_minimum_energy(dir_base, tar=tar)
        
        if fmin < neg_freq:
            
            ### If negative frequencies were round,
            if larger_sc == 0:
                msg = " %s : %.3f" % (POSSIBLE_STATUSES[1], fmin)
                logger.info(msg)
                return [True, "harm"]
            else:
                ### calculation with supercell
                flag = finish_larger_supercells(
                        dir_base, tar=tar, neg_freq=neg_freq)
                return flag
     
    ### check cubic
    flag_cube = finish_cubic(dir_base, tar=tar)
    if flag_cube == False:
        return [False, "cube"]
    
    ### check result
    flag = check_result(dir_base, tar=tar)
    if flag == False:
        return [False, "result"]

    return [True, "cube"]

def main(options):
    
    flag = finish_ak_calculation(
            options.directory, 
            neg_freq=options.neg_freq,
            vol_relax=options.vol_relax,
            larger_sc=options.larger_sc)
    
    if flag[0]:
        if "cube" in flag[1].lower():
            msg = " %s" % (POSSIBLE_STATUSES[2])
            logger.info(msg)
        elif flag[1].lower() in ["sc"]:
            msg = " %s" % (POSSIBLE_STATUSES["sc"])
            logger.info(msg)
    
if __name__ == '__main__':
    parser = OptionParser()
    
    parser.add_option("-d", "--directory", dest="directory", type="string",
                help="directory name for APDB")
    
    parser.add_option("--neg_freq", dest="neg_freq", type="float",
                default=-0.001, help="negtive frequency [-0.001]")
    
    parser.add_option("--vol_relax", dest="vol_relax", type="int",
                default=0, help="check volume relaxation (0.No, 1.Yes) [0]")
    
    parser.add_option("--larger_sc", dest="larger_sc", type="int",
                default=0, help="check larger supercell (0.No, 1.Yes) [0]")
    
    (options, args) = parser.parse_args()
    
    if options.directory.startswith("./"):
        options.directory = options.directory.replace("./", "")
    
    main(options)

