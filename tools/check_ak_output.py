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
#from auto_kappa.alamode.io import wasfinished_alamode
from auto_kappa.alamode.io import get_status

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

def file2str(filename, tar=None):
    if tar is None:
        if os.path.exists(filename) == False:
            return None
        else:
            lines = open(filename, 'r').readlines()
            return lines
    else:
        try:
            lines_tmp = tar.extractfile(filename).readlines()
            lines = [ll.decode('utf-8').replace("\n", "") for ll in lines_tmp]
            return lines
        except Exception:
            return None

def _read_poscar(filename, tar=None):
    """ Read POSCAR file to get a structure object
    
    Args
    -----
    filename : string
    
    tar : tarfile.TarFile
    
    """
    if tar is None:
        try:
            return ase.io.read(filename, format='vasp')
        except Exception:
            return None
    else:
        try:
            from pymatgen.io.vasp import Poscar
            content = tar.extractfile(filename).read().decode('utf-8')
            poscar = Poscar.from_str(content)
            return poscar.structure
        except Exception:
            #print(" Error: %s" % filename)
            return None
    return None

def finish_relaxation(dir_base, tar=None):
    """ 
    Args
    -----
    dir_base : string
        Name of a directory or .tar.gz file
    
    tar : tarfile.TarFile
    """
    #### easy check.
    dir_harm = dir_base + "/harm"
    if _exists(dir_harm, tar=tar):
        return True
    else:
        msg = " %s : relax" % (POSSIBLE_STATUSES[0])
        logger.info(msg)
        return False
    
    ###########################################################
    outdir = dir_base + "/" + output_directories["relax"]
    struct_orig = None
    for j in range(2):
        
        ### get initial structure
        if j == 0:
            filename = outdir + "/full-1/POSCAR.orig"
        else:
            filename = outdir + "/POSCAR.orig"
        
        ### for tar.gz
        try:
            struct_orig = _read_poscar(filename, tar=tar)
        except Exception:
            pass
        
        if struct_orig is not None:
            break
    
    if struct_orig is None:
        msg = " %s : relax" % (POSSIBLE_STATUSES[0])
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
    
def get_status_harmonic(dir_base, tar=None):
    """ Check calculations for harmonic properties """
    
    ### check log files
    dir_bandos = dir_base + "/" + output_directories['harm']['bandos']
    statuses = {}
    for propt in ["band", "dos"]:
        logfile = dir_bandos + "/" + propt + ".log"
        statuses[propt] = get_status(logfile, tar=tar)
        
    if (statuses["band"][0].lower() == "finished" and
            statuses["dos"][0].lower() == "finished"):
        status = "Finished"
        comment = ""
    
    elif (statuses["band"][0].lower() == "error" or
            statuses["dos"][0].lower() == "error"):
        status = "Error"
        comment = "%s\n%s" % (statuses["band"][1], statuses["dos"][1])
        print(dir_base, comment)
    
    else:
        status = "NotYet"
        comment = "%s\n%s" % (statuses["band"][1], statuses["dos"][1])
    
    return [status, comment]

def get_status_cubic(dir_base, tar=None):
    """ Check calculations for cubic properties and return its status """
    
    dir_cube = dir_base + "/cube"
    
    ### check
    if _exists(dir_cube, tar=tar) == False:
        msg = " %s : %s" % (POSSIBLE_STATUSES[0], dir_cube)
        return ["NotYet", msg]
    
    logfiles = []
    
    ### FC3
    log_fc3 = dir_cube + "/force_fd/fc3.log"
    log_lasso = dir_cube + "/lasso/lasso.log"
    if tar is None:
        for ll in [log_fc3, log_lasso]:
            if os.path.exists(ll):
                logfiles.append(ll)
    else:
        for ll in [log_fc3, log_lasso]:
            if ll in tar.getnames():
                logfiles.append(ll)
    
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
    
    for dd in dirs:
        fn = dd + "/kappa.log"
        if os.path.exists(fn):
            logfiles.append(fn)
    
    ###
    flag_fin = True
    flag_error = False
    msg = ""
    for i, logfile in enumerate(logfiles):
        
        ### get status
        status = get_status(logfile, tar=tar)
        
        if status[0].lower() == "error":
            #print(dir_base, status[1])
            flag_error = True
            msg += status[1]
            if i <= len(logfiles):
                msg += ", "

        if status[0].lower() != "finished":
            flag_fin = False
            msg += status[1]
            if i <= len(logfiles):
                msg += ", "
        
        ### check kl and kl_coherent files
        if logfile in ["kappa.log"]:
            for suffix in [".kl", ".kl_coherent"]:
                line = os.path.dirname(logfile) + "/*" + suffix
                fns = glob.glob(line)
                
                _is_error = False
                if len(fns) == 0:
                    _is_error = True
                else:
                    try:
                        dump = np.genfromtxt(fns[0])
                        if np.isnan(dump).any():
                            _is_error = True
                    except Exception:
                        _is_error = True
                
                if _is_error:
                    flag_error = True
                    msg += "Error while reading %s" % (suffix.replace(".", ""))
    
    if flag_error:
        return ["error", msg]
    else:
        if flag_fin:
            return ["finished", msg]
        else:
            return ["NotYet", msg]

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
    status_harm = get_status_harmonic(dir_base, tar=tar)
    
    if status_harm[0].lower() == "finished":
        ### check minimum frequency
        fmin = get_minimum_energy(dir_base, tar=tar)
        if fmin < neg_freq:
            msg = " %s : %.3f" % (POSSIBLE_STATUSES[1], fmin)
            logger.info(msg)
            return [True, "sc_harm"]
    
    else:
        if status_harm[0].lower() == "error":
            return [True, "sc_harm_error"]
        else:
            return [False, "sc_harm"]
    
    ### >>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    ### cubic for the supercell
    status_cube = get_status_cubic(dir_base, tar=tar)
    
    if status_cube[0].lower() == "error":
        return [True, "sc_cube_error"]
    elif status_cube[0].lower() == "notyet":
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
    
    if len(dirs) == 0:
        return [False, "sc_NotStarted"]
    
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

def _memory_over(logfile, tar=None):
    """ """
    lines = file2str(logfile, tar=tar)
    
    flag_single = False
    flag_stop = False
    for ll in lines:
        if "processes per node" in ll.lower() and "=> 1" in ll.lower():
            flag_single = True
        if "Error : ALAMODE job was not finished properly":
            flag_stop = True
    ##
    if flag_single and flag_stop:
        return True
    else:
        return False

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
    else:
        if _memory_over(logfile, tar=tar):
            return [True, "memory"]
    
    ### check structure optimization
    if finish_relaxation(dir_base, tar=tar) == False:
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
    status_harm = get_status_harmonic(dir_base, tar=tar)
    
    if status_harm[0].lower() == "finished":    
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
                out = finish_larger_supercells(
                        dir_base, tar=tar, neg_freq=neg_freq)
                return out
    
    elif status_harm[0].lower() == "error":
            return [True, "harm_error"]
    elif status_harm[0].lower() == "notyet":
            return [False, "harm"]
    
    ### check cubic
    status_cube = get_status_cubic(dir_base, tar=tar)
    
    if status_cube[0].lower() == "error":
        return [True, "cube_error"]
    elif status_cube[0].lower() == "notyet":
        return [False, "cube"]
    
    ### check result
    flag = check_result(dir_base, tar=tar)
    if flag == False:
        return [False, "result"]
    
    return [True, "cube"]

#def main(options):
#    
#    flag = finish_ak_calculation(
#            options.directory, 
#            neg_freq=options.neg_freq,
#            vol_relax=options.vol_relax,
#            larger_sc=options.larger_sc)
#    
#    if flag[0]:
#        if "cube" in flag[1].lower():
#            msg = " %s" % (POSSIBLE_STATUSES[2])
#            logger.info(msg)
#        elif flag[1].lower() in ["sc"]:
#            msg = " %s" % (POSSIBLE_STATUSES["sc"])
#            logger.info(msg)
#    
#if __name__ == '__main__':
#    parser = OptionParser()
#    
#    parser.add_option("-d", "--directory", dest="directory", type="string",
#                help="directory name for APDB")
#    
#    parser.add_option("--neg_freq", dest="neg_freq", type="float",
#                default=-0.001, help="negtive frequency [-0.001]")
#    
#    parser.add_option("--vol_relax", dest="vol_relax", type="int",
#                default=0, help="check volume relaxation (0.No, 1.Yes) [0]")
#    
#    parser.add_option("--larger_sc", dest="larger_sc", type="int",
#                default=0, help="check larger supercell (0.No, 1.Yes) [0]")
#    
#    (options, args) = parser.parse_args()
#    
#    if options.directory.startswith("./"):
#        options.directory = options.directory.replace("./", "")
#    
#    main(options)

