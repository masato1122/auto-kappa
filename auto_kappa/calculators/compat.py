#
# compat.py
#
# Ensures compatibility
#
# Copyright (c) 2025 Masato Ohnishi
#
# This file is distributed under the terms of the MIT license.
# Please see the file 'LICENCE.txt' in the root directory
# or http://opensource.org/licenses/mit-license.php for information.
#
import sys
import os
import glob
import fnmatch
import shutil
import subprocess

from auto_kappa import output_files

import logging
logger = logging.getLogger(__name__)

def remove_old_kappa_data(directories):
    
    ## Remove kappa files in "result" directory
    dir_result = directories['result']
    if os.path.exists(dir_result):
        _check_result_kappa(dir_result)
    
    ## Remove "kappa" directories
    for key in directories['cube']:
        if key.startswith('kappa'):
            dir_kappa = directories['cube'][key]
            _check_kappa_directory(dir_kappa)


def _check_result_kappa(dir_result):
    """ Check files in "result" directory 
    """
    fns = glob.glob(dir_result + '/*')
    for path_abs in fns:
        fn = path_abs.split('/')[-1]
        
        remove = False
        
        ## csv files
        for kind in ['tau', 'kspec']:
            if fnmatch.fnmatch(fn, f"{kind}_*.csv"):
                remove = True
                
        ## figures
        for kind in ['kappa', 'cumu', 'lifetime', 'scat_rates', 'kappa']:
            if fnmatch.fnmatch(fn, f"fig_{kind}*.png"):
                remove = True
        
        if remove:
            
            relative_path = os.path.relpath(path_abs, os.getcwd())
            os.remove(relative_path)
            
            msg = f' Remove {relative_path}'
            logger.info(msg)
    

def _check_kappa_directory(dir_kappa):
    """ Check "kappa" directories
    """
    line = dir_kappa + "_*"
    dirs = glob.glob(line)
    if len(dirs) == 0:
        return
    
    ## Make "box" directory
    dir_box = os.path.dirname(dirs[0]) + "/box"
    os.makedirs(dir_box, exist_ok=True)
    
    ## Move "kappa_*" directories
    logger.info("")
    for dir_each in dirs:
        
        relative_path = os.path.relpath(dir_each, os.getcwd())
        name = os.path.basename(dir_each)
        
        ## Set the distination directory name
        count = 1
        while True:
            dir_dist = f"{dir_box}/{name}_{count}"
            if not os.path.exists(dir_dist):
                break
            count += 1
        
        cmd = f"mv {relative_path} {dir_dist}"
        subprocess.run(cmd, shell=True)
        
        msg = f' Remove {relative_path}'
        logger.info(msg)
