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
    
def check_log_yaml(dir_name, tol_zero=-1e-3):
    """
    Return
    --------
    0 : not yet finished
    1 : finished at harmonic
    2 : finished at cubic
    """
    
    emin = get_minimum_energy(dir_name)
    if emin is None:
        return 0
    
    ##
    if emin < tol_zero:
        return 1
    else:
        figname = dir_name + "/result/fig_kappa.png"
        if os.path.exists(figname):
            return 2
        else:
            return 0

def main(options):

    flag = check_log_yaml(options.dir_apdb)
    
    ##print(options.dir_apdb, end=" ")
    if flag == 0:
        print("NotYet")
    elif flag == 1:
        print("Finished_harmonic")
    elif flag == 2:
        print("Finished_cubic")
    else:
        print("Error")
    
    
if __name__ == '__main__':
    parser = OptionParser()
    
    parser.add_option("-d", "--dir_apdb", dest="dir_apdb", type="string",
            help="directory name for APDB")
    
    (options, args) = parser.parse_args()
    main(options)

