#
# utils.py
#
# Interface for ALAMODE
#
# Copyright (c) 2025 Masato Ohnishi
#
# This file is distributed under the terms of the MIT license.
# Please see the file 'LICENCE.txt' in the root directory
# or http://opensource.org/licenses/mit-license.php for information.
#
import os
import random
import subprocess
import shutil

from auto_kappa.io import AlmInput
from auto_kappa.units import AToBohr

def get_num_disp_patterns(supercell, norder=1, tolerance=1e-5, command_alm='alm'):
    """ Get the number of displacement patterns using ALAMODE.
    
    Args
    ------
    supercell : ase Atoms object
        Supercell of the crystal structure
    
    norder : int
        Order of the force constants
    
    tolerance : float
        Tolerance for the symmetry operations with the unit of Angstrom
    """
    return None

    # # Make an AlmInput object
    # alminp = AlmInput.from_structure(
    #     supercell,
    #     mode='suggest',
    #     norder=norder,
    #     tolerance=tolerance * AToBohr
    #     )
    
    # ## Move to a temporary directory
    # cwd = os.getcwd()
    # workdir = os.path.join(cwd, "tmp_%d" % random.randint(0, 1e9))
    # os.makedirs(workdir, exist_ok=True)
    # os.chdir(workdir)
    
    # filename = "suggest.in"
    # logfile = "suggest.log"
    # alminp.to_file(filename=filename)
    # cmd = f"{command_alm} {filename} > {logfile}"
    # try:
    #     subprocess.run(cmd, shell=True)
    # except:
    #     print("Error: ALAMODE failed to run.")
    #     os.chdir(cwd)
    #     shutil.rmtree(workdir)
    #     return None
    
    # ##
    # lines = open(logfile).readlines()
    # num_patterns = None
    # for line in lines:
    #     if "Number of disp. patterns for" in line:
    #         n = int(line.split()[-1])
    #         if norder == 1:
    #             if "HARMONIC" in line:
    #                 num_patterns = n
    #         else:
    #             word = "ANHARM%d" % (norder + 1)
    #             if word in line:
    #                 num_patterns = n
    
    # ## Go back to the original directory
    # os.chdir(cwd)
    
    # ## Remove the temporary directory
    # shutil.rmtree(workdir)
    
    # return num_patterns
