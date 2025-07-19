#
# scph.py
#
# This script helps to conduct SCPH calculations.
#
# Copyright (c) 2024 Masato Ohnishi
#
# This file is distributed under the terms of the MIT license.
# Please see the file 'LICENCE.txt' in the root directory
# or http://opensource.org/licenses/mit-license.php for information.
#
# import sys
import os
import os.path
import numpy as np
import subprocess

import logging
logger = logging.getLogger(__name__)

def calculate_high_order_force_constants(
    almcalc, calculator, order=5, frac_nrandom=None, disp_temp=500,):
    """ Calculate high-order (up to 6th-order) force constants finally to
    calculate 4th order FCs for phonon renormalization.
    """
    ### calculate forces for high-order FCs
    almcalc.write_alamode_input(propt='suggest', order=order)
    almcalc.run_alamode(propt='suggest', order=order)
    almcalc.calc_forces(
            order=order, calculator=calculator,
            frac_nrandom=frac_nrandom,
            temperature=disp_temp,
            output_dfset=2,
            )
    
    ### calculate anharmonic force constants
    for propt in ['cv', 'lasso']:
        almcalc.write_alamode_input(propt=propt, order=order)
        almcalc.run_alamode(propt, order=order, neglect_log=False)
    
def conduct_scph_calculation(almcalc, order=6, temperatures=100*np.arange(1,11)):
    """ Conduct SCPH calculation and obtain effective harmonic FCs 
    """
    msg =  "\n Conduct SCPH calculation"
    msg += "\n ========================"
    logger.info(msg)
    
    ### parameters for SCPH
    tmin = np.min(temperatures)
    tmax = np.max(temperatures)
    dt = temperatures[1] - temperatures[0]
    
    ### SCPH calculation
    propt = "scph"
    almcalc.write_alamode_input(propt=propt, tmin=tmin, tmax=tmax, dt=dt)
    almcalc.run_alamode(propt, order=order, neglect_log=False)
    
    ### Create effectvie harmonic FCs
    for temp in temperatures:
        _create_effective_harmonic_fcs(almcalc, temperature=temp)

def _create_effective_harmonic_fcs(
        almcalc, temperature=300, workdir=None):
    """ Create effective harmonic FCs """
    
    prefix = almcalc.prefix
    command = almcalc.commands['alamode']['dfc2']
    
    ### get file names
    xml_orig = os.path.relpath(
            almcalc.out_dirs["higher"]["lasso"] + "/%s.xml" % prefix,
            almcalc.out_dirs["higher"]["scph"])
    xml_new = "%s_%dK.xml" % (prefix, temperature)
    file_corr = "%s.scph_dfc2" % prefix
    
    ### change directory
    dir_cur = os.getcwd()
    if workdir is None:
        workdir = almcalc.out_dirs['higher']['scph']
    os.chdir(workdir)
    
    ### run the job
    cmd = "%s %s %s %s %d" % (command, xml_orig, xml_new, file_corr, temperature)
    
    logfile = "dfc2.log"
    file_err = "std_err.txt"
    with open(logfile, 'w') as f, open(file_err, "w", buffering=1) as f_err:
        proc = subprocess.Popen(
                cmd, shell=True, env=os.environ, stdout=f, stderr=f_err)
    
    ### back to the original directory
    os.chdir(dir_cur)
    return 0

def set_parameters_scph(
        inp, primitive=None, deltak=0.01, kdensities=[30, 10], **kwargs):
    """ Set ALAMODE parameters for SCPH calculation.
    
    Args
    =====
    
    inp : auto_kappa.io.AnphonInput
    
    deltak : integer, 0.01
    
    kdensity : array of float,
        kdensity for SCPH calculation and interpolated kmesh
    
    """
    from auto_kappa.structure.crystal import get_automatic_kmesh
    
    scph_params = {
            ### general
            "tmin": 100,
            "tmax": 1000,
            "dt"  : 100,
            ### scph
            #"kmesh_scph": [1,1,1],
            "kmesh_interpolate": [1,1,1],
            "self_offdiag": 1,    ## not default (0)
            "mixalpha": 0.1,      ## default
            "maxiter": 2000,      ## double of default
            "tol_scph": 1e-10,    ## default
            }
    
    inp.set_kpoint(deltak=deltak)
    
    ### kmesh_scph
    kmesh_scph = get_automatic_kmesh(
            primitive, reciprocal_density=kdensities[0], dim=inp.dim)
    
    ### kmesh_interpolate
    kmesh_int = get_automatic_kmesh(
            primitive, reciprocal_density=kdensities[1], dim=inp.dim)
    
    ### kmesh_scph should be equal to or a multiple of the number of
    ### kmesh_interpolate in the same direction.
    ratio = np.asarray(kmesh_scph) / np.asarray(kmesh_int)
    ratio = (ratio + 0.5*np.ones(3)).astype(int)
    ratio = np.where(ratio < 1, 1, ratio)
    kmesh_scph = list(np.asarray(kmesh_int) * ratio)
    
    ### set parameters
    scph_params["kmesh_scph"] = kmesh_scph
    scph_params["kmesh_interpolate"] = kmesh_int
    
    ###
    scph_params.update(kwargs)
    
    ### update!
    inp.update(scph_params) 
    
