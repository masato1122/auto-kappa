# -*- coding: utf-8 -*-
import sys
import os, os.path
import math
import numpy as np
import pandas as pd
import glob

from auto_kappa.alamode.log_parser import (
        get_eigenvalues_from_logfile,
        get_minimum_frequency_from_logfile)

import logging
logger = logging.getLogger(__name__)

def calculate_pes(
        primitive, kpoint, 
        outdir="./pes", fcsxml="FC2xml", 
        command={"mpirun": "mpirun", "nprocs": 2, "anphon": "anphon"},
        nac=0, vasp_xml=None):
    """ Calculate potential energy surface """
    
    ### calculate eigenvector
    dir_evec = outdir + "/evec"
    calculate_evec(
            primitive, kpoint, 
            outdir=dir_evec, fcsxml=fcsxml,
            nac=nac, vasp_xml=vasp_xml, command=command)
     
    ### make vibration

    exit()
    return None

def calculate_evec(
        prim, kpoint, outdir="./evec", fcsxml=None, command=None,
        vasp_xml=None, file_anphon="evec.in", nac=0):
    """ Calculate eiven vector at the given kpoint """
    
    from auto_kappa.io.vasp import write_born_info
    from auto_kappa.io.alm import AnphonInput
    
    ### move to the target directory
    cwd = os.getcwd()
    os.makedirs(outdir, exist_ok=True)
    os.chdir(outdir)

    ### vasprun.xml for Born effective charge
    if vasp_xml is None:
        nac = 0
    else:
        if os.path.exists(vasp_xml) == False:
            nac = 0
            vasp_xml = None
    
    ### make BORNINFO
    borninfo = "BORNINFO"
    if vasp_xml is not None:
        write_born_info(vasp_xml, outfile=borninfo)
    
    ### make AnphonInput obj
    inp = AnphonInput.from_structure(
            prim, mode="phonons", kpmode=0,
            fcsxml=os.path.relpath(fcsxml),
            nonanalytic=nac,
            borninfo=borninfo)
    inp.update({"printevec": 1})
    inp.set_kpoint(kpoints=[kpoint])
    inp.to_file(filename=file_anphon)
    
    ### run anphon
    from auto_kappa.alamode.almcalc import run_alamode
    logfile = "evec.log"
    status = run_alamode(
            file_anphon, logfile, workdir=".", 
            mpirun=command["mpirun"],
            nthreads=command["nprocs"],
            command=command["anphon"])
    
    ### back to the initial directory
    os.chdir(cwd)
    
    return status

def get_representative_kpoint(
        log_band=None, log_dos=None, log_commensurate=None, 
        negative_freq=None, epsilon=None):
    """ Get representative kpoint. If negative frequencies cannot be found,
    return Gamma point as a representative kpoint. """

    ###
    for log in [log_band, log_dos, log_commensurate]:
        if os.path.exists(log) == False:
            msg = "\n Warning: cannot find %s" % log
            logger.warning(msg)
            return None
    
    ### get representative kpoint
    out = get_representative_kpoint_with_negative_freq(
            log_band=log_band,
            log_dos=log_dos,
            log_commensurate=log_commensurate)
    
    if out[0] is None:
        kp_type = "band"
        kpoint = np.zeros(3)
    else:
        kp_type = out[0]
        kpoint = out[1]
    
    ### print info
    msg = "\n Representative k-point : " + "%.3f " * 3 % tuple(kpoint)
    if np.linalg.norm(kpoint) > epsilon:
        if kp_type == "band":
            try:
                sym_points = get_symmetry_points_for_kpoint(
                        kpoint, filename=log_band, tol=epsilon)
                if sym_points is not None:
                    if isinstance(sym_points, str):
                        msg += "\n k-point is %s" % sym_points
                    else:
                        msg += "\n k-point is located between %s and %s (%.2f)" % (
                                sym_points[0], sym_points[1], sym_points[2])
            except Exception:
                pass
        elif kp_type == "commensurate":
            msg += "\n k-point is a commensurate point."
        elif kp_type == "dos":
            msg += "\n k-point is one of the grid-points for DOS."
    
    logger.info(msg)
    return [kp_type, kpoint]

def get_representative_kpoint_with_negative_freq(
        log_band=None, log_dos=None, 
        log_commensurate=None, epsilon=1e-5,
        negative_freq=-1e-3):
    """ Get representative negative eigen mode.
    1. If w < 0 exist at k = 0, return k = 0.
    2. If w < 0 exist at commensurate points, return the corresponding kpoint.
    3. If w < 0 is found only in DOS, return the corresponding kpoint.
    
    Return
    -------
    string :
        "band"        : kpoint on phonon dispersion
        "commensurate": kpoint of the commensurate points
        "dos"         : kpoint on k-mesh for DOS

    array : 
        kpoint for a negative frequency
    """
    def _get_kpoint_of_fmin(kpoints, frequencies):
        """ Get kpoint of the minimum frequency
        Args
        -------
        kpoints : shape=(nk, 3)
        
        frequencies : shape=(nk, nbands)
        
        Return
        ---------
        out : array with two components or None
        
        out[0] : array, shape=(3)
            k-point

        out[1] : float
            minimum frequency
        """
        ### check each kpoint
        ks_neg = []
        fs_neg = []
        for ik in range(len(kpoints)):
            
            keach = kpoints[ik]
            feach = frequencies[ik]
            
            ### check the minimum frequency
            fmin = np.min(feach)
            if fmin < negative_freq:
                if np.linalg.norm(keach) < epsilon:
                    return [np.zeros(3), np.min(feach)]
                else:
                    ks_neg.append(keach)
                    fs_neg.append(fmin)
        ##
        if len(ks_neg) == 0:
            return None
        else:
            imin = np.argmin(fs_neg)
            return [ks_neg[imin], fs_neg[imin]]
    
    ks1, es1 = get_eigenvalues_from_logfile(log_band)
    ks2, es2 = get_eigenvalues_from_logfile(log_commensurate)
    ks3, es3 = get_eigenvalues_from_logfile(log_dos)
    
    kp_type = None
    kmin = None
    fmin = 10000.
    
    ### phonon band
    kfmin1 = _get_kpoint_of_fmin(ks1, es1)
    if kfmin1 is not None:
        
        if kfmin1[1] < fmin:
            kp_type = "band"
            kmin = kfmin1[0]
            fmin = kfmin1[1]
        
        ### case 1: w < 0 at k = 0
        if np.linalg.norm(kfmin1[0]) < epsilon:
            return "band", kfmin1[0]
    
    ### commensurate points
    kfmin2 = _get_kpoint_of_fmin(ks2, es2)
    if kfmin2 is not None:
        
        if kfmin2[1] < fmin:
            kp_type = "commensurate"
            kmin = kfmin2[0]
            fmin = kfmin2[1]
        
        ### case 2: w < 0 at commensurate points
        return "commensurate", kfmin2[0]
    
    ### DOS
    kfmin3 = _get_kpoint_of_fmin(ks3, es3)
    if kfmin3 is not None:

        if kfmin3[1] < fmin:
            kp_type = "dos"
            kmin = kfmin3[0]
            fmin = kfmin3[1]
    
    ### case 3: return kpoint with the min. freq.
    return kp_type, kmin

def get_symmetry_points_for_kpoint(kpoint, filename=None, tol=1e-5):
    """ Search the given kpoint in the log-file of band (band.log) and return 
    symmetry points which the given kpoint is between. """
    
    from auto_kappa.alamode.log_parser import get_kpath
    kpaths = get_kpath(filename) 
    
    for ip, kline in enumerate(kpaths):
        
        keys = list(kline.keys())
        k0 = kline[keys[0]]
        k1 = kline[keys[1]]
        
        if np.linalg.norm(k0 - kpoint) < tol:
            return keys[0]
        if np.linalg.norm(k1 - kpoint) < tol:
            return keys[1]
        
        ### check distances 
        d0k = np.linalg.norm(kpoint - k0)
        d1k = np.linalg.norm(kpoint - k1)
        d01 = np.linalg.norm(k1 - k0)
        if abs(d01 - d0k - d1k) < tol:
            frac = d0k / d01
            return [keys[0], keys[1], frac]
    
    return None


#log_band = "./bandos/band.log"
#log_dos = "./bandos/dos.log"
#log_com = "./evec/evec_commensurate.log"
#
#kp_type, kmin = get_representative_kpoint_with_negative_frequency(
#        log_band=log_band, 
#        log_dos=log_dos, 
#        log_commensurate=log_com)
#
#print(kp_type, kmin)
#
#if kp_type == "band":
#    sym_points = get_symmetry_points_for_kpoint(kmin, filename=log_band)
#    print(sym_points)

