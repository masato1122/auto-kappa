# -*- coding: utf-8 -*-
import math
import numpy as np
import pandas as pd

from auto_kappa.alamode.log_parser import (
        get_eigenvalues_from_logfile,
        get_minimum_frequency_from_logfile)

def get_representative_kpoint_with_negative_frequency(
        logfile_band=None, logfile_dos=None, 
        logfile_commensurate=None, tol_k=1e-5):
    """ Get representative negative eigen mode.
    1. If w < 0 exist at k = 0, return k = 0.
    2. If w < 0 exist at commensurate points, return the corresponding kpoint.
    3. If w < 0 is found only in DOS, return the corresponding kpoint.
    
    Return
    -------
    integer : 
        0. Gamma point
        1. kpoint on phonon dispersion
        2. kpoint of the commensurate points
        3. kpoint on k-mesh for DOS
    """
    def _get_kpoint_of_fmin(kpoints, frequencies):
        """ Get kpoint of the minimum frequency
        Args
        =======
        kpoints : shape=(nk, 3)
        frequencies : shape=(nk, nbands)
        """
        ### check each kpoint
        ks_neg = []
        fs_neg = []
        for ik in range(len(kpoints)):
            
            keach = kpoints[ik]
            feach = frequencies[ik]
            
            ### check the minimum frequency
            fmin = np.min(feach)
            if fmin < 0.:
                if np.linalg.norm(keach) < tol_k:
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
    
    ks1, es1 = get_eigenvalues_from_logfile(logfile_band)
    ks2, es2 = get_eigenvalues_from_logfile(logfile_commensurate)
    ks3, es3 = get_eigenvalues_from_logfile(logfile_dos)
    
    imod = None
    kmin = None
    fmin = 10000.
    
    ### phonon band
    kfmin1 = _get_kpoint_of_fmin(ks1, es1)
    if kfmin1 is not None:
        
        if kfmin1[1] < fmin:
            imod = 1
            kmin = kfmin1[0]
            fmin = kfmin1[1]
        
        ### case 1: w < 0 at k = 0
        if np.linalg.norm(kfmin1[0]) < tol_k:
            return 0, kfmin1[0]
    
    ### commensurate points
    kfmin2 = _get_kpoint_of_fmin(ks2, es2)
    if kfmin2 is not None:
        
        if kfmin2[1] < fmin:
            imod = 2
            kmin = kfmin2[0]
            fmin = kfmin2[1]
        
        ### case 2: w < 0 at commensurate points
        return 1, kfmin2[0]
    
    ### DOS
    kfmin3 = _get_kpoint_of_fmin(ks3, es3)
    if kfmin3 is not None:

        if kfmin3[1] < fmin:
            imod = 3
            kmin = kfmin2[0]
            fmin = kfmin2[1]
    
    ### case 3: return kpoint with the min. freq.
    return imod, kmin

def get_symmetry_points_for_kpoint(kpoint, filename=None, tol=1e-5):
    """ Search the given kpoint in the band file (*.bands) and return the
    symmetry points which the given kpoint is between. """

    from auto_kappa.alamode.log_parser import get_kpath
    kpaths = get_kpath(filename) 
    
    for ip, kline in enumerate(kpaths):
        
        keys = list(kline.keys())
        k0 = kline[keys[0]]
        k1 = kline[keys[1]]
        
        frac = (k1 - k0) / (kpoint - k0)
        
        if (abs(frac[0] - frac[1]) < 1e-5 and
                abs(frac[0] - frac[2]) < 1e-5):
            return keys
    
    return None


log_band = "./bandos/band.log"
log_dos = "./bandos/dos.log"
log_com = "./evec/evec_commensurate.log"

itype, kmin = get_representative_kpoint_with_negative_frequency(
        logfile_band=log_band, 
        logfile_dos=log_dos, 
        logfile_commensurate=log_com)

print(itype, kmin)

sym_points = get_symmetry_points_for_kpoint(kmin, filename=log_band)
print(sym_points)

# test

