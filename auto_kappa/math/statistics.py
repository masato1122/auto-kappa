#
# statistics.py
#
# This file calculates Bose-Einstein and Fermi-Dirac distributions.
#
# Copyright (c) 2022 Masato Ohnishi
#
# This file is distributed under the terms of the MIT license.
# Please see the file 'LICENCE.txt' in the root directory
# or http://opensource.org/licenses/mit-license.php for information.
#
import sys
import numpy as np
import auto_kappa.units as unit

import logging
logger = logging.getLogger(__name__)

EPS_DEFAULT = 1e-7

def get_statistical(ftype, temperature, ene_tmp, diff=0, eps=EPS_DEFAULT):
    """
    temperature : unit [K]
    energy : unit [J]
    """
    energy = np.zeros_like(ene_tmp)
    energy = ene_tmp.copy()
    energy /= unit.kb            # J/kb
    
    if ftype.lower() == "be":
        X = np.where(abs(temperature)<eps, 0.0, energy/temperature)
        X = np.where(abs(X)<eps, eps, X)
        if diff == 0:
            return np.where(abs(temperature)<eps or X<eps, 0.,
                    1/(np.exp(X)-1.0))
    elif ftype.lower() == "fd":
        X = np.where(abs(temperature)<eps, eps, energy/temperature)
        if diff == 0:
            return np.where(abs(temperature)<eps, 0., 1/(np.exp(X)+1.0))
    
    logger.getLogger(" ERROR")
    sys.exit()

def get_bose(temperature, frequency, diff=0, eps=EPS_DEFAULT):
    return get_statistical("BE", temperature, frequency, diff=diff, eps=eps)

def get_fd(temperature, energy, diff=0, eps=EPS_DEFAULT):
    return get_statistical("FD", temperature, energy, diff=diff, eps=eps)

def get_diffrential_statistics(energy, Ts, ftype, denom, eps=EPS_DEFAULT,
        min_energy=1e-8):
    """Calculage differential of Bose-Einstein distribution, dn/dT
    Args
    -----
    ftype : string
        f... (Fermi-Dirac) or b... (Bose-Einstein)
    denom : string
        e... (df/dE) or t... (df/dT)
    energy : float
        Frequency [J]
    T : float, shppe=(ndat)
        temperature [K]
    """
    if type(energy) is int or type(energy) is float:
        get_diffrential_statistics2(energy, Ts, ftype, denom, 
                eps=EPS_DEFAULT, min_energy=1e-8)
        if abs(energy*unit.JToCm) < min_energy:
            return 0.
    if ftype.lower()[0] == "f":
        FB = 1.
    else:
        FB = -1.
    
    T_tmp = np.where(Ts<eps, eps, Ts)
    kbTs = unit.kb * T_tmp
    Xs = energy / kbTs
    if FB == -1:
        torr = eps
        Xs = np.where(Xs<torr, torr, Xs)
    
    dnde = - 1. / kbTs / (np.exp(Xs) + FB*2. + np.exp(-Xs))
     
    if denom.lower()[0] == "e":
        return dnde
    else:
        dndt = - unit.kb * Xs * dnde
        return dndt
    
def get_diffrential_statistics2(energies, T, ftype, denom, eps=EPS_DEFAULT,
        min_energy=1e-8):
    if abs(T) < 1e-5:
        return 0.
    if ftype.lower()[0] == "f":
        FB = 1.
    else:
        FB = -1.
    
    kbT = unit.kb * T
    Xs = energies / kbT
    if FB == -1:
        torr = eps
        Xs = np.where(Xs<torr, torr, Xs)
    
    dnde = - 1. / kbT / (np.exp(Xs) + FB*2. + np.exp(-Xs))
    
    if denom.lower()[0] == "e":
        return dnde
    else:
        dndt = - unit.kb * Xs * dnde
        return dndt
    
    

