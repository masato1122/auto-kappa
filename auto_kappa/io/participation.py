# 
# participation.py
# 
# This script handles participation ratio.
# 
# Created on August 02, 2025
# Copyright (c) 2025 Masato Ohnishi
#
# This file is distributed under the terms of the MIT license.
# Please see the file 'LICENCE.txt' in the root directory
# or http://opensource.org/licenses/mit-license.php for information.
# 
import os
import numpy as np
import matplotlib.pyplot as plt

from auto_kappa.plot import get_customized_cmap
from auto_kappa.io.band import Band

import logging
logger = logging.getLogger(__name__)

class BandPR:
    """ Class to handle participation ratio (.band.pr file).
    If the band file exists in the same directory, it will also read 
    the band structure from ".bands" file.
    
    How to use
    ----------
    >>> from auto_kappa.io.participation import BandPR
    >>> filename = "./mp-149/harm/bandos/Si.band.pr"
    >>> band_pr = BandPR(filename)
    >>> band_pr.plot(ax)  # Plot participation ratio with band structure
    """
    def __init__(self, filename):
        
        if os.path.isabs(filename):
            filename = os.path.relpath(filename, os.getcwd())
        self._file_pr = filename
        self._file_bands = filename.replace('.band.pr', '.bands')
        self._kpoints = None
        self._prs = None
        self.band = None
        try:
            self._read_file_pr()
        except Exception as e:
            msg = f"\n Error reading data from {self._file_pr}: {e}"
            logger.warning(msg)
        try:
            self._read_file_bands()
        except Exception as e:
            msg = f"\n Error reading data from {self._file_bands}: {e}"
            logger.warning(msg)
    
    @property
    def file_pr(self):
        return self._file_pr
    @property
    def file_bands(self):
        return self._file_bands
    @property
    def kpoints(self):
        if self._kpoints is None:
            self._read_file_pr()
        return self._kpoints
    @property
    def prs(self):
        if self._prs is None:
            self._read_file_pr()
        return self._prs
    @property
    def nbands(self):
        return self.prs.shape[1]

    def _read_file_pr(self):
        self._kpoints, self._prs = parse_file_pr(self.file_pr)
        
    def _read_file_bands(self):
        self.band = Band(self.file_bands)
    
    def plot(self, ax, lw=None, cmap='rainbow', cbar_location='right', plot_G2G=False):
        
        assert len(self.kpoints) == len(self.band.kpoints), \
            "Number of kpoints in participation ratio and band structure must match."
        
        cmap = get_customized_cmap(100, color1='red', color2='blue')
        norm = plt.Normalize(0, 1)
        
        self.band.plot_with_weighted_colors(
            ax, self.prs, lw=lw, cmap=cmap, norm=norm, 
            cbar_location=cbar_location, plot_G2G=plot_G2G,
            clabel='Participation Ratio')


def parse_file_pr(filepath):
    """ Read participation ratio from .band.pr file
    """
    kpoints = []
    prs = []
    
    with open(filepath, 'r') as f:
        lines = f.readlines()
    
    for line in lines:
        line = line.strip()
        if not line:
            continue
        
        ## kpoint line
        if line.startswith("#") and "xk =" in line:
            kpoint = list(map(float, line.split("=")[1].strip().split()))
            kpoints.append(kpoint)
            prs.append([])
            continue
        
        ## PR data
        parts = line.split()
        if len(parts) == 3:
            pr_val = float(parts[2])
            prs[-1].append(pr_val)
    
    return np.array(kpoints), np.array(prs)
