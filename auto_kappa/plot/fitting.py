# -*- coding: utf-8 -*-
#
# fitting.py
#
# Helps to generate figure of the energy-volume relationship.
#
# Copyright (c) 2022 Masato Ohnishi
#
# This file is distributed under the terms of the MIT license.
# Please see the file 'LICENCE.txt' in the root directory
# or http://opensource.org/licenses/mit-license.php for information.
#
import numpy as np
import os
import sys

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from .initialize import (set_matplot, set_axis, set_legend)
import glob

import logging
logger = logging.getLogger(__name__)

def plot_fitting_result(bm, figname='fig.png', color='black',
        dpi=300, fontsize=7, fig_width=2.3, aspect=0.6, lw=0.5, ms=2.3):
    
    set_matplot(fontsize=fontsize)
    fig = plt.figure(figsize=(fig_width, aspect*fig_width))
    
    ### check fitting mode
    from pymatgen.analysis.eos import BirchMurnaghan
    if type(bm) is not BirchMurnaghan:
        msg = "\n Warning: fitting type of %s is not supported." % type(bm)
        logger.warning(msg)
        return None
    
    ax = plt.subplot()
    ax.set_xlabel('Volume (${\\rm \\AA}$)')
    ax.set_ylabel('Energy (eV)')
    
    ### plot original data
    xdat = bm.volumes
    ydat = bm.energies
    ax.plot(xdat, ydat, linestyle='None', lw=lw, 
            marker='o', markersize=ms, mec=color, mfc='none', mew=lw)
    
    ### plot fitting result
    nfit = 51
    dx = np.max(xdat) - np.min(xdat)
    xmin = np.min(xdat) - 0.01 * dx
    xmax = np.max(xdat) + 0.01 * dx
    xfit = np.linspace(xmin, xmax, nfit)
    yfit = bm.func(xfit)
    
    ax.plot(xfit, yfit, linestyle='--', c=color, lw=lw, marker=None)
    
    ### range of y-axis
    dy = np.max(yfit) - np.min(yfit)
    ymax = np.max(yfit) + 0.4 * dy
    ax.set_ylim(ymax=ymax)
    
    ### show text
    text = "Minimum energy = %.2f eV" % (bm.e0)
    text += "\nMinimum or reference volume = %.2f eV" % (bm.v0)
    text += "\nBulk modulus = %.2f ${\\rm eV/\\AA^3}$ GPa" % (bm.b0_GPa.real)
    text += "\nDerivative of bulk modulus wrt pressure = %.2f" % (bm.b1)
    ax.text(0.03, 0.97, text, fontsize=4, transform=ax.transAxes,
            horizontalalignment="left", verticalalignment="top",
            bbox=dict(facecolor='white', edgecolor='none', alpha=0.7, pad=0.0))
    
    set_axis(ax)
    fig.savefig(figname, dpi=dpi, bbox_inches='tight')
    msg = " Output %s" % figname
    logger.info(msg)
    return fig

