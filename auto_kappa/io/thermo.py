# 
# thermo.py
# 
# This script handles ".thermo" file
# 
# Created on August 19, 2025
# Copyright (c) 2025 Masato Ohnishi
#
# This file is distributed under the terms of the MIT license.
# Please see the file 'LICENCE.txt' in the root directory
# or http://opensource.org/licenses/mit-license.php for information.
# 
import os
import pandas as pd
import matplotlib.pyplot as plt

from auto_kappa.units import kb, JToEv, RyToJ, RyToEv
from auto_kappa.plot import set_axis, set_legend

import logging
logger = logging.getLogger(__name__)

class Thermo:
    """ Class for handling thermodynamic data from .thermo file
    
    How to use
    -----------
    >>> import matplotlib
    >>> import matplotlib.pyplot as plt
    >>> from auto_kappa.plot import set_matplot
    >>> from auto_kappa.io.thermo import Thermo
    >>> 
    >>> filename = "./mp-149/harm/bandos/Si.thermo"
    >>> 
    >>> fontsize = 7
    >>> fig_width = 2.5
    >>> aspect = 1.3
    >>> dpi = 600
    >>> figname = 'fig_thermo.png'
    >>> 
    >>> set_matplot(fontsize=fontsize)
    >>> fig = plt.figure(figsize=(fig_width, aspect*fig_width))
    >>> plt.subplots_adjust(hspace=0.05)
    >>> ax1 = plt.subplot(2, 1, 1)
    >>> ax2 = plt.subplot(2, 1, 2)
    >>> 
    >>> thermo = Thermo(filename)
    >>> thermo.plot(ax1, ax2)
    >>> 
    >>> fig.savefig(figname, dpi=dpi, bbox_inches='tight')
    """
    cmap = plt.get_cmap('tab10')
    
    def __init__(self, filename):
        self.filename = filename
        if os.path.exists(filename):
            self.set_thermo_data()
        else:
            self._data = None
        
    def set_thermo_data(self):
        # Read the .thermo file and return the data
        # "Cv/kB": specific heat at constant volume
        # "U[Ry]": internal energy
        # "S/kB": vibrational entropy
        # "U[Ry]": the Helmholtz free energy
        self._data = read_thermo_file(self.filename)
    
    def get_temperature(self):
        return self._data['T[K]'].values
    
    def get_specific_heat(self, unit='kB'):
        values = self._data['Cv/kB'].values
        return _convert_unit_from_kb(values, unit)

    def get_entropy(self, unit='kB'):
        values = self._data['S/kB'].values
        return _convert_unit_from_kb(values, unit)

    def get_internal_energy(self, unit='eV'):
        values = self._data['U[Ry]'].values
        return _convert_unit_from_Ry(values, unit)
    
    def get_free_energy(self, unit='eV'):
        values = self._data['F[Ry]'].values
        return _convert_unit_from_Ry(values, unit)

    def plot_specific_heat(self, ax, unit='meV',
                           color=cmap(0),
                           xlabel='Temperature (K)',
                           ylabel='Specific heat (meV)',
                           label='${\\rm C_v}$ (meV)'):
        xdat = self.get_temperature()
        ydat = self.get_specific_heat(unit=unit)
        _plot_data(ax, xdat, ydat, color=color, label=label)
        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)
    
    def plot_entropy(self, ax, unit='meV',
                     color=cmap(1),
                     xlabel='Temperature (K)',
                     ylabel='Vibrational entropy (meV)',
                     label='${\\rm S}$ (meV)'):
        xdat = self.get_temperature()
        ydat = self.get_entropy(unit=unit)
        _plot_data(ax, xdat, ydat, color=color, label=label)
        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)
    
    def plot_internal_energy(self, ax, unit='eV',
                             color=cmap(2),
                             xlabel='Temprature (K)',
                             ylabel='Internal energy (eV)',
                             label='U (eV)'):
        xdat = self.get_temperature()
        ydat = self.get_internal_energy(unit=unit)
        _plot_data(ax, xdat, ydat, color=color, label=label)
        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)
    
    def plot_free_energy(self, ax, unit='eV',
                         color=cmap(3),
                         xlabel='Temperature (K)',
                         ylabel='Free energy (eV)',
                         label='F (eV)'):
        xdat = self.get_temperature()
        ydat = self.get_free_energy(unit=unit)
        _plot_data(ax, xdat, ydat, color=color, label=label)
        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)
    
    def plot(self, ax1, ax2=None, unit1='kB', unit2='eV',
             title="Thermodynamic properties"):
        
        ax1.set_title(title)
        
        if unit1.lower() == 'kb':
            unit_label = "$k_B$"
        else:
            unit_label = unit1
        ylabel1 = "${\\rm C_v}$, S (%s)" % unit_label
        ax1.axhline(0, linestyle='-', color='grey', lw=0.3)
        self.plot_specific_heat(ax1, xlabel=None, ylabel=None, unit=unit1, label='Specific heat')
        self.plot_entropy(ax1, xlabel=None, ylabel=ylabel1, unit=unit1, label='Entropy')
        
        from matplotlib.ticker import MaxNLocator
        ax1.yaxis.set_major_locator(MaxNLocator(integer=True))
        
        ylabel2 = "Energy (%s)" % unit2
        ax2.axhline(0, linestyle='-', color='grey', lw=0.3)
        self.plot_internal_energy(ax2, xlabel=None, ylabel=None, unit=unit2, label='Internal energy')
        self.plot_free_energy(ax2, xlabel='Temperature (K)', ylabel=ylabel2, unit=unit2, label='Free energy')
        
        set_axis(ax1)
        set_axis(ax2)
        ax1.set_xticklabels([])
        
        set_legend(ax1, fontsize=6, alpha=0.5, loc='best')
        set_legend(ax2, fontsize=6, alpha=0.5, loc='best')


def _plot_data(ax, xdat, ydat, linestyle='-', lw=0.8, color='black', label=None):
    ax.plot(xdat, ydat, linestyle=linestyle, lw=lw, 
            color=color, marker=None, label=label)

def read_thermo_file(filename):    
    df = pd.read_csv(filename, sep='\s+', comment='#', header=None,
                     names=['T[K]', 'Cv/kB', 'S/kB', 'U[Ry]', 'F[Ry]'])
    return df

def _convert_unit_from_kb(values, unit):
    """ Convert the unit from kB to the specified unit """
    if unit.lower() == 'kb':
        return values
    elif unit == 'J':
        return values * kb
    elif unit == 'eV':
        return values * kb * JToEv
    elif unit == 'meV':
        return values * kb * JToEv * 1000
    else:
        logger.warning(f"\n Unknown unit: {unit}")
        return None

def _convert_unit_from_Ry(values, unit):
    """ Convert the unit from Ry to the specified unit """
    if unit.lower() == 'ry':
        return values
    elif unit == 'J':
        return values * RyToJ
    elif unit == 'eV':
        return values * RyToEv
    elif unit == 'meV':
        return values * RyToEv * 1000
    else:
        logger.warning(f"\n Unknown unit: {unit}")
        return None
