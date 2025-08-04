# 
# gruneisen.py
# 
# This script ...
# 
# Created on August 01, 2025
# Copyright (c) 2025 Masato Ohnishi
#
# This file is distributed under the terms of the MIT license.
# Please see the file 'LICENCE.txt' in the root directory
# or http://opensource.org/licenses/mit-license.php for information.
# 
import os
import sys

from auto_kappa.structure.crystal import get_automatic_kmesh
from auto_kappa.io.gruneisen import GruAll, Gruneisen
from auto_kappa.plot import make_figure

import logging
logger = logging.getLogger(__name__)

class GruneisenCalculator:

    def calculate_gruneisen_parameters(self, deltak=0.01, reciprocal_density=1500):
        """
        Calculate Grüneisen parameters.
        
        Args
        ----
        deltak : float
            Step size for k-point sampling of the phonon dispersion.
        
        reciprocal_density : int
            Number of k-points in the reciprocal space determining the k-mesh grids.
        
        """
        ## Start calculation
        line = "Calculate Grüneisen parameters:"
        msg  = f"\n {line}"
        msg += f"\n " + "-" * len(line)
        logger.info(msg)
        
        for mode in ['band', 'dos']:
            
            if mode == 'band':
                kwargs = {}
                
                fig_width = 2.6
                aspect = 0.7
                suffix = 'gruneisen'
                GruObj = Gruneisen
            
            elif mode == 'dos':
                kpts = get_automatic_kmesh(
                    self.primitive,
                    reciprocal_density=reciprocal_density, dim=self.dim)
                kwargs = {'kpts': kpts}
                
                fig_width = 2.3
                aspect = 0.9
                suffix = 'gru_all'
                GruObj = GruAll
            
            try:
                ### Make the working directory, write input file, and run the job
                self.write_alamode_input(propt=f'gruneisen_{mode}', **kwargs)
                self.run_alamode(propt=f'gruneisen_{mode}', ignore_log=False, verbose=False)
                
                ### Plot
                fig, axes = make_figure(ncols=1, nrows=1, fig_width=fig_width, aspect=aspect)
                ax = axes[0][0]
                
                ## file and figure names
                workdir = self.out_dirs['cube']['gruneisen']
                if workdir.startswith('/'):
                    workdir = os.path.relpath(workdir, os.getcwd())
                filename = f"{workdir}/{self.prefix}.{suffix}"
                figname = f"{workdir}/fig_{suffix}.png"
                
                obj = GruObj(filename)
                obj.plot(ax)
                fig.savefig(figname, dpi=600, bbox_inches='tight')
                
                msg = " Output : %s" % figname
                logger.info(msg)
                if mode == 'band':
                    logger.info("")
            
            except Exception as e:
                msg = f" Failed to calculate or plot Grüneisen parameters for {mode} mode."
                logger.error(msg)
                logger.exception(e)

