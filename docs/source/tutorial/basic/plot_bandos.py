#!/usr/bin/env python3
# Filename: plot_bandos.py
import os
import argparse
import glob

from auto_kappa.io import Band, Dos, BandPR
from auto_kappa.plot.initialize import set_matplot, set_axis

import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

def get_formula(base_dir):
    line = f"{base_dir}/harm/bandos/*.bands"
    fns = glob.glob(line)
    if len(fns) == 0:
        raise ValueError(" No band file found")
    fn = fns[0]
    return os.path.basename(fn).split(".")[0]

def main(options):
    
    fontsize = 7
    fig_width = 3.5
    aspect = 0.5
    dpi = 600
    
    ## Prepare figure frame
    set_matplot(fontsize=fontsize)
    fig = plt.figure(figsize=(fig_width, aspect*fig_width))
    plt.subplots_adjust(wspace=0.05)
    gs = gridspec.GridSpec(1, 2, width_ratios=[3,1])
    ax1 = fig.add_subplot(gs[0])
    ax2 = fig.add_subplot(gs[1])
    
    ## Get file names
    formula = get_formula(options.base_dir)
    # file_band = f"{options.base_dir}/harm/bandos/{formula}.bands"
    file_band_pr = f"{options.base_dir}/harm/bandos/{formula}.band.pr"
    file_dos = f"{options.base_dir}/harm/bandos/{formula}.dos"
    figname = "fig_bandos.png"
    
    ## Make Band and Dos objects
    # band = Band(file_band)
    band_pr = BandPR(file_band_pr)
    dos = Dos(file_dos)
    
    ## Plot band and DOS
    # band.plot(ax1, plot_G2G=False)
    band_pr.plot(ax1, plot_G2G=False, cbar_location='upper left')
    dos.plot(ax2, xlabel=None)
    
    ## Use same y limit
    ax2.set_ylim(ax1.get_ylim())
    
    ## Adjust figure
    ax2.tick_params(labelleft=False)
    set_axis(ax1)
    set_axis(ax2)
    
    fig.savefig(figname, dpi=dpi, bbox_inches='tight')
    print(" Output", figname)
    

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Input parameters')
    parser.add_argument('--base_dir', dest='base_dir', type=str,
                        default="./mp-22862", help="directory of auto-kappa output")
    args = parser.parse_args()    
    main(args)
    