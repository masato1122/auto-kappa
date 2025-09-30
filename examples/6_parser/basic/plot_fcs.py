#!/usr/bin/env python3
# Filename: plot_fcs.py
import os
import argparse
import glob

from auto_kappa.io import FCSxml
from auto_kappa.plot.initialize import set_matplot, set_axis

import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

def main(options):
    
    fontsize = 7
    fig_width = 2.8
    aspect = 1.0
    dpi = 600
    
    ## Prepare figure frame
    set_matplot(fontsize=fontsize)
    fig = plt.figure(figsize=(fig_width, aspect*fig_width))
    plt.subplots_adjust(hspace=0.05)
    axes = []
    for ifig in range(2):
        axes.append(plt.subplot(2, 1, ifig+1))
    
    ## Get file name
    from plot_bandos import get_formula
    formula = get_formula(options.base_dir)
    line = f"{options.base_dir}/cube/*/{formula}.xml"
    try:
        file_fc3 = glob.glob(line)[0]
    except IndexError:
        print(" No FCS file found")
        return
    
    ## Read FCS file
    fcs = FCSxml(file_fc3)
    
    ## Plot FCs
    fcs.plot_fc2(axes[0])
    fcs.plot_fc3(axes[1])
    axes[0].tick_params(labelbottom=False)
    
    ## Save figure
    figname = "fig_fcs.png"
    fig.savefig(figname, dpi=dpi, bbox_inches='tight')
    print(" Output", figname)
    
    print("FC2:", fcs.fc2.shape)
    print("FC3:", fcs.fc3.shape)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Input parameters')
    parser.add_argument('--base_dir', dest='base_dir', type=str,
                        default="./mp-22862", help="directory of auto-kappa output")
    args = parser.parse_args()    
    main(args)
    