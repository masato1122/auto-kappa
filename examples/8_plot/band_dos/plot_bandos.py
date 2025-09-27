
import argparse
from auto_kappa.io.band import Band
from auto_kappa.io.dos import Dos
from auto_kappa.plot.initialize import set_axis

import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

def main(options):
    
    fig_width = 4.0
    aspect = 0.5
    
    ## Prepare figure
    fig = plt.figure(figsize=(fig_width, aspect*fig_width))
    plt.subplots_adjust(wspace=0.05)
    gs = gridspec.GridSpec(1, 2, width_ratios=[3,1])
    ax1 = fig.add_subplot(gs[0])
    ax2 = fig.add_subplot(gs[1])
    
    ## File names
    file_band = f"{options.dir_bandos}/{options.formula}.bands"
    file_dos = f"{options.dir_bandos}/{options.formula}.dos"
    figname = "fig_bandos.png"
    
    ## Make Band and Dos objects
    band = Band(file_band)
    dos = Dos(file_dos)
    
    ## Plot data
    band.plot(ax1)
    dos.plot(ax2, xlabel=None)
    
    ## Set limits
    ax2.set_ylim(ax1.get_ylim())
    
    ## Adjust figure
    ax2.tick_params(labelleft=False)
    set_axis(ax1)
    set_axis(ax2)
    
    fig.savefig(figname, dpi=600, bbox_inches='tight')
    print(" Output", figname)


if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Input parameters')

    parser.add_argument('--dir_bandos', dest='dir_bandos', type=str,
            default="./bandos", help="directory name for bandos")
    
    parser.add_argument('--formula', dest='formula', type=str,
            default="Si", help="formula")
    
    args = parser.parse_args()
    
    main(args)

