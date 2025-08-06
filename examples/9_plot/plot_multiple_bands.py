
from auto_kappa.plot import make_figure, set_matplot, set_legend
from auto_kappa.plot.initialize import prepare_two_axes
from auto_kappa.io.band import Band
from auto_kappa.io.dos import Dos

import matplotlib.pyplot as plt

def main():
    
    fontsize = 7
    fig_width = 3.0
    aspect = 0.7
    lw = 0.5
    
    ## Prepare figure
    set_matplot(fontsize=fontsize)
    fig = plt.figure(figsize=(fig_width, aspect*fig_width))
    ax1, ax2 = prepare_two_axes(ratio="2:1")
    cmap = plt.get_cmap('tab10')
    
    ## Prepare filenames and labels
    filenames = {}
    filenames["2x2x2"] = "./mp-21276-2x2x2/harm/bandos/PbS.bands"
    filenames["3x3x3"] = "./mp-21276-3x3x3/harm/bandos/PbS.bands"
    
    ## Plot bands and DOS
    bands = {}
    doses = {}
    for i, key in enumerate(filenames):
        fn = filenames[key]
        
        bands[key] = Band(fn)                          # Load .bands file
        doses[key] = Dos(fn.replace(".bands", ".dos")) # Load .dos file
        
        bands[key].plot(ax1, color=cmap(i), lw=lw, plot_G2G=True, label=key)
        doses[key].plot(ax2, color=cmap(i), lw=lw, xlabel=None, 
                        plot_pdos=False, show_legend=False)
    
    set_legend(ax1, fontsize=6, loc='lower left', loc2=(0.0, 1.0), ncol=1)
    
    figname = 'fig_bands.png'
    fig.savefig(figname, dpi=600, bbox_inches='tight')
    print(" Output", figname)

if __name__ == '__main__':
    main()
