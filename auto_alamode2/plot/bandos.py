# -*- coding: utf-8 -*-
import os.path
import numpy as np
from optparse import OptionParser
from .alamode.band import Band
from .alamode.dos import Dos
from .alamode.participation import Participation
from .initialize import (
        set_matplot, set_axis, set_spaces, 
        get_both_axis, set_legend
        )

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import cm

def file_check(file):
    if os.path.exists(file) == False:
        return False
    return True

def get_label(FILE):
    line = FILE.split("/")
    line = line[ len(line)-1 ]
    line = line.split(".")
    label = line[ len(line)-2 ]
    return label

def _get_pdos(dos):
    """ 
    dos : Dos object in mytool
    """
    nene = len(dos.dos_atom)
    nel = len(dos.nat_el)
    pdos = np.zeros((nel, nene))
    i1 = 0
    for iel in range(nel):
        i0 = i1
        i1 = i0 + dos.nat_el[iel]

        ratio = 1.
        #ratio = float(dos.nat_el[iel])

        pdos[iel, :] = (
                np.sum(dos.dos_atom[:,i0:i1], axis=1) / ratio)
    return pdos

def _plot_pdos(ax, frequencies, pdos, lw=0.8, labels=None):
    
    cmap = plt.get_cmap("tab10")
    nelements = len(pdos)
    for ie in range(nelements):
        ax.plot(pdos[ie,:], frequencies,
                linestyle='-', c=cmap(ie),
                lw=lw, label=labels[ie])

def set_colorbar(sc, ax):
    
    from mpl_toolkits.axes_grid1.inset_locator import inset_axes
    axins = inset_axes(ax,
            width="50%", height="3%", 
            loc='upper center',
            bbox_to_anchor=[0, -0.03, 1, 1],
            bbox_transform=ax.transAxes,
            borderpad=0
            )
    set_axis(axins)
    cb = plt.colorbar(sc, cax=axins, orientation='horizontal')
    sc.set_clim([0.0, 1.0])

def plot_bandos(figname, band, dos, options, pdos=False, wspace=0.05,
        legend_loc='upper right'):
    """
    figname : string
        output figure name
    band : Band class
        .kpoints : array, float, shape=(nk)
        .frequencies : ndarray, float, shape=(nk,nbands)
    dos : Dos class
        .frequencies : array, float, shape=(nfreq)
        .dos : array, float, shape=(nfreq)
    """
    # --- get participation ratio
    prfile = None
    pr_ratio = 0
    if options.pr_ratio == 1:
        prfile1 = "{:s}.pr".format(options.prefix)
        prfile2 = "{:s}.band.pr".format(options.prefix)
        if file_check(prfile1):
            prfile = prfile1
        elif file_check(prfile2):
            prfile = prfile2
        
        if prfile is not None:
            pr_ratio = Participation(file=prfile)
        
    ylabel = conv_unit(options.unit, band, dos)

    global plt
    set_matplot(fontsize=7)
    fig = plt.figure(
            figsize=(options.fig_width, options.fig_aspect*options.fig_width))
    plt.subplots_adjust(left=0.14, bottom=0.15, right=0.98, top=0.98,
            wspace=wspace)
    
    lw = options.lw
    # ---- ax1 : band, ax2 : dos
    x2label = "DOS (a.u.)"
    ylabel = "Frequency ({:s})".format(ylabel)
    fmin, fmax = np.min(band.frequencies), np.max(band.frequencies)
    frange = [fmin - (fmax-fmin) * 0.05, fmax + (fmax-fmin) * 0.05]
    if options.y0 is not None:
        frange[0] = options.y0
    if options.y1 is not None:
        frange[1] = options.y1
    
    ax1, ax2 = get_both_axis(frange, ylabel, band.ksym, band.label, x2label)
    
    def set_yticks(options, ax):
        yticks = None
        if options.yticks is not None:
            yticks = options.yticks
        myticks = None
        if options.myticks is not None:
            myticks = options.myticks
        return set_axis(ax, yticks=yticks, myticks=myticks)
    ax1 = set_yticks(options, ax1)
    ax2 = set_yticks(options, ax2)
    
    if options.pr_ratio == 0:
        for ib in range(band.nbands):
            ax1.plot(band.kpoints, band.frequencies[:,ib], 
                    marker="None", c=options.col, 
                    mew=lw, ms=1, mfc='none', lw=lw,
                    zorder=2)
    else:
        # --- coloring fllowing the participation ratio
        cdict = {'red':   ((0.0, 1.0, 1.0), (1.0, 0.0, 0.0)),
                 'green': ((0.0, 0.0, 0.0), (1.0, 1.0, 1.0)),
                 'blue':  ((0.0, 0.0, 0.0), (1.0, 0.0, 0.0))}
        cmap = matplotlib.colors.LinearSegmentedColormap(
                'my_colormap',cdict,256)
        ## get all data
        nk = band.nk
        nb = band.nbands
        xdat = np.zeros((nk,nb))
        for ib in range(band.nbands):
            xdat[:,ib] = band.kpoints[:]
        xdat = xdat.reshape(nk*nb)
        ydat = band.frequencies.reshape(nk*nb)
        cdat = pr_ratio.ratio.reshape(nk*nb)
        #cdat[np.argmin(cdat)] = 0.0
        isort = np.argsort(cdat)[::-1]
        sc = ax1.scatter(
                xdat[isort], 
                ydat[isort],
                c=cdat[isort],
                cmap=cmap,
                marker=".", s=0.3, lw=lw, zorder=2)
        ##
        if options.colorbar == 1:
            set_colorbar(sc, ax=ax1)
     
    if dos is not None:
        ax2.plot(dos.dos, dos.frequencies, 
                marker="None", c=options.col, mew=lw, ms=1, mfc='none', lw=lw,
                zorder=2, label='total')
        
        ##xmin = np.min(dos.dos)
        xmin = 0.
        if options.maxdos is None:
            xmax = np.max(dos.dos)
        else:
            xmax = options.maxdos
        xlength = xmax - xmin
        ax2.set_xlim(xmin - 0.1*xlength, xmax + 0.1*xlength)

        ##
        if pdos:
            pdos = _get_pdos(dos)
            _plot_pdos(
                    ax2, dos.frequencies, pdos,
                    lw=lw*0.5, labels=dos.elements)

            if options.legend == 1:
                set_legend(ax2, fs=6, alpha=0.5, loc=legend_loc)
    else:
        ax2.set_xlim(-0.1, 1.1)
        
    # --- read the second band
    if options.prefix2 is not None:
        bfile2 = "{:s}.bands".format(options.prefix2)
        dfile2 = "{:s}.dos".format(options.prefix2)
        file_check(bfile2)
        band2 = Band(file=bfile2)
        for ib in range(band2.nbands):
            ax1.plot(band2.kpoints, band2.frequencies[:,ib], 
                    marker="None", c=options.col2, mew=lw,
                    ms=1, mfc='none', lw=lw*0.8, zorder=1)
        if file_check(dfile2):
            dos2 = Dos(file=dfile2)
            ax2.plot(dos2.dos, dos2.frequencies, 
                    marker="None", c=options.col2,
                    mew=lw, ms=1, mfc='none', lw=lw,
                    zorder=1)
        
        # --- reset axis range
        if np.min(band2.frequencies) < frange[0] and options.y0 is None:
            frange[0] = np.min(band2.frequencies)
        if np.max(band2.frequencies) > frange[1] and options.y1 is None:
            frange[1] = np.max(band2.frequencies)
        
        ax1.set_ylim(frange[0], frange[1])
        ax2.set_ylim(frange[0], frange[1])

        ##
        if file_check(dfile2) and dos is not None:
            idx_plot = np.where(
                    frange[0] < dos2.frequencies[idx_plot] < frange[1])[0]
            xmax = np.max(dos2.dos[idx_plot])
            ax2.set_ylim(xmax=xmax*1.2)
            print(idx_plot)
    
    plt.savefig(figname, dpi=options.dpi, bbox_inches='tight')
    plt.close()
    print(" Output", figname)
    return fig

def conv_unit(unit, band, dos):
    
    from ..units import CmToHz, CmToEv

    if unit.lower() == "thz":
        unit_conv = CmToTHz
        ylabel = "THz"
    elif unit.lower() == "mev":
        unit_conv = CmToEv * 1e3
        ylabel = "meV"
    elif unit.lower() == "cm" or unit.lower() == "kayser":
        unit_conv = 1.
        ylabel = "$\\rm{cm^{-1}}$"
    else:
        print("Error: {:s} is not defined.".format(unit))
        sys.exit()
    band.frequencies *= unit_conv
    if dos is not None:
        dos.frequencies *= unit_conv
    return ylabel

def main(bfile, dfile, options):
    band = Band(file=bfile)
    if file_check(dfile):
        dos = Dos(file=dfile)
    else:
        dos = None
        options.plot_pdos = 0
    
    if options.plot_pdos != 0 and dos.elements is not None:
        plot_pdos = True
    else:
        plot_pdos = False
     
    LABEL = get_label(bfile)
    if options.figname is None:
        figname = "fig_{:s}_bandos.png".format(LABEL)
    else:
        figname = options.figname
    plot_bandos(
            figname, band, dos, options, 
            pdos=plot_pdos, wspace=options.wspace,
            legend_loc=options.legend_loc
            )

if __name__ == '__main__':
    parser = OptionParser()
    
    parser.add_option("--prefix", dest="prefix", type="string", 
            help="prefix (prefix.band, prefix.dos will be read.)")
    
    parser.add_option("--figname", dest="figname", type="string", 
            help="figure name")
    
    parser.add_option("--plot_pdos", dest="plot_pdos", type="int", 
            default=1, help="0. not plot PDOS, 1. plot PDOS")
    
    parser.add_option("--prefix2", dest="prefix2", type="string", 
            help="[optional] prefix for the second phonon band (prefix2.band "
            "and prefix2.dos will be read when they are exist.)")
    parser.add_option("--pr_ratio", dest="pr_ratio", type="int", 
            default=0,
            help="Flag for participation ratio (0.off, 1.on)")
    
    parser.add_option("--y0", dest="y0", type="float", 
            help="Minimum frequency")
    parser.add_option("--y1", dest="y1", type="float", 
            help="Maximum frequency")
    
    parser.add_option("--maxdos", dest="maxdos", type="float", 
            help="Maximum DOS")
    
    parser.add_option("--lw", dest="lw", type="float", 
            default=0.3, help="line width of the figure (default: 0.3)")
    parser.add_option("--yticks", dest="yticks", type="float", 
            default=100, help="ticks of y-axis (default: 100)")
    parser.add_option("--myticks", dest="myticks", type="int", 
            default=2, help="mticks of y-axis (default: 2)")
    
    parser.add_option("--fig_width", dest="fig_width", type="float", 
            default=3.3, help="width of figure (default: 3.3)") 
    parser.add_option("--fig_aspect", dest="fig_aspect", type="float", 
            default=1.0, help="aspect ratio of figure (default: 1.0)") 
    
    parser.add_option("--dpi", dest="dpi", type="int", 
            default=300, 
            help="dpi: resolution of the figure (default: 300)")
    
    parser.add_option("--col", dest="col", type="string", 
            default="blue", help="line color (default: blue)")
    parser.add_option("--col2", dest="col2", type="string", 
            default="grey", help="line color of 2nd band (default: grey)")
    
    parser.add_option("--wspace", dest="wspace", type="float", 
            default=0.05, help="wspace")
    
    parser.add_option("--unit", dest="unit", type="string", 
            default="cm", 
            help="unit of frequency, cm, THz, or meV. (default: cm(^-1))")
    
    parser.add_option("--colorbar", dest="colorbar", type="int", 
            default=1, help="color bar")
    
    parser.add_option("--legend", dest="legend", type="int", 
            default=1, help="legend")
    
    parser.add_option("--legend_loc", dest="legend_loc", type="string", 
            default='best', help="legend location")
    
    (options, args) = parser.parse_args()
    if options.prefix is None:
        print("Input prefix")
        sys.exit()
    
    prefix = options.prefix
    BFILE = "{:s}.bands".format(prefix)
    DFILE = "{:s}.dos".format(prefix)
    if file_check(BFILE) is False:
        print("Cannot find ", BFILE)
        exit()

    main(BFILE, DFILE, options)

