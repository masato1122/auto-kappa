# -*- coding: utf-8 -*-
#
# bandos.py
#
# Usuful functions to plot phonon band structure (.bands) and DOS (.dos)
# calculated with Alamode.
#
# Copyright (c) 2022 Masato Ohnishi
#
# This file is distributed under the terms of the MIT license.
# Please see the file 'LICENCE.txt' in the root directory
# or http://opensource.org/licenses/mit-license.php for information.
#
import os.path
import numpy as np
from optparse import OptionParser
from auto_kappa.plot.alamode.band import Band
from auto_kappa.plot.alamode.dos import Dos
from auto_kappa.plot.alamode.participation import Participation
from auto_kappa.plot.initialize import (
        set_matplot, set_axis, set_spaces, 
        get_both_axis, set_legend
        )

import matplotlib
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

def _get_band_files(directory, prefix):
    dir1 = directory + '/' + prefix
    extensions = ['bands', 'dos', 'pr', 'band.pr']
    filenames = {}
    for ext in extensions:
        filename = directory+'/'+prefix+'.'+ext
        if os.path.exists(filename):
            filenames [ext] = filename
    return filenames

def plot_bandos(directory='.', prefix=None, figname=None,
        directory2=None, prefix2=None,
        fig_width=4.0, fig_aspect=0.5,
        ymin=None, ymax=None, xmax2=None, yticks=None, myticks=None,
        fig_labels=[None, None],
        lw=0.4, lw2=0.6, wspace=0.05, 
        dpi=300, col='blue', col2='grey',
        unit='cm', legend_loc='best',
        plot_dos=True, plot_pdos=True, plot_pr=True,
        plot_dos2=False, 
        show_legend=True, show_colorbar=True
        ):
    """ Plot phonon dispersion and DOS. If .band.pr file (participation ratio) 
    is located in the given directory and ``plot_pr`` is True, it is shown with
    colors on phonon dispersion.
    
    Args
    -----
    directory : string, default=.
    prefix : string, default None
        {directory}/{prefix}.*** are searched.
    
    figname : string, default=None
        outputfigure name
    
    directory2 : string, default None
    prefix2 : string, default None
        Prefix for the second band dispersion
    
    fig_width : double, default 3.5
        figure width
    fig_aspect : double, default 0.5
        figure aspect

    ymin : double, default None
    ymax : double, default None
    xmax2 : double, default None
    
    yticks : double, default None
    myticks : int, default None

    lw : double, default 0.5
    lw2 : double
        line width for the second phonon dispersion
    wspace : double, default 0.05

    dpi : double, default 300
    col : default blue
    col2 : default grey
    unit : default cm
    legend_loc : default best

    plot_dos, plot_pdos, plot_pr : default True
    show_legend, show_colorbar : default True
    
    How to Use
    ------------
    >>> plot_bandos(directory='.', prefix='Si',
    >>>             figname='fig_bandos.png',
    >>>             plot_pr=True)
     
    """
    ### anphon file names
    filenames = _get_band_files(directory, prefix)
    
    ### check keys
    #for key in params:
    #    if key not in _get_default_figure_parameters():
    #        print(" Attention: %s is not defined." % key)
    
    ### get band
    from auto_kappa.plot.alamode.band import Band
    band = None
    if 'bands' not in filenames:
        print("")
        print(" Error: fn_band must be given.")
        print("")
        return None
    else:
        if os.path.exists(filenames['bands']) == False:
            print(" %s does not exist." % fn_band)
        else:
            band = Band(filename=filenames['bands'])
     
    ### get DOS
    from auto_kappa.plot.alamode.dos import Dos
    dos = None
    if 'dos' in filenames and plot_dos:
        if os.path.exists(filenames['dos']) == False:
            print(" %s does not exist." % fn_dos)
        else:
            dos = Dos(filename=filenames['dos'])
    
    ### min and max frequencies
    f0 = np.amin(band.frequencies)
    f1 = np.amax(band.frequencies)
    df = f1 - f0
    
    if ymin is None:
        fmin = f0 - df * 0.05
    else:
        fmin = ymin
    
    if ymax is None:
        if ('pr' in filenames or 'band.pr' in filenames) and plot_pr:
            fmax = f1 + df * 0.2
        else:
            fmax = f1 + df * 0.05
    else:
        fmax = ymax
    
    ### get participation ratio
    pr_ratio = None
    for ext in ['band.pr', 'pr']:
        if ext in filenames and plot_pr:
            prfile = filenames[ext]
            pr_ratio = Participation(file=prfile)
    
    ylabel = conv_unit(unit, band, dos)
    
    ##
    if fig_labels[0] is None:
        fig_labels[0] = prefix
    
    if fig_labels[1] is None:
        fig_labels[1] = prefix2
    
    global plt
    
    ### set figure
    set_matplot(fontsize=7)
    fig = plt.figure(figsize=(fig_width, fig_aspect*fig_width))
    
    plt.subplots_adjust(wspace=wspace)
    
    ax1, ax2 = get_both_axis(ratio='3:1')
    
    x2label = "DOS (a.u.)"
    ylabel = "Frequency (%s)" % ylabel
    ax1.set_ylabel(ylabel)
    ax2.set_xlabel(x2label)
    
    ### set ticks and labels
    def set_xticks_labels(ax, kmax, ksym, labels):
        
        dk_all = kmax
        
        label_mod = []
        for i, label in enumerate(labels):
            
            exception = False
            if i < len(labels)-1:
                dk_each = ksym[i+1] - ksym[i]
                fw_each = dk_each / dk_all
                if "|" in label and fw_each < 0.1:
                    names = label.split('|')
                    label_mod.append("${\\rm ^{%s}/_{%s}}$" % (
                        names[0], names[1]))
                    exception = True
            
            if exception == False:
                label_mod.append("${\\rm %s}$" % label.replace("G", "\\Gamma"))
        
        ###
        ax.set_xticks(ksym)
        ax.set_xticklabels(label_mod)
    
    set_xticks_labels(
            ax1, band.kpoints[-1],
            band.ksym, band.label)
    
    ax1 = set_axis(ax1, yticks=yticks, myticks=myticks)
    ax2 = set_axis(ax2, yticks=yticks, myticks=myticks)
    
    if pr_ratio is None:
        
        if prefix2 is not None:
            lab = fig_labels[0]
        else:
            lab = None

        _plot_bands(ax1, band.kpoints, band.frequencies, band.label, 
                col=col, lw=lw, zorder=2, label=lab)
    
    else:
        # --- coloring fllowing the participation ratio
        cdict = {
                'red':   ((0.0, 1.0, 1.0), (1.0, 0.0, 0.0)),
                'green': ((0.0, 0.0, 0.0), (1.0, 1.0, 1.0)),
                'blue':  ((0.0, 0.0, 0.0), (1.0, 0.0, 0.0))
                }
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
        isort = np.argsort(cdat)[::-1]
        sc = ax1.scatter(
                xdat[isort], 
                ydat[isort],
                c=cdat[isort],
                cmap=cmap,
                marker=".", s=0.3, lw=lw, zorder=2)
        ##
        if show_colorbar:
            set_colorbar(sc, ax=ax1)
    
    ## zero-line
    ax1.axhline(0, ls='-', lw=0.2, c='grey')
    
    ### plot DOS
    if dos is not None:
        ax2.plot(dos.dos, dos.frequencies, 
                marker="None", c=col, mew=lw, ms=1, mfc='none', lw=lw,
                zorder=2, label='total')
        
        xmin = 0.
        if xmax2 is None:
            xmax = np.max(dos.dos)
        else:
            xmax = xmax2
        dx = xmax - xmin
        ax2.set_xlim(xmin - 0.05*dx, xmax + 0.05*dx)
        
        ### plot PDOS
        if plot_pdos:
            
            pdos = _get_pdos(dos)
            
            if len(pdos) > 1:
                _plot_pdos(
                        ax2, dos.frequencies, pdos,
                        lw=lw*0.5, labels=dos.elements
                        )
                if show_legend:
                    set_legend(ax2, fs=6, alpha=0.5, loc=legend_loc)
    
    else:
        ax2.set_xlim(-0.1, 1.1)
        
    ### Second band and DOS
    if directory2 is not None and prefix2 is not None:
        
        filenames2 = _get_band_files(directory2, prefix2)
        
        if 'bands' in filenames2:
            band2 = Band(filename=filenames2['bands'])
            _plot_bands(ax1, band2.kpoints, band2.frequencies, band2.label, 
                    col=col2, lw=lw2, zorder=1, label=fig_labels[1])
        
        if 'dos' in filenames2 and plot_dos2:
            dos2 = Dos(filename=filenames2['dos'])
            ax2.plot(dos2.dos, dos2.frequencies, 
                    marker="None", c=col2,
                    mew=lw, ms=1, mfc='none', lw=lw,
                    zorder=1)
        
    ###        
    ax1.set_ylim(fmin, fmax)
    ax2.set_ylim(fmin, fmax)
    
    if prefix2 is not None:
        set_legend(ax1, fs=6, alpha=0.5, loc='best')
        
    if figname is not None:
        plt.savefig(figname, dpi=dpi, bbox_inches='tight')
        plt.close()
        print(" Output", figname)
    return fig

def _plot_bands(ax, ks_tmp, frequencies, xlabels, col='blue', lw=0.5, zorder=10,
        label=None):

    kpoints = ks_tmp.copy()
    
    nbands = len(frequencies[0])
    for ib in range(nbands):

        idx_zero = np.where(np.diff(kpoints) < 1e-10)[0]
        
        koffset = 0.
        for isec in range(len(idx_zero)+1):
            
            if isec == 0:
                i0 = 0
                i1 = idx_zero[0] + 1
            elif isec < len(idx_zero):
                i0 = idx_zero[isec-1] + 1
                i1 = idx_zero[isec] + 1
            else:
                i0 = idx_zero[isec-1] + 1
                i1 = len(kpoints)
            
            if ib == 0 and isec == 0:
                lab = label
            else:
                lab = None

            ax.plot(kpoints[i0:i1], frequencies[i0:i1,ib], 
                    marker="None", c=col,
                    mew=lw, ms=1, mfc='none', lw=lw,
                    zorder=zorder, label=lab)

def conv_unit(unit, band, dos):
    
    from auto_kappa.units import CmToHz, CmToEv

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


#def main(bfile, dfile, options):
#    band = Band(file=bfile)
#    if file_check(dfile):
#        dos = Dos(file=dfile)
#    else:
#        dos = None
#        options.plot_pdos = 0
#    
#    if options.plot_pdos != 0 and dos.elements is not None:
#        plot_pdos = True
#    else:
#        plot_pdos = False
#     
#    LABEL = get_label(bfile)
#    if options.figname is None:
#        figname = "fig_{:s}_bandos.png".format(LABEL)
#    else:
#        figname = options.figname
#    plot_bandos(
#            figname, band, dos, options, 
#            pdos=plot_pdos, wspace=options.wspace,
#            legend_loc=options.legend_loc
#            )
#
#if __name__ == '__main__':
#    parser = OptionParser()
#    
#    parser.add_option("--prefix", dest="prefix", type="string", 
#            help="prefix (prefix.band, prefix.dos will be read.)")
#    
#    parser.add_option("--figname", dest="figname", type="string", 
#            help="figure name")
#    
#    parser.add_option("--plot_pdos", dest="plot_pdos", type="int", 
#            default=1, help="0. not plot PDOS, 1. plot PDOS")
#    
#    parser.add_option("--prefix2", dest="prefix2", type="string", 
#            help="[optional] prefix for the second phonon band (prefix2.band "
#            "and prefix2.dos will be read when they are exist.)")
#    parser.add_option("--pr_ratio", dest="pr_ratio", type="int", 
#            default=0,
#            help="Flag for participation ratio (0.off, 1.on)")
#    
#    parser.add_option("--y0", dest="y0", type="float", 
#            help="Minimum frequency")
#    parser.add_option("--y1", dest="y1", type="float", 
#            help="Maximum frequency")
#    
#    parser.add_option("--maxdos", dest="maxdos", type="float", 
#            help="Maximum DOS")
#    
#    parser.add_option("--lw", dest="lw", type="float", 
#            default=0.3, help="line width of the figure (default: 0.3)")
#    parser.add_option("--yticks", dest="yticks", type="float", 
#            default=100, help="ticks of y-axis (default: 100)")
#    parser.add_option("--myticks", dest="myticks", type="int", 
#            default=2, help="mticks of y-axis (default: 2)")
#    
#    parser.add_option("--fig_width", dest="fig_width", type="float", 
#            default=3.3, help="width of figure (default: 3.3)") 
#    parser.add_option("--fig_aspect", dest="fig_aspect", type="float", 
#            default=1.0, help="aspect ratio of figure (default: 1.0)") 
#    
#    parser.add_option("--dpi", dest="dpi", type="int", 
#            default=300, 
#            help="dpi: resolution of the figure (default: 300)")
#    
#    parser.add_option("--col", dest="col", type="string", 
#            default="blue", help="line color (default: blue)")
#    parser.add_option("--col2", dest="col2", type="string", 
#            default="grey", help="line color of 2nd band (default: grey)")
#    
#    parser.add_option("--wspace", dest="wspace", type="float", 
#            default=0.05, help="wspace")
#    
#    parser.add_option("--unit", dest="unit", type="string", 
#            default="cm", 
#            help="unit of frequency, cm, THz, or meV. (default: cm(^-1))")
#    
#    parser.add_option("--colorbar", dest="colorbar", type="int", 
#            default=1, help="color bar")
#    
#    parser.add_option("--legend", dest="legend", type="int", 
#            default=1, help="legend")
#    
#    parser.add_option("--legend_loc", dest="legend_loc", type="string", 
#            default='best', help="legend location")
#    
#    (options, args) = parser.parse_args()
#    if options.prefix is None:
#        print("Input prefix")
#        sys.exit()
#    
#    prefix = options.prefix
#    BFILE = "{:s}.bands".format(prefix)
#    DFILE = "{:s}.dos".format(prefix)
#    if file_check(BFILE) is False:
#        print("Cannot find ", BFILE)
#        exit()
#
#    main(BFILE, DFILE, options)

