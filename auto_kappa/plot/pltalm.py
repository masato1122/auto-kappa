# -*- coding: utf-8 -*-
#
# pltalm.py
#
# Helps to generate figures of phonon properties calculated with Alamode.
#
# Copyright (c) 2022 Masato Ohnishi
#
# This file is distributed under the terms of the MIT license.
# Please see the file 'LICENCE.txt' in the root directory
# or http://opensource.org/licenses/mit-license.php for information.
#
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from .initialize import (set_matplot, set_axis, set_legend)
import glob

def plot_kappa(df, figname='fig_kappa.png',
        dpi=300, fontsize=7, fig_width=2.3, aspect=0.9, lw=0.5, ms=2.0):

    cmap = plt.get_cmap("tab10")
    set_matplot(fontsize=fontsize)
    fig = plt.figure(figsize=(fig_width, aspect*fig_width))

    ax = plt.subplot()
    ax.set_xlabel('T (K)')
    ax.set_ylabel('${\\rm \\kappa (Wm^{-1}K^{-1})}$')

    markers = ['+', 'x', 'v', 'o']
    for ik, key in enumerate(['kxx', 'kyy', 'kzz', 'kave']):

        xdat = df['temperature'].values
        ydat = df[key].values
        ax.plot(xdat, ydat, linestyle='None', lw=lw,
                marker=markers[ik], markersize=ms,
                mfc='none', mew=lw, label=key
                )

    set_axis(ax, xscale='log', yscale='log')
    set_legend(ax, fs=6)

    fig.savefig(figname, dpi=dpi, bbox_inches='tight')
    print()
    print(" Output", figname)
    return fig

def _get_log_range(vmin, vmax, space=0.05):

    vmin_log = np.log10(vmin)
    vmax_log = np.log10(vmax)
    dv_log = vmax_log - vmin_log
    v0_log = vmin_log - dv_log * space
    v1_log = vmax_log + dv_log * space
    v0 = np.power(10, v0_log)
    v1 = np.power(10, v1_log)
    return [v0, v1]

def plot_lifetime(dfs, figname='fig_lifetime.png', xscale='linear',
        dpi=300, fontsize=7, fig_width=2.3, aspect=0.9, lw=0.3, ms=1.3):
    
    cmap = plt.get_cmap("tab10")
    set_matplot(fontsize=fontsize)
    fig = plt.figure(figsize=(fig_width, aspect*fig_width))

    ax = plt.subplot()
    ax.set_xlabel('Frequency (${\\rm cm^{-1}}$)')
    ax.set_ylabel('Lifetime (ps)')
    
    xlim = [100, -100]
    ylim = [100, -100]
    for it, key in enumerate(dfs):

        df = dfs[key]
        
        xorig = df['frequency'].values
        yorig = df['lifetime'].values
        
        ### get only available data
        idx_pos = np.where((xorig > 0.) & (yorig > 1e-5))[0]
        xdat = xorig[idx_pos]
        ydat = yorig[idx_pos]
        
        label = "%d K" % int(key)
        ax.plot(xdat, ydat, linestyle='None', lw=lw, 
            marker='o', markersize=ms, mew=lw, mec=cmap(it), mfc='None',
            label=label)
        
        ### axes range
        idx_sort = np.argsort(xdat)
        idx = np.where(xdat > 1e-5)[0]
        xmin = np.min(xdat[idx_sort[5:]])
        xmax = np.max(xdat)
        xlim[0] = min(xlim[0], xmin)
        xlim[1] = max(xlim[1], xmax)
        
        idx_sort = np.argsort(ydat)
        n = len(idx_sort)
        ymin = np.min(ydat[idx_sort[5:]])
        ymax = np.max(ydat)
        ylim[0] = min(ylim[0], ymin)
        ylim[1] = max(ylim[1], ymax)
    
    ax.set_ylim(_get_log_range(ylim[0], ylim[1], space=0.05))
    if xscale == 'log':
        ax.set_xlim(_get_log_range(xlim[0], xlim[1], space=0.05))
    else:
        xlim[0] = 0.
        xlim[0] = -0.05 * xlim[1]
        xlim[1] = xlim[1] * 1.05
        ax.set_xlim(xlim)

    set_axis(ax, xscale=xscale, yscale='log')
    set_legend(ax, fs=6, loc='best', alpha=0.5)
    
    fig.savefig(figname, dpi=dpi, bbox_inches='tight')
    print("")
    print(" Output", figname)

def plot_scattering_rates(frequencies, scat_rates, labels, 
        figname='fig_scat_rates.png', 
        dpi=300, fontsize=7, fig_width=2.3, aspect=0.9, lw=0.3, ms=1.3):
    
    cmap = plt.get_cmap("tab10")
    set_matplot(fontsize=fontsize)
    fig = plt.figure(figsize=(fig_width, aspect*fig_width))

    ax = plt.subplot()
    ax.set_xlabel('Frequency (${\\rm cm^{-1}}$)')
    ax.set_ylabel('Scattering rate (${\\rm ps^{-1}}$)')
    
    markers = ['o', '^', 's', 'd', 'v']
    xlim = [100, -100]
    ylim = [100, -100]
    for ik, key in enumerate(labels):

        xorig = frequencies.copy()
        yorig = scat_rates[key].copy()
        label = labels[key]
        
        ### get only available data
        idx_pos = np.where((xorig > 0.) & (yorig > 1e-5))[0]
        xdat = xorig[idx_pos]
        ydat = yorig[idx_pos]
        
        ax.plot(xdat, ydat, linestyle='None', lw=lw, 
            marker=markers[ik%len(markers)], 
            markersize=ms, mew=lw, mec=cmap(ik), mfc='None',
            label=label
            )
        
        ### axes range
        idx_sort = np.argsort(xdat)
        xmin = np.min(xdat[idx_sort[5:]])
        xmax = np.max(xdat[idx_sort[:-5]])
        
        idx_sort = np.argsort(ydat)
        n = len(idx_sort)
        ymin = np.min(ydat[idx_sort[5:]])
        ymax = np.max(ydat[idx_sort[:-5]])

        xlim[0] = min(xmin, xlim[0])
        xlim[1] = max(xmax, xlim[1])
        ylim[0] = min(ymin, ylim[0])
        ylim[1] = max(ymax, ylim[1])
        
    
    ax.set_xlim(_get_log_range(xlim[0], xlim[1], space=0.05))
    ax.set_ylim(_get_log_range(ylim[0], ylim[1], space=0.05))
    
    set_axis(ax, xscale='log', yscale='log')
    set_legend(ax, fs=6, loc='best', alpha=0.5)
    
    fig.savefig(figname, dpi=dpi, bbox_inches='tight')
    print("")
    print(" Output", figname)
    return fig

def plot_cumulative_kappa(dfs, 
        figname='fig_kcumu.png', xlabel=None, ylabel=None, xscale='linear',
        dpi=300, fontsize=7, fig_width=2.3, aspect=0.9, lw=0.8, ms=2.0):
    """ 
    dfs : dict of DataFrame
        key is temperature
    """
    cmap = plt.get_cmap("tab10")
    set_matplot(fontsize=fontsize)
    fig = plt.figure(figsize=(fig_width, aspect*fig_width))
        
    ax = plt.subplot()
    ax2 = ax.twinx()

    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    lab_kappa = "${\\rm \\kappa_{lat}}$"
    unit2 = "${\\rm Wm^{-1}K^{-1}}$"
    ax2.set_ylabel('Cumulative %s (%s)' % (lab_kappa, unit2))
    
    ylim1 = [0, -100]
    ylim2 = [0, -100]
    for it, key in enumerate(dfs):
        
        label = "%d K" % int(key)
        
        idx = np.where(dfs[key]['xdat'].values > 1e-3)[0]
        xdat = dfs[key]['xdat'].values[idx]
        ydat = dfs[key]['kspec'].values[idx]
        ydat2 = dfs[key]['kcumu'].values[idx]
        
        ax.plot(xdat, ydat, linestyle='-', c=cmap(it),
                lw=0.2, marker=None, zorder=10+it)
        ax.fill_between(xdat, 0, ydat, color=cmap(it), alpha=0.5, zorder=1)
        
        ax2.plot(xdat, ydat2, linestyle='-', lw=lw, marker=None, c=cmap(it),
                label=label)
        ##
        ylim1[1] = max(ylim1[1], np.max(ydat))
        ylim2[1] = max(ylim2[1], np.max(ydat2))
        
    dy1 = ylim1[1] - ylim1[0]
    dy2 = ylim2[1] - ylim2[0]
    ylim1[0] -= dy1 * 0.05
    ylim1[1] += dy1 * 0.05
    ylim2[0] -= dy2 * 0.05
    ylim2[1] += dy2 * 0.05
    
    ax.set_ylim(ylim1)
    ax2.set_ylim(ylim2)
    
    set_axis(ax, xscale=xscale)
    set_axis(ax2, xscale=xscale)
    
    ## legend
    h1, l1 = ax.get_legend_handles_labels()
    h2, l2 = ax2.get_legend_handles_labels()
    leg = ax2.legend(h1+h2, l1+l2,
            fontsize=6, labelspacing=0.3, handlelength=1.0,
            fancybox=False, edgecolor='black', loc='best')
    leg.get_frame().set_alpha(0.5)
    leg.get_frame().set_linewidth(0.2)
    
    ax.tick_params(labelright=False, right=False, which='both')
    ax2.tick_params(labelleft=False, left=False, which='both')
    
    fig.savefig(figname, dpi=dpi, bbox_inches='tight')
    print()
    print(" Output", figname)
    return fig

def plot_cvsets(directory='.', figname='fig_cvsets.png',
        dpi=300, fontsize=7, fig_width=2.0, aspect=0.9, lw=0.5, ms=0.5,
        ):
    
    cmap = plt.get_cmap("tab10")
    set_matplot(fontsize=fontsize)
    fig = plt.figure(figsize=(fig_width, aspect*fig_width))
        
    ax = plt.subplot()
    ax.set_xlabel('${\\rm \\alpha}$ (-)')
    ax.set_ylabel('Error (-)')
    
    count = 0
    flag_label = True
    ymax = 0.0
    for i in range(10):
        
        if i == 0:
            line = "%s/*.cvscore" % (directory)
        else:
            line = "%s/*.cvset%d" % (directory, i)
        fns = glob.glob(line)
        
        if len(fns) == 0:
            continue
        
        fn = fns[0]
        print(" Read", fn)
        data = np.genfromtxt(fn)
        
        xdat = data[:,0]
        if i == 0:
            ydat1 = data[:,1]
            ydat2 = data[:,3]
            col = 'black'
            zorder = 10
        else:
            ydat1 = data[:,1]
            ydat2 = data[:,2]
            col = cmap(count)
            zorder = i
        
        if flag_label:
            lab1 = "fitting"
            lab2 = "validation"
            flag_label = False
        else:
            lab1 = None
            lab2 = None
        
        ax.plot(xdat, ydat1, linestyle='--', lw=lw, c=col,
            marker='o', markersize=ms, 
            mew=lw, mec=col, mfc='None',
            label=lab1, zorder=zorder,
            )
        
        ax.plot(xdat, ydat2, linestyle='-', lw=lw, c=col,
            marker='^', markersize=ms, 
            mew=lw, mec=col, mfc='None',
            label=lab2, zorder=zorder
            )
        
        ### alpha
        if i == 0:
            alpha = _get_recommended_alpha(fn)
            ax.axvline(alpha, lw=lw, c='black')
            a = int(np.log10(alpha))
            b = alpha / np.power(10, float(a))
            text = "${\\rm %.2fx10^{%d}}$" % (b, a)
            ax.text(alpha, np.max(ydat1) * 0.9, text, 
                    fontsize=6, transform=ax.transData,
                    horizontalalignment="left", 
                    verticalalignment="center"
                    )
         
        ymax = max(ymax, np.max(ydat1))
        ymax = max(ymax, np.max(ydat2))
        if i > 0:
            count += 1
    
    ## y-lim
    ax.set_ylim(ymin=-0.05*ymax, ymax=ymax*1.03)
    
    set_axis(ax, xscale='log')
    set_legend(ax, fs=6, alpha=0.5)
    
    if figname is not None:
        fig.savefig(figname, dpi=dpi, bbox_inches='tight')
        print()
        print(" Output", figname)
    return fig

def _get_recommended_alpha(filename):
    lines = open(filename, 'r').readlines()
    for ll in lines:
        if "Minimum CVSCORE" in ll:
            data = float(ll.split()[-1])
            return data
    return None


