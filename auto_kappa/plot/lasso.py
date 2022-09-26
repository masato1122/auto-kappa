# -*- coding: utf-8 -*-
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from .initialize import (set_matplot, set_axis, set_legend)
import glob

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
        print(" Output", figname)
    return fig

def _get_recommended_alpha(filename):
    lines = open(filename, 'r').readlines()
    for ll in lines:
        if "Minimum CVSCORE" in ll:
            data = float(ll.split()[-1])
            return data
    return None


