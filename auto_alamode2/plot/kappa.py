# -*- coding: utf-8 -*-
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from .initialize import (set_matplot, set_axis, set_legend)

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
    
    set_axis(ax, xformat='log', yformat='log')
    set_legend(ax, fs=6)
    
    fig.savefig(figname, dpi=dpi, bbox_inches='tight')
    print(" Output", figname)
    return fig

