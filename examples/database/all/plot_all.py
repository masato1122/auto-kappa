# -*- coding: utf-8 -*-
import numpy as np
from optparse import OptionParser
import glob
from pymatgen.core.composition import Composition

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from auto_alamode2.plot.initialize import (set_matplot, set_axis, set_legend)

def get_all_data():
    
    all_data = []
    line = '../mp-*/cube/kappa/*.kl'
    fns = glob.glob(line)
    for fn in fns:
        if 'boundary' in fn:
            continue
        words = fn.split('/')
        mpid = words[1]
        comp = words[-1].split('.')[0]
        print(mpid, comp)
        all_data.append({
            'mpid': mpid,
            'composition': comp,
            'data': np.genfromtxt(fn)
            })
    return all_data

def plot_all_data(all_data, figname='fig_all_kappas.png', 
        dpi=300, fontsize=7, fig_width=2.3, aspect=0.9, lw=0.5, ms=0.5):
    
    cmap = plt.get_cmap("tab10")
    set_matplot(fontsize=fontsize)
    fig = plt.figure(figsize=(fig_width, aspect*fig_width))
        
    ax = plt.subplot()
    ax.set_xlabel('T (K)')
    ax.set_ylabel('${\\rm \\kappa_{lat} (Wm^{-1}K^{-1})}$')
    
    for ii, each in enumerate(all_data):
        
        comp = Composition(each['composition']).as_dict()
        lab_comp = ""
        for cc in comp:
            if int(comp[cc]) == 1:
                lab_comp += cc
            else:
                lab_comp += "%s_%d" % (cc, int(comp[cc]))
        label = "${\\rm %s}$ (%s)" % (lab_comp, each['mpid'])
        
        xdat = each['data'][:,0]
        ydat = np.zeros(len(xdat))
        for j in range(3):
            ydat += each['data'][:,1+4*j] / 3.
        ##
        icol = ii % 10
        ax.plot(xdat, ydat, linestyle='-', c=cmap(icol),
                lw=lw, marker='o', markersize=ms,
                mfc='none', mew=lw, label=label)
    
    ##
    set_axis(ax, xformat='log', yformat='log')
    set_legend(ax, fs=5, loc='upper left', loc2=[1.0,1.0])
    
    fig.savefig(figname, dpi=dpi, bbox_inches='tight')
    print(" Output", figname)
    return fig


def main(options):

    all_data = get_all_data()
    
    plot_all_data(all_data)


if __name__ == '__main__':
    parser = OptionParser()
    parser.add_option("-f", "--filename", dest="filename", type="string",
        help="input file name")
    (options, args) = parser.parse_args()
    #file_check(options.filename)
    main(options)
