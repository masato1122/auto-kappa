# -*- coding: utf-8 -*-
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from ..initialize import set_axis

class Band():
    def __init__(self, filename=None):
        """Band
        Variables
        -----------
        unit : string
            unit of frequency
        nk : integer
            # of k points
        nbands : integer
            # of bands
        kpoints : array, float, shape=(nk)
            k-points
        frequencies : ndarray, float, shape=(nk, nbands)
            frequencies
        label : array, string, shape=(nsym)
            Labels for symmetry points
        ksym : array, float, shape=(nsym)
            k-points for symmetry points
        nsym : integer
            # of symmetry points
        """
        self.unit = "cm^1"
        self.nk = None
        self.nbands = None
        # -- eigenvalues and vectors
        self.kpoints = None
        self.freqencies = None
        # -- symmetry points
        self.label = None
        self.ksym = None
        if filename is not None:
            self.read_bfile(filename)
    
    def set_nk_nbands(self, bfile):
        self.nk, self.nbands = get_nk_nbands(bfile)
    
    def set_symmetry_points(self, bfile):
        self.label, self.ksym = get_symmetry_points(bfile)
    
    def set_eigen(self, bfile):
        self.kpoints, self.frequencies = (
                get_eigen(bfile, self.nk, self.nbands))
     
    def read_bfile(self, bfile):
        """Read band file
        """
        self.set_nk_nbands(bfile)
        self.set_symmetry_points(bfile)
        self.set_eigen(bfile)
    
    #def plot_band(self, ax, lw=0.8, color='blue',
    #        ylabel="Frequency (cm$^{-1}$)",
    #        pr=None, fs_ticks=None
    #        ):
    #    """
    #    pr : Participation object
    #    """
    #    
    #    ## set ticks
    #    ax.set_xticks(self.ksym)
    #    if fs_ticks is not None:
    #        ax.set_xticklabels(self.label, fontsize=fs_ticks)
    #    else:
    #        ax.set_xticklabels(self.label)
    #    
    #    ## plot band
    #    if pr is None:
    #        for ib in range(self.nbands):
    #            ax.plot(self.kpoints, 
    #                    self.frequencies[:,ib], 
    #                    linestyle='-', lw=lw, 
    #                    c=color,
    #                    marker='None'
    #                    )
    #    else:
    #        cdict = {'red':   ((0.0, 1.0, 1.0), (1.0, 0.0, 0.0)),
    #                 'green': ((0.0, 0.0, 0.0), (1.0, 1.0, 1.0)),
    #                 'blue':  ((0.0, 0.0, 0.0), (1.0, 0.0, 0.0))}
    #        cmap = matplotlib.colors.LinearSegmentedColormap(
    #                'my_colormap', cdict, 256)
    #        ##
    #        kdump = np.zeros((self.nk, self.nbands))
    #        fdump = np.zeros((self.nk, self.nbands))
    #        cdump = np.zeros((self.nk, self.nbands))
    #        for ib in range(self.nbands):
    #            kdump[:,ib] = self.kpoints[:]
    #            fdump[:,ib] = self.frequencies[:,ib]
    #            cdump[:,ib] = pr.ratio[:,ib]
    #        kdump = kdump.reshape(self.nk*self.nbands)
    #        fdump = fdump.reshape(self.nk*self.nbands)
    #        cdump = cdump.reshape(self.nk*self.nbands)
    #        isort = np.argsort(cdump)[::-1]
    #        sc = ax.scatter(
    #                kdump[isort], fdump[isort], c=cdump[isort],
    #                cmap=cmap, marker='.', s=0.15, lw=lw, zorder=2)
    #        
    #        ## color bar
    #        from mpl_toolkits.axes_grid1.inset_locator import inset_axes
    #        axins = inset_axes(
    #                ax, width="50%", height="3%", loc="upper center",
    #                bbox_to_anchor=(0.0, -0.01, 1, 1),
    #                bbox_transform=ax.transAxes
    #                )
    #        [axins.spines[k].set_visible(False) for k in axins.spines]
    #        axins.set_facecolor([1,1,1,0.5])
    #        axins.tick_params(
    #                axis = 'both',
    #                left = False,
    #                top = False,
    #                right = False,
    #                bottom = False,
    #                labelleft = False,
    #                labeltop = False,
    #                labelright = False,
    #                labelbottom = False
    #                )
    #        sc.set_clim([0,1])
    #        cb = plt.colorbar(
    #                sc, cax=axins, orientation='horizontal', 
    #                ticks=[0,0.5,1])
    #        ## backgroud
    #        set_axis(axins)
    #
    #    ##
    #    ax.set_ylabel(ylabel)
    #    set_axis(ax)

def get_nk_nbands(bfile):
    """Get nk and nbands from band file
    Parameters
    ------------
    bfile : string
        band file name of ALAMODE
    
    Returns
    ----------
    nk, nband : integer
        # of k-points and bands
    """
    nline = sum(1 for line in open(bfile))
    nk = nline - 3
    ifs = open(bfile, "r")
    for i in range(4):
        line = ifs.readline()
    data = line.split()
    nbands = len(data) - 1
    return nk, nbands

def get_symmetry_points(bfile):
    """Read symmetry points from band file
    Returns
    ---------
    label : string
        label for symmetry points
    ksym : double
        |k| for symmetry points
    """
    ifs = open(bfile, "r")
    nline = sum(1 for line in open(bfile))
    line = ifs.readline(); data1 = line.split()
    line = ifs.readline(); data2 = line.split()
    label_tmp = []
    kpoints = []
    for i in range(len(data1)-1):
        label_tmp.append(data1[1+i])
        kpoints.append(float(data2[1+i]))
    label = []
    nticks = len(label_tmp)
    for it in range(nticks):
        lab0 = label_tmp[it]
        lab = label_tmp[it]
        if it != 0 and it != nticks-1:
            if abs(kpoints[it] - kpoints[it+1]) < 1e-5:
                lab = "%s|%s" % (lab0, label_tmp[it+1])
            if abs(kpoints[it] - kpoints[it-1]) < 1e-5:
                lab = "%s|%s" % (label_tmp[it-1], lab0)
        if "gamma" in lab.lower() or 'G' in lab:
            label.append("G")
        else:
            label.append("%s"%(lab))
    
    ### adjust
    lab_new = []
    knew = []
    lab_new.append(label[0])
    knew.append(kpoints[0])
    for ii in range(1,len(label)):
        if (abs(kpoints[ii] - kpoints[ii-1]) < 1e-5 and
                label[ii] == label[ii-1]):
            pass
        else:
            lab_new.append(label[ii])
            knew.append(kpoints[ii])
    return lab_new, knew

def get_eigen(bfile, nk, nbands):
    """Get kpoints and eigenvalues from band file
    Parameters
    ----------
    bfile : string
        band file name
    nk, nbnads : integer
        # of kpoints and bands
    Returns
    --------
    kpoints : float, shape=(nk)
        k-points
    freqs : float, shape=(nk, nbands)
        frequencies
    """
    ifs = open(bfile, "r")
    nline = sum(1 for line in open(bfile))
    kpoints = np.zeros(nk)
    freqs = np.zeros((nk, nbands))
    count = 0
    for il in range(nline):
        line = ifs.readline()
        data = line.split()
        if len(data) == 0:
            continue
        if line[0] == "#":
            continue
        kpoints[count] = float(data[0])
        for ib in range(nbands):
            freqs[count,ib] = float(data[1+ib])
        count += 1
    ifs.close()
    return kpoints, freqs

