from auto_kappa.plot.bandos import plot_bandos
from auto_kappa.plot.initialize import set_axis, set_legend

dir_band = "./bandos"
prefix = "Si"
figname = "fig_bandos.png"

fig = plot_bandos(
        directory=dir_band, prefix=prefix,
        directory2=dir_band, prefix2=prefix,
        col2='red',
        figname=figname)

