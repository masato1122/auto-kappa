from auto_kappa.plot.bandos import plot_bandos

dir_band = "./bandos"
prefix = "Si"
figname = "fig_bandos.png"

plot_bandos(
        directory=dir_band, prefix=prefix,
        directory2=dir_band, prefix2=prefix,
        col2='red',
        figname=figname)

