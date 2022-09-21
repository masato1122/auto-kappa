==============
Useful Tips
==============


Plot phonon dispersion and DOS
=================================

``plot_bandos`` helps to plot the phonon dispersion and DOS.
Different files such as .bands, .dos, .band.pr are supposed to be in the same directory.

.. code-block:: python

    from auto_alamode2.plot.bandos import plot_bandos

    plot_bandos(directory='.', prefix='Si',
                figname='fig_bandos.png')
..



Make DFSET file for Alamode
=============================

``get_dfset`` helps to make a DFSET file containing displacements and forces extracted from many vasprun.xml files.
Many directories which contain a vasprun.xml file are supposed to be located under the given directory.

.. code-block:: python

    from auto_alamode2.io.vasp import get_dfset
    
    directory = './mp-149/harm/force'
    offset_xml = directory + '/prist/vasprun.xml'
    outfile = 'DFSET'
    disps, force = get_dfset(directory, offset_xml=offset_xml, outfile=outfile)
..


