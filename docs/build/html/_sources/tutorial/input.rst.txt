==============================
How to Prepare Input Scripts
==============================

You can find example scripts in example directory of aiida-alamode package.

ALM Input Script
===========================

.. code-block:: python

    from aiida_alamode.io import AlmInput

    filename = 'POSCAR-supercell'
    dfset = 'dfset.txt'
    
    ### mode "suggest"
    alminp = AlmInput.from_structure_file(
            filename,
            mode='suggest',
            norder=2
            )
    alminp.to_file()
    #alminp.to_file(filename='suggest.in')
    
    ### mode "optimize"
    ### Note: Cutoff should be given for high-order IFCs.
    alminp = AlmInput.from_structure_file(
            filename,
            mode='optimize',
            norder=2,
            cutoff={'Si-Si': [None, None]},
            dfset='dfset.txt'
            )
    alminp.to_file()


BORNINFO for ALM
===========================

VASP
----------

.. code-block:: python
    
    from aiida_alamode.io.vasp import write_born_info
    filename = 'vasprun.xml'
    write_born_info(filename, outfile='BORNINFO')
    
    ## old version
    #from aiida_alamode.io.vasp import Vasprun
    #vasprun = Vasprun('../vasprun.xml')
    #vasprun.write_born_info(filename='BORNINFO')

QuantumEspresso
-----------------

.. code-block:: python
    
    from aiida_alamode.io.qe import ***




ANPHON Input Script
=========================
 
AnphonInput class can be used with a similar manner as that for AlmInput.

.. code-block:: python
    
    from aiida_alamode.io import AnphonInput
    
    filename = 'POSCAR-primitive'
    fcsxml = 'IFCs_anharm.xml'
    nac = 1; borninfo = 'BORNINFO'
    #nac = 0; borninfo = None
    
    nk = 10
    deltak = 0.2
    
    ### to calculate band
    anpinp = AnphonInput.from_structure_file(
            filename,
            mode='phonons',
            kpmode=1,
            fcsxml=fcsxml,
            nonanalytic=nac, borninfo=borninfo
            )
    anpinp.set_kpoint(deltak=0.01)
    anpinp.to_file(filename='band.in')
    
    ### to calculate DOS
    anpinp = AnphonInput.from_structure_file(
            filename,
            mode='phonons',
            kpmode=2,
            fcsxml=fcsxml,
            nonanalytic=nac, borninfo=borninfo
            )
    #anpinp.set_kpoint(deltak=deltak)
    anpinp.update({'kpts':[nk, nk, nk]})
    anpinp.to_file(filename='dos.in')

    ### to calculate thermal conductivity
    anpinp = AnphonInput.from_structure_file(
            filename,
            mode='RTA',
            kpmode=2, kpts=[nk, nk, nk],
            fcsxml=fcsxml,
            nonanalytic=nac, borninfo=borninfo,
            isotope=2
            )
    anpinp.to_file(filename='RTA.in')


