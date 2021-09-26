=================
Tutorial
=================

How to Make Input Scripts
===========================

You can find example scripts in example directory of aiida-alamode package.


ALM Input Script
-------------------

.. code-block:: python

    from aiida_alamode.io import AlmInput

    filename = 'POSCAR-supercell'
    dfset = 'dfset.txt'
    
    ## mode "suggest"
    alminp = AlmInput.from_structure_file(
            filename,
            mode='suggest',
            norder=2
            )
    alminp.to_file()
    #alminp.to_file(filename='filename.in')
    
    ## mode "optimize"
    ## Note: Cutoff should be given for high-order IFCs.
    alminp = AlmInput.from_structure_file(
            filename,
            mode='optimize',
            norder=2,
            cutoff={'Si-Si': [None, None]},
            dfset='dfset.txt'
            )
    alminp.to_file()



ANPHON Input Script
----------------------
  
AnphonInput class can be used with a similar manner as that for AlmInput.

.. code-block:: python
    
    from aiida_alamode.io import AnphonInput
    
    filename = 'POSCAR-supercell'
    fcsxml = 'IFCs.xml'
    
    ## to calculate band
    anpinp = AnphonInput.from_structure_file(
            filename,
            mode='phonons',
            kpmode=1,
            fcsxml=fcsxml
            )
    anpinp.set_kpoint(deltak=0.01)
    anpinp.to_file()

    ## to calculate DOS
    anpinp = AnphonInput.from_structure_file(
            filename,
            mode='phonons',
            kpmode=2,
            fcsxml=fcsxml
            )
    anpinp.set_kpoint(deltak=0.1)
    #anpinp.update({'kpts':[10, 10, 10]})
    anpinp.to_file()

    ## to calculate thermal conductivity
    anpinp = AnphonInput.from_structure_file(
            filename,
            mode='RTA',
            kpmode=2,
            fcsxml=fcsxml
            )
    anpinp.set_kpoint()
    anpinp.to_file()



