========================================
Input Scripts for ALAMODE Calculations
========================================

.. contents:: Contents
   :local:
   :depth: 1


Preparation of ALAMODE input files
==================================

Force constants
---------------

To suggest structures for a brute force method:

.. literalinclude:: ./examples/mksuggest.py

To calculate force constants with calculated forces:

.. literalinclude:: ./examples/mkoptimize.py

Phonon dispersion
------------------

.. literalinclude:: ./examples/mkband.py

.. Eigenvalues at commensurate points
.. ----------------------------------
.. 
.. .. include:: ./examples/mkcommensurate.py


Thermal conductivity
--------------------

.. literalinclude:: ./examples/mkkappa.py


BORNINFO file for ALM
======================

from vasprun.xml
------------------

.. code-block:: python
    
    from auto_kappa.io.born import BORNINFO
    file_vasp = 'vasprun.xml'
    file_born = 'FC2.xml'
    born = BORNINFO(file_vasp, file_born=file_born)
    born.write(outfile='BORNINFO')

