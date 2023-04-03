===========================
Input Scripts for ALAMODE
===========================

.. You can find example scripts in example directory of aiida-alamode package.

How to prepare input scripts for ALAMODE
=========================================

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
    
    from auto_kappa.io.vasp import write_born_info
    filename = 'vasprun.xml'
    write_born_info(filename, outfile='BORNINFO')
    

