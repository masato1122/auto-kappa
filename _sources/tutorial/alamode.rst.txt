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

