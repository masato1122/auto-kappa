Auto-KAPPA
============

Version 0.1
---------------

Auto-KAPPA is a software for automated calculations of thermal conductivity.
The first version works with VASP and ALAMODE while, in furture versions,
other software such as QuantumEspress, Phono3py, and ShenBTE, are available,
which enables an easy comparison of results obtained with different software and
use of different methods such as SCPH method, 4-phonon scattering, iterative or direct solution of BTE, etc.

Requirements
-------------

* VASP
* Alamode >= 1.4.x
* Phonopy
* ASE
* Pymatgen
* Spglib
* Custodian


References
-----------

ALAMODE:

- T. Tadano, Y. Gohda, and S. Tsuneyuki, J. Phys.: Condens. Matter 26, 225402 (2014).


To Do
--------

- Plot kappa (kappa-fd-1, kappa-lasso-1, ...)
- Get log for kappa

