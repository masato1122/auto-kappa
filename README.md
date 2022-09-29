Auto-KAPPA
============

Version 0.0
---------------

Auto-KAPPA is a software for automated calculations of thermal conductivity.
The first version works with VASP and ALAMODE while, in furture versions,
other software such as QuantumEspress, Phono3py, and ShenBTE, are available,
which enables an easy comparison of results obtained with different software and
use of different methods such as SCPH method, 4-phonon scattering, iterative or direct solution of BTE, etc.

Requirements
-------------

* VASP
* Alamode   >= 1.4.1
* Phonopy   >= 2.15.1
* ASE       >= 3.22.1
* Pymatgen  >= 2022.3.24
* Spglib    >= 1.16.5
* Custodian >= 2022.5.26


References
-----------

ALAMODE:

- T. Tadano, Y. Gohda, and S. Tsuneyuki, J. Phys.: Condens. Matter 26, 225402 (2014).

To Do
------

* Primitive cell with an atom may lead to an error of the treatment of crystal structure 
such as Born effective charge. (which one?)

* Check non-analytical term (e.g. 5504, 5637): 
Use dataset (displacements and forces) in Phonondb and calculate phonon dispersion with Alamode.
Use dataset (displacements and forces) in Apdb and calculate phonon dispersion with Phonopy.

Make phonopy_disp.yaml from DFSET of Alamode.

* How were the supercell size and k-mesh decided for Phonondb?
Number of atoms? How can k-mesh be obtained from Materials Project?

Additional
------------

* Claculate kappa with different k-points and extrapolate?

* Output analysis conditions, fitting error, etc.? (with json?)

* Compression of results

* Automated creation of the homepage with Sphinx

Done
------

* Relaxation caluclations for many materials were not finished in Chariot.

* EDIFFG is modified.

