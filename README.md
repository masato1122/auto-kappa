Auto-kappa
============

Version 0.4.0
---------------

Auto-kappa is an automation software for calculating anharmonic phonon properties, 
including thermal conductivity, using VASP and ALAMODE.

Requirements
-------------

* VASP >= 6.3.x
* Alamode >= 1.4.x
* Phonopy
* ASE
* Pymatgen
* Spglib
* Custodian

Installation
-------------

Please follow these steps to install the package.

1. git clone https://github.com/masato1122/auto-kappa.git
2. cd ./auto-kappa
3. sh install.sh

How to Use (minimum)
---------------------

You can perform a simple calculation by following the steps below. 
Please refer to the manual for details.

1. Set the ``VASP_PP_PATH`` variable to allow ASE read pseudopotential files of VASP: 
([Pseudopotential with ASE](https://wiki.fysik.dtu.dk/ase/ase/calculators/vasp.html#pseudopotentials))

ASE reads pseudopotential files in ``${VASP_PP_PATH}/potpaw_PBE/{element name}``.

2. Prepare a structure file: for example, ``POSCAR.Si``
3. ``akrun --file_structure POSCAR.Si --material_name Si``.

Citation
---------

Please cite the following paper, as well as related papers in the references, when using auto-kappa.

- M. Ohnishi et al., "Database and deep-learning scalability of anharmonic phonon properties by automated brute-force first-principles calculations", 
[arXiv:2504.21245](https://arxiv.org/abs/2504.21245) (2025).

References
-----------

- **ALAMODE:** T. Tadano, Y. Gohda, and S. Tsuneyuki, J. Phys.: Condens. Matter 26, 225402 (2014).

- **ALAMODE (SCP):** T. Tadano and S. Tsuneyuki, Phys. Rev. B 92, 054301 (2015).

- **Phonopy:** A. Togo and I. Tanaka, Scr. Mater., 108, 1-5 (2015).

- **Spglib:** A. Togo, K. Shinohara, and I. Tanaka, Sci. technol. adv. material, Meth. 4, 1 (2025).

- **Pymatgen** and **Custodian:** S. P. Ong, et al., Comp. Mater. Sci. 68, 314-319 (2013).

- **ASE:** A. H. Larsen, et al., J. Phys.: Cond. Matter 29, 273002 (2017).

