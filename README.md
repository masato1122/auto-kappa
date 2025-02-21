Auto-kappa
============

Version 0.3
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

Please follow the these steps to view the manual.

1. git clone https://github.com/masato1122/auto-kappa.git
2. cd ./auto-kappa
3. sh install.sh

How to Use (minimum)
---------------------

You can perform a simple calculation by following the steps below. 
Please refer to the manual for details.

1. Set the ``VASP_PP_PATH`` variable to allow ASE read pseudopotential files of VASP: 
([maunal HP for ASE][https://wiki.fysik.dtu.dk/ase/ase/calculators/vasp.html#pseudopotentials]).

ASE reads pseudopotential files in ``${VASP_PP_PATH}/potpaw_PBE/{element name}``.

2. Prepare a structure file: for example, ``POSCAR.Si``
3. ``akrun --file_structure POSCAR.Si --material_name Si``.

Citation
---------

Please cite the following paper, as well as papers in the references, when using Auto-kappa.

- M. Ohnishi, et al., ... (2025).

References
-----------

- **ALAMODE:** T. Tadano, Y. Gohda, and S. Tsuneyuki, J. Phys.: Condens. Matter 26, 225402 (2014).

- **ALAMODE (SCP):** T. Tadano and S. Tsuneyuki, Phys. Rev. B 92, 054301 (2015).

- **Phonopy:** A. Togo and I. Tanaka, Scr. Mater., 108, 1-5 (2015).

- **Spglib:** A. Togo and I. Tanaka, arXiv:1808.01590 (2018).

- **Pymatgen** and **Custodian:** S. P. Ong, et al., Comp. Mater. Sci. 68, 314-319 (2013).

- **ASE:** A. H. Larsen, et al., J. Phys.: Cond. Matter 29, 273002 (2017).

