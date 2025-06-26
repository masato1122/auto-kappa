Auto-kappa
============

Version 0.4.0
---------------

Auto-kappa is an automation tool for calculating anharmonic phonon properties
—including thermal conductivity and mode-dependent phonon lifetimes—using VASP and ALAMODE.


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

Follow these steps to install the package:

1. git clone https://github.com/masato1122/auto-kappa.git
2. cd ./auto-kappa
3. sh install.sh

After installation, ensure that the ``akrun`` command is available.
You can view a description of the input parameters by running ``akrun -h``.

Preparation
--------------

You can perform a simple calculation following the steps below. 
Please refer to the manual for details.

1. Set the ``VASP_PP_PATH`` environment variable so that ASE can locate VASP pseudopotential files:
([Pseudopotential with ASE](https://wiki.fysik.dtu.dk/ase/ase/calculators/vasp.html#pseudopotentials))

ASE expects the pseudopotential files in: ``${VASP_PP_PATH}/potpaw_PBE/{element name}``.

2. Prepare a structure file, e.g., ``POSCAR.Si``
3. Run the command: ``akrun --file_structure POSCAR.Si --material_name Si``.

Some options
-------------

You can view a description of the options using ``akrun -h``. 
Frequently used commands are listed below.

- file_structure: Structure file name (Different formats are acceptable such as POSCAR, cif, etc.)

- material_name: Name of the output directory

- nprocs: Number of process for the calculation [2]

- mpirun: MPI command [mpirun]

- command_{VASP/alm/anphon}: Command to run VASP, alm, and anphon [vasp, alm, anphon]

- volume_relaxation: Perform relaxation calculations using the Birch-Murnaghan equation of state

- analyze_with_larger_supercell : Use a larger supercell when imaginary frequencies appear


Citation
---------

If you use auto-kappa, please cite the following paper, along with any related papers listed in the references:

- M. Ohnishi et al., "Database and deep-learning scalability of anharmonic phonon properties by automated brute-force first-principles calculations", 
[arXiv:2504.21245](https://arxiv.org/abs/2504.21245) (2025).

References
-----------

- **ALAMODE:** T. Tadano, Y. Gohda, and S. Tsuneyuki, J. Phys.: Condens. Matter 26, 225402 (2014).

- **ALAMODE (SCP):** T. Tadano and S. Tsuneyuki, Phys. Rev. B 92, 054301 (2015).

- **Phonopy:** A. Togo and I. Tanaka, Scr. Mater., 108, 1-5 (2015).

- **Spglib:** A. Togo, K. Shinohara, and I. Tanaka, Sci. technol. adv. material, Meth. 4, 1 (2025).

- **Pymatgen** and **Custodian:** S. P. Ong et al., Comp. Mater. Sci. 68, 314-319 (2013).

- **ASE:** A. H. Larsen et al., J. Phys.: Cond. Matter 29, 273002 (2017).

