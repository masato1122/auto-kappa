===========================
Input Parameters of akrun
===========================

You can see descriptions for the input parameters of ``akrun`` with the following command.

.. code-block:: shell
    
    akrun -h

You can find the options below.
[\*\*\*] denotes the default value.


Parameters to be modified for each material
============================================

.. glossary::

    -d DIRECTORY, --directory=DIRECTORY
        Directory of Phonondb in which different files such as POSCAR-unitcell, phonopy.conf, 
        KPOINTS-{relax, nac, force}, etc. are located. ``--directory`` or ``--file_structure``
        must be given. When both of them are given, ``--directory`` takes priority over 
        ``--file_structure``.

    --file_structure=FILE_STRUCTURE
        Structure file name. Different formats such as POSCAR and cif which can be read by 
        ``Structure.from_read`` module in ``Pymatgen.core.structure``.

    --material_name=MATERIAL_NAME
        Material name which is used as the output directory [out].

        
Parameters depending on the environment
========================================

.. glossary::
    
    -n NCORES, --ncores=NCORES
        Number of cores used for the calculation [2]
    
    --anphon_para=ANPHON_PARA
        Parallel mode of anphon: "mpi" for MPI, or "omp" for OpenMP [omp]
        
    --mpirun=MPIRUN
        MPI command, for example, ``mpirun``, ``mpiexec``, etc. [mpirun]
    
    --command_vasp=COMMAND_VASP
        Command to run VASP [vasp]

    --command_anphon=COMMAND_ANPHON
        Command to run anphon [anphon]

    --command_alm=COMMAND_ALM
        Command to run alm [alm]
    
 
Parameters for the calculation conditions
==========================================

.. glossary::
    
    --cutoff3=CUTOFF3     
        Cutoff length for cubic force constants (FCs) with the unit of angstrom [4.3]
    
    --nmax_suggest=NMAX_SUGGEST
        Maximum number of suggested patterns for cubic FCs [100]
        
    --frac_nrandom=FRAC_NRANDOM
        ``Npattern * Natom / Nfc3``, where ``Npattern`` is the number of generated random 
        displacement patterns, ``Natom`` is the number of atoms in a supercell, and ``Nfc3``
        is the number of cubic FCs [1.0]
    
    --random_disp_temperature=RANDOM_DISP_TEMPERATURE
        Temperature for random displacement [500]. Atoms are displaced randomly based the 
        normal coordinates for the LASSO calculation. Tools implemented in Alamode are used.
    
    --magnitude2=MAGNITUDE2
        Magnitude of random displacement for cubic FCs [0.03]
    
    --negative_freq=NEGATIVE_FREQ
        Threshold of negative frequency [-0.001]. If the minimum value of frequencies in 
        phonon dispersion or phonon DOS, the calculation is stopped.
    
    --volume_relaxation=VOLUME_RELAXATION
        Relaxation with the relationship between volume and energy. (0.off or 1.on) [0]
        If it is "1", the relaxation is conducted by calculating the total energy with 
        different volumes and fitting to Birch-Murnaghan equation of states.
    
    --relaxed_cell=RELAXED_CELL
        Cell type used for the relaxation calculation [None] (primitive, unit/conventional, or
        supercell). For the restart calculation, the same type as the previous calculation 
        is used while, for the new calculation, the conventional cell is used. It is 
        recommended to use the conventional cell.


Parameters for k-mesh densities and size of supercells
========================================================

.. glossary::
    
    --k_length=K_LENGTH
        Length to determine k-mesh. [20]
        The way in VASP manual is used.
    
    --max_natoms=MAX_NATOMS
        Maximum limit of the number of atoms in the supercell for harmonic FCs [150]

    --max_natoms3=MAX_NATOMS3
        Maximum limit of the number of atoms in the supercell for cubic FCs [None]
    

Parameters that may be modified for test calculations
======================================================

.. glossary::

    --restart=RESTART
        The calculation will restart (1) or will be overwritten if the output directory exsits (0). [1]
        When the directory does not exist, a new calculation will be conducted.
        
    --verbose=VERBOSITY
        Verbosity (0 or 1) [0].

    --ignore_log=IGNORE_LOG
        Ignore log (1) or not (0) [0]. If it is "0", some Alamode calculations will be performed even if 
        they have been already done.
