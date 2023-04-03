===========================
Input Parameters of akrun
===========================

You can see the input parameters of ``akrun`` command with the following command.

.. code-block:: shell
    
    akrun -h

You can find the following commands.


Parameters I
=============

You have to modify the following parameters for each material.
[***] denotes the default value.

.. glossary::

    -d, --directory
        Directory of Phonondb in which different files such as POSCAR-unitcell, phonopy.conf, 
        KPOINTS-{relax, nac, force}, etc. are located. ``--directory`` or ``--file_structure``
        must be given. When both of them are given, ``--directory`` takes priority over 
        ``--file_structure``.

    --file_structure
        Structure file name. Different formats such as POSCAR and cif which can be read by 
        ``from_read`` module in ``Pymatgen.core.structure``.
    
    --material_name
        Material name which is used as the output directory [out].

Parameters II
===============
    
You need to modify the following parameters depending on your environment.

.. glosarry::
    
    -n, --ncores
        Number of cores in the node used for the calculation [2]
    
    --mpirun
        MPI command, for example, ``mpirun``, ``mpiexec``, etc. [mpirun]
    
    --command_vasp
        command to run VASP [vasp]

    --command_anphon
        command to run anphon [anphon]

    --command_alm
        command to run alm [alm]
    




