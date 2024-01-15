#
# ak_parser.py
#
# Parser of akrun command
#
# Copyright (c) 2022 Masato Ohnishi
#
# This file is distributed under the terms of the MIT license.
# Please see the file 'LICENCE.txt' in the root directory
# or http://opensource.org/licenses/mit-license.php for information.
#
from optparse import OptionParser

def get_parser():
    
    parser = OptionParser()
    
    ### parameters that need to be modified for each calculation
    parser.add_option("-d", "--directory", dest="directory", type="string",
            default=None, help=""
            "Directory name for data of Phonondb "
            "where POSCAR-unitcell, phonopy.conf, KPOINTS-**, etc. "
            "are located. "
            "``--directory`` or ``--file_structure`` must be given. "
            "When both of them are given, ``--directory`` takes "
            "priority over ``--file_structure``."
            )
     
    parser.add_option("--file_structure", dest="file_structure", type="string",
            default=None, 
            help="Structure file name. Different kinds of format such as "
            "POSCAR and cif file can be read by ``Structure.from_file`` module "
            "in ``Pymatgen.core.structure``.")
    
    parser.add_option("--material_name", dest="material_name", type="string",
            default="./out", 
            help="Material name which is used as the output directory [./out]")
    
    ### Parameters that need to be modified depending on the environment
    parser.add_option("-n", "--ncores", dest="ncores", type="int",
            default=2, help="Number of cores used for the calculation [2]")
    
    parser.add_option("--anphon_para", 
            dest="anphon_para", type="string", default="mpi", 
            help="Parallel mode for \"anphon\". "
            "Input value should be \"mpi\" for MPI and "
            "\"omp\" for OpenMP [mpi]. "
            "While \"mpi\", which is the default value, is "
            "faster for most cases, "
            "if the system is complex and requires a large memory, "
            "it is recommended to use \"omp\". ")
     
    parser.add_option("--mpirun", dest="mpirun", type="string",
            default="mpirun", help="MPI command [mpirun]")
    
    parser.add_option("--command_vasp", 
            dest="command_vasp", type="string", default="vasp", 
            help="Command to run vasp [vasp]")
    
    parser.add_option("--command_anphon", 
            dest="command_anphon", type="string", default="anphon", 
            help="Command to run anphon [anphon]")
    
    parser.add_option("--command_alm", 
            dest="command_alm", type="string", default="alm", 
            help="Command to run alm [alm]")
    
    parser.add_option("--nonanalytic", 
            dest="nonanalytic", type="int", default=2, 
            help="NONANALYTIC tag for Anphon calculation [2]. "
            "While the default value is 2, the value is "
            "modified to 0 when NAC is not considered.")
    
    ### Parameters for the calculation condictions
    parser.add_option("--cutoff3", dest="cutoff3", type="float",
            default=4.3, 
            help="Cutoff length for cubic force constants with the unit of "
            "angstrom [4.3]")
    
    parser.add_option("--nmax_suggest", 
            dest="nmax_suggest", type="int", default=100, 
            help="Maximum number of suggested patterns for cubic FCs [100]")
    
    parser.add_option("--frac_nrandom", 
            dest="frac_nrandom", type="float", default=1.,
            help="``Npattern * Natom / Nfc3``, where Npattern is the number "
            "of generated random displacement patterns, Natom is the number "
            "of atoms in a supercell, and Nfc3 is the number of FC3 [1.0]",
            )
    
    parser.add_option("--random_disp_temperature", 
            dest="random_disp_temperature", type="float",
            default=500., 
            help="temperature for random displacement [500]")
    
    parser.add_option("--magnitude2", 
            dest="magnitude2", type="float", default=0.03, 
            help="magnitude of random displacement for cubic FCs [0.03]")
    
    parser.add_option("--negative_freq", dest="negative_freq", type="float",
            default=-0.001, help="threshold of negative frequency [-0.001]")
            
    parser.add_option("--volume_relaxation", 
            dest="volume_relaxation", type="int", default=0,
            help="strict relaxation with the relation between energy and "
            "volume (0.off or 1.on) [0]")
    
    parser.add_option("--relaxed_cell", 
            dest="relaxed_cell", type="string", default=None,
            help="Cell type used for the relaxation calculation [None]. "
            "For the restarted calculation, the same type "
            "as that for the previous calculation is used "
            "while the conventional cell "
            "is used for the new calculation. "
            "The value should be primitive (p), conventional (c), "
            "or unitcell (u), where \"c\" and \"u\" are the same."
            )
    
    ### Parameters to determine k-mesh densities and size of supercells
    parser.add_option("--k_length", dest="k_length", type="float", 
            default=20, help=
            "Length to determine k-mesh [20]. "
            "The following equation is used to determine the k-mesh; "
            "N = max(1, int(k_length * |a|^* + 0.5))."
            )
    
    ### options for supercell
    parser.add_option("--max_natoms", dest="max_natoms", type="int", 
            default=150, 
            help="Initial maximum limit of the number of atoms in the "
            "supercells [150]. When negative frequencies exist "
            "at kpoints except for the commensurate points, "
            "the maximum limit for harmonic FCs "
            "will be increased with the step of "
            "\"delta_max_natoms\". Note that the maximum limit for "
            "cubic FCs is not changed during the simlation. ")
    
    ### three options for larger supercell
    parser.add_option("--analyze_with_larger_supercell", 
            dest="analyze_with_largersc", type="int", default=0,
            help="Analyze harmonic properties with larger supercells if "
            "negative frequencies exist. (0.no, 1.yes) [0]"
            )
    
    parser.add_option("--delta_max_natoms", dest="delta_max_natoms", 
            type="int", default=50, 
            help="Increasing interval of the maximum limit of the number of "
            "atoms in the supercell to calculate FC2 [50]. "
            "This options is used only when ``--analyze_with_largersc`` == 1.")
    
    parser.add_option("--max_loop_for_largesc", dest="max_loop_for_largesc", 
            type="int", default=1, 
            help="Number of the maximum loop for larger supercell [1]."
            "This options is used only when ``--analyze_with_largersc`` == 1.")
    
    ### parameters for NSW
    parser.add_option("--nsw_params", dest="nsw_params", 
            type="string", default="100:10:20", 
            help="Parameters which determine NSW for relaxation calculations "\
                    "[100:10:20]. "\
                    "\"{nsw_init}:{nsw_diff}:{nsw_min}\": " \
                    "NSW = min(``nsw_min``, "\
                    "``nsw_init`` - ``nsw_diff`` * ``num_errors``), where "\
                    "``num_errors`` is the number of errors.")
    
    #### parameters for amin
    #parser.add_option("--amin", dest="amin", 
    #        type="float", default=None, 
    #        help=("AMIN parameter for VASP [None]: If the length of a lattice "
    #            "vector exceeds 5 nm, AMIN of the given value is set for the "
    #            "VASP job.")
    #        )
    
    #### calculate potential energy sruface
    #parser.add_option("--pes", dest="pes",
    #        type="int", default=0, 
    #        help="Calculate potential energy surface (PES) "
    #        "with respect to a phonon mode with a negative frequency. [0] "
    #        "PES is "
    #        "0. not calculated, "
    #        "1. calculated only for larger supercells, or "
    #        "2. calculated for both of small and larger supercells. "
    #        "A representative k-point is chosen to calculate the PES: "
    #        "Gamma point, commensurate point, or grid-point for DOS."
    #        )
    
    #parser.add_option("--max_natoms3", dest="max_natoms3", type="int", 
    #        default=None, 
    #        help="This options is invalid! PLEASE DO NOT USE this option."\
    #                "Maximum limit of the number of atoms in the supercell for "\
    #                " FC3 [None].")
    
    ### Parameters that need to be changed for test calculations 
    parser.add_option("--restart", dest="restart", type="int",
            default=1,
            help="The calculation will restart (1) or will NOT restart (0) "
            "when the directory exsits. [1]")
            
    parser.add_option("--verbosity", dest="verbosity", type="int",
            default=1, help="verbosity [0]")
    
    parser.add_option("--neglect_log", dest="neglect_log", type="int",
            default=0, help="neglect log (1) or not (0) [0]")
    
    parser.add_option("--harmonic_only", dest="harmonic_only", type="int",
            default=0, 
            help="Calculate harmonic properties only (0.No, 1.Yes) [0]")
    
    parser.add_option("--max_relax_error", dest="max_relax_error", type="int",
            default=500, 
            help="Maximum number of errors for relaxation calculations with "
            "VASP. Set this option if the number of error is too many. [500]")
    
    ### author
    parser.add_option("--authors", dest="authors", type="string",
            default=None, help="authors' name (A^1, B^1, C^2)")
    parser.add_option("--affiliations", dest="affiliations", type="string",
            default=None, help="affiliation (1. Univ. 1, 2. Univ. 2)")
    
    (options, args) = parser.parse_args()
    
    return options

