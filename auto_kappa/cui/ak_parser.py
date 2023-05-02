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
    line = "Directory of Phonondb where POSCAR-unitcell, phonopy.conf, "
    line += "KPOINTS-**, etc. are located. --directory or --file_structure "
    line += "must be given. When both of them are given, --directory takes "
    line += "priority over --file_structure."
    parser.add_option("-d", "--directory", dest="directory", type="string",
            default=None, help=line)
     
    parser.add_option("--file_structure", dest="file_structure", type="string",
            default=None, 
            help="Structure file name. Different kinds of format such as "\
                    "POSCAR and cif file are accepted.")

    parser.add_option("--material_name", dest="material_name", type="string",
            default="./out", 
            help="Material name which is used as the output directory [out]")
    
    ### Parameters that need to be modified depending on the environment
    parser.add_option("-n", "--ncores", dest="ncores", type="int",
            default=2, help="Number of cores used for the calculation [2]")
    
    parser.add_option("--anphon_para", 
            dest="anphon_para", type="string", default="mpi", 
            help="Parallel mode for \"anphon\". " \
                    "Input value should be \"mpi\" for MPI and \"omp\" for OpenMP. "\
                    "While \"mpi\", which is the default value, is "\
                    "faster for most cases, "\
                    "if the system is complex and requires a large memory, "\
                    "use \"omp\". "\
                    "It is also recommended to check if \"anphon\" is running "\
                    "properly with MPI before starting "\
                    "a high-throughput calculation. [mpi]")
     
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
            help="NONANALYTIC tag for Anphon calculation. "\
                    "The default value is 0 if NAC is not considered. "\
                    "while this is set to be 2 if NAC is considered.")
    
    ### Parameters for the calculation condictions
    parser.add_option("--cutoff3", dest="cutoff3", type="float",
            default=4.3, 
            help="Cutoff length for cubic force constants with the unit of "\
                    "angstrom: [4.3]")
    
    parser.add_option("--nmax_suggest", 
            dest="nmax_suggest", type="int", default=100, 
            help="Maximum number of suggested patterns for cubic FCs [100]")
    
    parser.add_option("--frac_nrandom", 
            dest="frac_nrandom", type="float", default=1.,
            help="``Npattern * Natom / Nfc3``, where Npattern is the number "\
            "of generated random displacement patterns, Natom is the number "\
            "of atoms in a supercell, and Nfc3 is the number of FC3 [1.0]",
            )
    
    parser.add_option("--random_disp_temperature", 
            dest="random_disp_temperature", type="float",
            default=500., 
            help="temperature for random displacement [500]")
    
    parser.add_option("--magnitude2", 
            dest="magnitude2", type="float", default=0.03, 
            help="magnitude of random displacement for FC3 [0.03]")
    
    parser.add_option("--negative_freq", dest="negative_freq", type="float",
            default=-0.001, help="threshold of negative frequency [-0.001]")
            
    parser.add_option("--volume_relaxation", 
            dest="volume_relaxation", type="int", default=0,
            help="strict relaxation with the relation between energy and "\
                    "volume (0.off or 1.on) [0]")
    
    parser.add_option("--relaxed_cell", 
            dest="relaxed_cell", type="string", default=None,
            help="Cell type used for the relaxation calculation [None]. "\
                    "For the restarted calculation, the same type "\
                    "as that for the previous calculation is used "\
                    "while, for the new calculation, the conventional cell "
                    "is used."
            )
    
    ### Parameters to determine k-mesh densities and size of supercells
    parser.add_option("--k_length", dest="k_length", type="float", 
            default=20, help="Length to determine k-mesh. [20]")
    
    parser.add_option("--max_natoms", dest="max_natoms", type="int", 
            default=150, 
            help="Maximum limit of the number of atoms in the supercell "\
                    "for FC2 [150]. For the automation calculation, it is NOT"\
                    "recommended to modify this option.")
    
    parser.add_option("--max_natoms3", dest="max_natoms3", type="int", 
            default=None, 
            help="This options is invalid! PLEASE DO NOT USE this option."\
                    "Maximum limit of the number of atoms in the supercell for "\
                    " FC3 [None].")
    
    ### Parameters that need to be changed for test calculations 
    parser.add_option("--restart", dest="restart", type="int",
            default=1,
            help="The calculation will restart (1) or will NOT restart (0) "\
                    "when the directory exsits. [1]")
            
    parser.add_option("--verbosity", dest="verbosity", type="int",
            default=1, help="verbosity [0]")
    
    parser.add_option("--neglect_log", dest="neglect_log", type="int",
            default=0, help="neglect log (1) or not (0) [0]")
    
    parser.add_option("--harmonic_only", dest="harmonic_only", type="int",
            default=0, 
            help="calculate harmonic properties only (0.No, 1.Yes) [0]")
    
    (options, args) = parser.parse_args()
    
    return options

