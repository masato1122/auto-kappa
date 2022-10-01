
from optparse import OptionParser

def get_parser():
    
    parser = OptionParser()

    ### parameters which need to be modified for each.
    parser.add_option("-d", "--directory", dest="directory", type="string",
            default="../mp-149", help="directory of phonondb")
    
    parser.add_option("--material_name", dest="material_name", type="string",
            default="mp-149", 
            help="material name such as mateiral ID which is used as "\
                    "the name of directory [mp-149]")
    
    parser.add_option("--restart", dest="restart", type="int",
            default=1,
            help="The calculation will restart (1) or will NOT restart (0) "\
                    "when the directory exsits. [1]")
            
    parser.add_option("--mpirun", dest="mpirun", type="string",
            default="mpirun", help="MPI command [mpirun]")
    
    parser.add_option("-n", "--ncores", dest="ncores", type="int",
            default=2, help="ncores [2]")
    
    parser.add_option("--verbosity", dest="verbosity", type="int",
            default=1, help="verbosity [0]")
    
    parser.add_option("--neglect_log", dest="neglect_log", type="int",
            default=1, help="neglect log (1) or not (0) [1]")
    
    ### parameters which may not need to be changed.
    parser.add_option("--cutoff3", dest="cutoff3", type="float",
            default=4.3, help="cutoff3, unit=Ang [4.3]")
    
    parser.add_option("--nmax_suggest", 
            dest="nmax_suggest", type="int", default=200, 
            help="Maximum number of suggested patterns for cubic FCs [200]")

    parser.add_option("--frac_nrandom", 
            dest="frac_nrandom", type="float", default=0.02,
            help="Ratio of the number of generated patterns with random "\
                    "displacement to the number for the suggested patterns "
                    "with ALM [0.02]")
    
    parser.add_option("--command_vasp", 
            dest="command_vasp", type="string", default="vasp", 
            help="command to run vasp [vasp]")
    
    parser.add_option("--command_anphon", 
            dest="command_anphon", type="string", default="anphon", 
            help="command to run anphon [anphon]")
    
    parser.add_option("--command_alm", 
            dest="command_alm", type="string", default="alm", 
            help="command to run alm [alm]")
    
    ### parameters which do not need to be changed.
    parser.add_option("--nagative_freq", dest="negative_freq", type="float",
            default=-0.001, help="threshold of negative frequency [-0.001]")
            
    parser.add_option("--random_disp_temperature", 
            dest="random_disp_temperature", type="float",
            default=500., 
            help="temperature for random displacement [500]")
    
    (options, args) = parser.parse_args()
    
    return options


