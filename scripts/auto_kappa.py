#!/usr/bin/env python
import numpy as np
from optparse import OptionParser

from auto_alamode2.cui.scripts import start_aa2
from auto_alamode2.io.phonondb import Phonondb
from auto_alamode2.apdb import ApdbVasp
from auto_alamode2.almcalc import AlmCalc
from auto_alamode2 import output_directories

start_aa2()

def main(options):
    
    ### Read data of phonondb
    ### phonondb is used to obtain structures (primitive, unit, and super cells)
    ### and k-points.
    phdb = Phonondb(options.directory, mpid=options.mpid)

    pmat = phdb.primitive_matrix
    smat = phdb.scell_matrix
    
    ### Set output directoreis
    out_dirs = {}
    for k1 in output_directories.keys():
        values1 = output_directories[k1]
        if type(values1) == str:
            out_dirs[k1] = './' + options.mpid + '/' + values1
        else:
            out_dirs[k1] = {}
            for k2 in values1.keys():
                values2 = values1[k2]
                out_dirs[k1][k2] = './' + options.mpid + '/' + values2
    
    ### command to run VASP jobs
    command = {
            'mpirun': 'mpirun', 'nprocs': options.ncores, 
            'nthreads': 1, 'vasp': 'vasp'
            }
    
    ### Set ApdbVasp object
    apdb = ApdbVasp(
            phdb.get_unitcell(format='ase'), 
            primitive_matrix=pmat,
            scell_matrix=smat,
            command=command)
    
    ### Relaxation calculation
    mode = 'relax'
    apdb.run_vasp(mode, 
            out_dirs[mode],
            phdb.get_kpoints(mode=mode).kpts[0],
            print_params=True
            )
    
    ### Born effective charge
    if phdb.nac == 1:
        mode = 'nac'
        apdb.run_vasp(mode,
                out_dirs[mode],
                phdb.get_kpoints(mode=mode).kpts[0],
                print_params=True
                )
    
    ### Set AlmCalc
    command = {
            'mpirun': 'mpirun', 'nprocs': 1, 
            'nthreads': options.ncores, 'anphon': 'anphon'
            }
    almcalc = AlmCalc(
            apdb.primitive,
            mpid=options.mpid,
            primitive_matrix=pmat,
            scell_matrix=smat,
            cutoff2=-1, cutoff3=options.cutoff3, 
            nbody=[2,3,3,2], mag=0.01,
            nac=phdb.nac,
            command=command,
            verbosity=options.verbosity
            )
    
    ### ASE calculator for forces
    mode = 'force'
    calc_force = apdb.get_calculator(mode,
            out_dirs['harm'][mode],
            phdb.get_kpoints(mode=mode).kpts[0]
            )
    
    ### calculate forces for harmonic FCs
    almcalc.calc_forces(1, calc_force)
    almcalc.calc_harmonic_force_constants()
    
    ### calculate band 
    almcalc.write_anphon_input(propt='band')
    almcalc.run_anphon(propt='band', force=True)
    
    ### calculate DOS
    almcalc.write_anphon_input(propt='dos')
    almcalc.run_anphon(propt='dos')
    almcalc.plot_bandos()
    
    ### Check negative frequency
    if almcalc.frequency_range[0] < options.negative_freq:
        print("")
        print(" Negative eigenvalues were found. Stop the calculation.")
        print(" Minimum frequency : %.2f" % (almcalc.frequency_range[0]))
        print("")
        exit()
    
    ### calculate forces for cubic FCs
    mode = 'force'
    almcalc.calc_forces(
            2, calc_force, 
            nmax_suggest=options.nmax_suggest,
            frac_nrandom=options.frac_nrandom,
            temperature=options.random_disp_temperature,
            )
    
    almcalc.calc_anharm_force_constants()
    
    ### calculate kappa
    almcalc.write_anphon_input(propt='kappa', kpts=[15,15,15])
    almcalc.run_anphon(propt='kappa')
    almcalc.plot_kappa()

if __name__ == '__main__':
    
    parser = OptionParser()

    ### parameters which need to be modified for each.
    parser.add_option("-d", "--directory", dest="directory", type="string",
            default="../mp-149", help="directory of phonondb")
    
    parser.add_option("--mpid", dest="mpid", type="string",
            default="mp-149", 
            help="material ID, which is used for the name of directory [mp-149]")
    
    parser.add_option("-n", "--ncores", dest="ncores", type="int",
            default=2, help="ncores [2]")
    
    parser.add_option("--verbosity", dest="verbosity", type="int",
            default=1, help="verbosity [0]")
    
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
    
    ### parameters which do not need to be changed.
    parser.add_option("--nagative_freq", dest="negative_freq", type="float",
            default=-0.001, help="threshold of negative frequency [-0.001]")
            
    parser.add_option("--random_disp_temperature", 
            dest="random_disp_temperature", type="float",
            default=500., 
            help="temperature for random displacement [500]")
            
    (options, args) = parser.parse_args()

    main(options)

