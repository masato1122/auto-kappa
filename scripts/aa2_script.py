# -*- coding: utf-8 -*-
import numpy as np
from optparse import OptionParser

from auto_alamode2.io.phonondb import Phonondb
from auto_alamode2.apdb import ApdbVasp
from auto_alamode2.almcalc import AlmCalc
from auto_alamode2 import output_directories

def main(options):
    
    ### Read data of phonondb
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
    
    ### command to run a VASP job
    vasp_cmd = {'mpirun': 'mpirun', 'nprocs': options.nprocs, 'vasp': 'vasp'}
    
    ### Set ApdbVasp object
    apdb = ApdbVasp(
            phdb.get_unitcell(format='ase'), 
            primitive_matrix=pmat,
            scell_matrix=smat,
            command=vasp_cmd)
    
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
    almcalc = AlmCalc(
            apdb.primitive,
            mpid=options.mpid,
            primitive_matrix=pmat,
            scell_matrix=smat,
            cutoff2=-1, cutoff3=options.cutoff3, 
            nbody=[2,3], mag=0.01,
            nac=phdb.nac
            )
    
    ### ASE calculator for forces
    mode = 'force'
    calc_force = apdb.get_calculator(mode,
            out_dirs['harm'][mode],
            phdb.get_kpoints(mode=mode).kpts[0]
            )

    ### calculate forces for harmonic IFCs
    almcalc.calc_forces(1, calc_force)
    almcalc.calc_harmonic_force_constants()

    ##### calculate band and DOS
    almcalc.write_anphon_input(propt='band')
    almcalc.run_anphon(
            propt='band', force=True, 
            nprocs=1, nthreads=options.nprocs,
            )
    ##
    almcalc.write_anphon_input(propt='dos')
    almcalc.run_anphon(propt='dos',
            nprocs=1, nthreads=options.nprocs,
            )
    #almcalc.plot_bandos()
    
    ### Check negative frequency
    if almcalc.frequency_range[0] < -0.001:
        print("")
        print(" Negative eigenvalues were found. Stop the calculation.")
        print("")
        exit()
    
    ### calculate forces for cubic IFCs
    mode = 'force'
    almcalc.calc_forces(
            2, calc_force, 
            nmax_suggest=options.nmax_suggest,
            frac_nrandom=options.frac_nrandom,
            temperature=500.,
            )
    
    almcalc.calc_anharm_force_constants()
    
    ### calculate kappa
    almcalc.write_anphon_input(propt='kappa')
    almcalc.run_anphon(propt='kappa', nprocs=1, nthreads=options.nprocs)
    almcalc.plot_kappa()

if __name__ == '__main__':
    
    parser = OptionParser()

    parser.add_option("-d", "--directory", dest="directory", type="string",
            default="../mp-149", help="directory of phonondb")

    parser.add_option("--mpid", dest="mpid", type="string",
            default="mp-149", help="material ID [mp-149]")

    parser.add_option("--cutoff3", dest="cutoff3", type="float",
            default=4.3, help="cutoff3, unit=Ang [4.3]")
    
    parser.add_option("-n", "--nprocs", dest="nprocs", type="int",
            default=2, help="nprocs [2]")
            
    parser.add_option("--nmax_suggest", 
            dest="nmax_suggest", type="int", default=200, 
            help="Maximum number of suggested patterns for cubic FCs [200]")

    parser.add_option("--frac_nrandom", 
            dest="frac_nrandom", type="float", default=0.02,
            help="Ratio of the number of generated patterns with random "\
                    "displacement to the number for the suggested patterns "
                    "with ALM [0.02]")

    (options, args) = parser.parse_args()

    main(options)

