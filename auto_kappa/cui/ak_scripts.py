#!/home/ohnishi/.conda/envs/alm/bin/python -u
import numpy as np

from ..io.phonondb import Phonondb
from ..apdb import ApdbVasp
from ..alamode.almcalc import AlmCalc
from .. import output_directories
from .ak_parser import get_parser

def start_autokappa():
    """ Print the logo.
    Font: stick letters
    """
    print("\n\
                  ___  __                __   __       \n\
         /\  |  |  |  /  \ __ |__/  /\  |__) |__)  /\  \n\
        /~~\ \__/  |  \__/    |  \ /~~\ |    |    /~~\ \n\
        \n")

def end_autokappa():
    print("\n\
         ___       __  \n\
        |__  |\ | |  \ \n\
        |___ | \| |__/ \n\
        \n")

def main():
    
    options = get_parser()

    start_autokappa()
    
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
            'mpirun': options.mpirun, 
            'nprocs': options.ncores, 
            'nthreads': 1, 
            'vasp': options.command_vasp,
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
            'mpirun': options.mpirun, 'nprocs': 1, 
            'nthreads': options.ncores, 
            'anphon': options.command_anphon,
            'alm': options.command_alm,
            }
    almcalc = AlmCalc(
            apdb.primitive,
            base_directory=options.mpid,
            restart=options.restart,
            primitive_matrix=pmat,
            scell_matrix=smat,
            cutoff2=-1, cutoff3=options.cutoff3, 
            mag=0.01,
            nac=phdb.nac,
            command=command,
            ncores=options.ncores,
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
    almcalc.write_alamode_input(propt='band')
    almcalc.run_alamode(propt='band', force=True)
    
    ### calculate DOS
    almcalc.write_alamode_input(propt='dos')
    almcalc.run_alamode(propt='dos')
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
    
    exit()

    ### calculate anharmonic force constants
    if almcalc.lasso:
        from auto_alamode2.io.vasp import get_dfset
        for propt in ['cv', 'lasso']:
            almcalc.write_alamode_input(propt=propt)
            almcalc.run_alamode(propt)
    else: 
        almcalc.calc_anharm_force_constants()
    
    ### calculate kappa
    almcalc.write_alamode_input(propt='kappa', kpts=[15,15,15])
    almcalc.run_alamode(propt='kappa')
    almcalc.plot_kappa()
    
    end_autokappa()

