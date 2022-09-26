import numpy as np
import datetime

from ..io.phonondb import Phonondb
from ..apdb import ApdbVasp
from ..alamode.almcalc import AlamodeCalc
from .. import output_directories
from .ak_parser import get_parser

def start_autokappa():
    """ Print the logo.
    Font: stick letters
    """
    from ..version import __version__
    print("\n\
                  ___  __                __   __       \n\
         /\  |  |  |  /  \ __ |__/  /\  |__) |__)  /\  \n\
        /~~\ \__/  |  \__/    |  \ /~~\ |    |    /~~\ \n\
        \n\
        ver. %s\n\
        \n" % (__version__))
    time = datetime.datetime.now()
    print(" Start at", time.strftime("%m/%d/%Y %H:%M:%S"))
    print("")

def print_times(times):
    ttot = times['total'].seconds
    print("")
    print("")
    print(" Calculation times:")
    print(" ==================")
    print("")
    for key in times:
        if key == 'total':
            continue
        tt = float(times[key].seconds)
        percentage = 100.*tt/ttot
        print(" %15s (sec): %13.3f (%.1f%%)" % (key, tt, percentage))
    print("")

def end_autokappa():
    
    time = datetime.datetime.now()
    print("\n\
         ___       __  \n\
        |__  |\ | |  \ \n\
        |___ | \| |__/ \n\
        \n")
    print(" at", time.strftime("%m/%d/%Y %H:%M:%S"), "\n\n")

def main():
    
    times = {}
    
    options = get_parser()
    
    start_autokappa()
    
    ### Read data of phonondb
    ### phonondb is used to obtain structures (primitive, unit, and super cells)
    ### and k-points.
    phdb = Phonondb(options.directory)
    
    pmat = phdb.primitive_matrix
    smat = phdb.scell_matrix
    
    ### Set output directories
    out_dirs = {}
    for k1 in output_directories.keys():
        values1 = output_directories[k1]
        if type(values1) == str:
            out_dirs[k1] = './' + options.material_name + '/' + values1
        else:
            out_dirs[k1] = {}
            for k2 in values1.keys():
                values2 = values1[k2]
                out_dirs[k1][k2] = './' + options.material_name + '/' + values2
    
    ### command to run VASP jobs
    command_vasp = {
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
            command=command_vasp)
    
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
    commands = {
            'alamode': {
                'mpirun': options.mpirun, 
                'nprocs': 1, 
                'nthreads': options.ncores, 
                'anphon': options.command_anphon,
                'alm': options.command_alm,
                },
            'vasp': command_vasp,
            }
    
    almcalc = AlamodeCalc(
            apdb.primitive,
            material_name=options.material_name,
            restart=options.restart,
            primitive_matrix=pmat,
            scell_matrix=smat,
            cutoff2=-1, cutoff3=options.cutoff3, 
            mag=0.01,
            nac=phdb.nac,
            commands=commands,
            verbosity=options.verbosity
            )
    
    ### ASE calculator for forces
    mode = 'force'
    calc_force = apdb.get_calculator(mode,
            kpts=phdb.get_kpoints(mode=mode).kpts[0]
            )
    
    t11 = datetime.datetime.now()
    
    ##### suggest and creat structures for harmonic FCs
    ## ver.1: with ALM library
    #almcalc.calc_forces(order=1, calculator=calc_force)
    #
    ## ver.2: with alm command
    almcalc.write_alamode_input(propt='suggest', order=1)
    almcalc.run_alamode(propt='suggest', order=1)
    almcalc.calc_forces(order=1, calculator=calc_force)
    
    t12 = datetime.datetime.now()
    times['harm_forces'] = t12 - t11
    
    ### calculate forces for harmonic FCs
    almcalc.write_alamode_input(propt='fc2')
    almcalc.run_alamode(propt='fc2', neglect_log=True)

    ### calculate band
    almcalc.write_alamode_input(propt='band')
    almcalc.run_alamode(propt='band', neglect_log=True)

    ### calculate DOS
    almcalc.write_alamode_input(propt='dos')
    almcalc.run_alamode(propt='dos', neglect_log=True)
    almcalc.plot_bandos()

    ### Check negative frequency
    if almcalc.frequency_range[0] < options.negative_freq:
        print("")
        print(" Negative eigenvalues were found. Stop the calculation.")
        print(" Minimum frequency : %.2f" % (almcalc.frequency_range[0]))
        print("")
        exit()
    
    t13 = datetime.datetime.now()
    times['harm_alamode'] = t13 - t12
    
    ### calculate forces for cubic FCs
    ## ver.1: with ALM library
    #mode = 'force'
    #almcalc.calc_forces(
    #        2, calc_force, 
    #        nmax_suggest=options.nmax_suggest,
    #        frac_nrandom=options.frac_nrandom,
    #        temperature=options.random_disp_temperature,
    #        )
    ## ver.2: with alm command
    almcalc.write_alamode_input(propt='suggest', order=2)
    almcalc.run_alamode(propt='suggest', order=2)
    almcalc.calc_forces(order=2, calculator=calc_force)
    
    t21 = datetime.datetime.now()
    times['cubic_forces'] = t21 - t13
    
    ### calculate anharmonic force constants
    if almcalc.lasso:
        from ..io.vasp import get_dfset
        for propt in ['cv', 'lasso']:
            almcalc.write_alamode_input(propt=propt)
            almcalc.run_alamode(propt)
    else: 
        ## ver.1 : with ALM library
        ##almcalc.calc_anharm_force_constants()
        ## ver.2: with alm command
        almcalc.write_alamode_input(propt='fc3')
        almcalc.run_alamode(propt='fc3', neglect_log=True)
    
    t22 = datetime.datetime.now()
    times['cubic_fcs'] = t22 - t21
    
    ### calculate kappa
    almcalc.write_alamode_input(propt='kappa', kpts=[15,15,15])
    almcalc.run_alamode(propt='kappa')
    almcalc.plot_kappa()
    
    t23 = datetime.datetime.now()
    times['kappa'] = t23 - t22
    times['total'] = t23 - t11
    

    ### END of calculations
    print_times(times)
        
    end_autokappa()
    
