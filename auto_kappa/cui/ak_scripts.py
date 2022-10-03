#
# ak_script.py
#
# This file is akrun command user interface.
#
# Copyright (c) 2022 Masato Ohnishi
#
# This file is distributed under the terms of the MIT license.
# Please see the file 'LICENCE.txt' in the root directory
# or http://opensource.org/licenses/mit-license.php for information.
#
import numpy as np
import datetime

from auto_kappa.io.phonondb import Phonondb
from auto_kappa.apdb import ApdbVasp
from auto_kappa.alamode.almcalc import AlamodeCalc
from auto_kappa import output_directories
from auto_kappa.cui.ak_parser import get_parser
from auto_kappa.alamode.log_parser import AkLog

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
    
    ###############################
    if options.neglect_log == 1:
        neglect_log = True
    elif options.neglect_log == 0:
        neglect_log = False
    ###############################

    ### calculate forces for harmonic FCs
    almcalc.write_alamode_input(propt='fc2')
    almcalc.run_alamode(propt='fc2', neglect_log=neglect_log)
    
    ### calculate band
    almcalc.write_alamode_input(propt='band')
    almcalc.run_alamode(propt='band', neglect_log=neglect_log)

    ### calculate DOS
    almcalc.write_alamode_input(propt='dos')
    almcalc.run_alamode(propt='dos', neglect_log=neglect_log)
    almcalc.plot_bandos()
    
    ### eigenvalues at commensurate points
    almcalc.write_alamode_input(propt='evec_commensurate')
    almcalc.run_alamode(propt='evec_commensurate', neglect_log=neglect_log)
    
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
    times['anharm_forces'] = t21 - t13
    
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
        almcalc.run_alamode(propt='fc3', neglect_log=neglect_log)
    
    t22 = datetime.datetime.now()
    times['anharm_fcs'] = t22 - t21
    
    ### calculate kappa
    almcalc.write_alamode_input(propt='kappa', kpts=[15,15,15])
    almcalc.run_alamode(propt='kappa')
    
    ### analyze phonons
    print()
    print()
    print(" Plot anharmonic properties:")
    print(" ---------------------------")
    almcalc.plot_kappa()
    almcalc.plot_lifetime(temperatures="300:500")
    almcalc.plot_scattering_rates(temperature=300., grain_size=1000.)
    almcalc.plot_cumulative_kappa(
            temperatures="300:500", wrt='frequency', xscale='linear')
    almcalc.plot_cumulative_kappa(
            temperatures="300:500", wrt='mfp', xscale='log')
     
    ### get time
    t23 = datetime.datetime.now()
    times['kappa'] = t23 - t22
    times['total'] = t23 - t11
    
    ### output log.yaml and fig_times.png
    log = AkLog(options.material_name)
    log.write_yaml()
    log.plot_times()
    
    ### END of calculations
    print_times(times)
    
    end_autokappa()
    
