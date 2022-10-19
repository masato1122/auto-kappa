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
import os
import numpy as np
import datetime

from auto_kappa.io.phonondb import Phonondb
from auto_kappa.apdb import ApdbVasp
from auto_kappa.alamode.almcalc import AlamodeCalc
from auto_kappa import output_directories
from auto_kappa.cui.ak_parser import get_parser
from auto_kappa.alamode.log_parser import AkLog
from auto_kappa.structure.crystal import get_automatic_kmesh

def start_autokappa():
    """ Print the logo.
    Font: stick letters
    """
    from auto_kappa.version import __version__
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

def print_options(options):

    msg = "\n"
    msg += " Input parameters:\n"
    msg += " =================\n"
    print(msg)
    opt_dict = eval(str(options))
    for key in opt_dict:
        if opt_dict[key] is not None:
            print("", key.ljust(20), " : ", opt_dict[key])
    print("")

def _modify_directory_name(mater_name, propt=None):
    
    dir1 = mater_name + "/cube/%s" % propt
    dir2 = mater_name + "/cube/%s_fd" % propt
    if os.path.exists(dir1) == False:
        return None
    if os.path.exists(dir2) == False:
        os.rename(dir1, dir2)    
        print("")
        print(" Modify directory name: %s => %s" % (dir1, dir2))
        return None
    ##
    for i in range(1, 100):
        dir2_new = mater_name + "/cube/%s_fd-backup%d" % (propt, i)
        if os.path.exists(dir2_new) == False:
            os.rename(dir2, dir2_new)
            print("")
            print(" Modify directory name: %s => %s" % (dir2, dir2_new))
            break
    ##
    os.rename(dir1, dir2)
    print("")
    print(" Modify directory name: %s => %s" % (dir1, dir2))

def main():
    
    times = {}
    
    start_autokappa()
    
    options = get_parser()
    print_options(options)

    ### Read data of phonondb
    ### phonondb is used to obtain structures (primitive, unit, and super cells)
    ### and k-points.
    phdb = Phonondb(options.directory)
    
    ### get required data from Phonondb
    pmat = phdb.primitive_matrix
    smat = phdb.scell_matrix
    unitcell = phdb.get_unitcell(format='ase')
    kpts_for_relax = phdb.get_kpoints(mode='relax').kpts[0]
    kpts_for_force = phdb.get_kpoints(mode='force').kpts[0]
    nac = phdb.nac
    if nac == 1:
        kpts_for_nac = phdb.get_kpoints(mode='nac').kpts[0]
    else:
        kpts_for_nac = None

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
            unitcell,
            primitive_matrix=pmat,
            scell_matrix=smat,
            command=command_vasp,
            )
    
    ### Relaxation calculation
    mode = 'relax'
    apdb.run_vasp(mode, 
            out_dirs[mode],
            kpts_for_relax,
            print_params=True
            )
    
    ### Born effective charge
    if nac == 1:
        mode = 'nac'
        apdb.run_vasp(mode,
                out_dirs[mode],
                kpts_for_nac,
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
            magnitude=0.01,
            magnitude2=0.03,
            nac=nac,
            commands=commands,
            verbosity=options.verbosity
            )
    
    ### ASE calculator for forces
    mode = 'force'
    calc_force = apdb.get_calculator(mode,
            kpts=kpts_for_force,
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
        
        log = AkLog(options.material_name)
        log.write_yaml()
        log.plot_times()
        
        print("")
        print(" Negative eigenvalues were found. Stop the calculation.")
        print(" Minimum frequency : %.2f" % (almcalc.frequency_range[0]))
        print("")
        exit()
    
    t13 = datetime.datetime.now()
    times['harm_alamode'] = t13 - t12
    
    ##############################
    ##                          ##
    ##  Start FC3 calculation   ##
    ##                          ##
    ##############################
    
    ## modify directory name automatically which has been created with a old
    ## version
    _modify_directory_name(almcalc.material_name, propt='force')
    _modify_directory_name(almcalc.material_name, propt='kappa')
    
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
    almcalc.calc_forces(order=2, calculator=calc_force,
            nmax_suggest=options.nmax_suggest,
            frac_nrandom=options.frac_nrandom,
            output_dfset=2,
            )
    
    t21 = datetime.datetime.now()
    times['anharm_forces'] = t21 - t13
    
    ### calculate anharmonic force constants
    #if almcalc.lasso:
    if almcalc.fc3_type == 'lasso':
        from auto_kappa.io.vasp import get_dfset
        for propt in ['cv', 'lasso']:
            order = 2
            almcalc.write_alamode_input(propt=propt, order=order)
            almcalc.run_alamode(propt, order=order, neglect_log=neglect_log)
    else: 
        ## ver.1 : with ALM library
        ##almcalc.calc_anharm_force_constants()
        ## ver.2: with alm command
        almcalc.write_alamode_input(propt='fc3')
        almcalc.run_alamode(propt='fc3', neglect_log=neglect_log)
    
    t22 = datetime.datetime.now()
    times['anharm_fcs'] = t22 - t21
    
    ### calculate kappa
    kpts = get_automatic_kmesh(almcalc.primitive)
    almcalc.write_alamode_input(propt='kappa', order=2, kpts=kpts)
    almcalc.run_alamode(propt='kappa', neglect_log=neglect_log)
    
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
    
