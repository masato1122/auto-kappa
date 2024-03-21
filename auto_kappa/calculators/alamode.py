#
# alamode.py
#
# This script helps to conduct phonon analyses using ALAMODE and VASP.
#
# Copyright (c) 2022 Masato Ohnishi
#
# This file is distributed under the terms of the MIT license.
# Please see the file 'LICENCE.txt' in the root directory
# or http://opensource.org/licenses/mit-license.php for information.
#
import sys
import os
import os.path
import numpy as np
import glob
import shutil

from auto_kappa.apdb import ApdbVasp
from auto_kappa import output_directories
from auto_kappa.alamode.almcalc import AlamodeCalc, should_rerun_alamode
from auto_kappa.alamode.log_parser import AkLog
from auto_kappa.structure.crystal import get_automatic_kmesh
from auto_kappa.cui.suggest import klength2mesh
from auto_kappa.cui import ak_log

import logging
logger = logging.getLogger(__name__)

def analyze_phonon_properties(
        almcalc, calc_force=None, negative_freq=-1e-3, 
        material_name=None, neglect_log=False, 
        #
        harmonic_only=False, calc_kappa=True,
        #
        nmax_suggest=None, frac_nrandom=1.0,
        #
        params_nac={'apdb': None, 'kpts': None},
        kdensities_for_kappa=[500, 1000, 1500],
        ##
        scph=0, disp_temp=500., frac_nrandom_higher=0.34,
        ):
    """ Analyze phonon properties
    
    Args
    =====
    
    almcalc : auto_kappa.alamode.almcalc.AlamodeCalc
    
    calc_force : ase.calculators.vasp.vasp.Vasp

    negative_freq : float   
        threshold for negative frequency
    
    harmonic_only : int, 0
        0 : stop the calculation after harmonic properties are calculated.
    
    calc_kappa : bool
        Calculate thermal conductivity (True) or not (False)

    nmax_suggest : int
        Maximum number of suggested partterns for finite displacement method. If
        the number exceeds this value, LASSO is used for cubic FCs.

    frac_nrandom : int
        Coefficient to determine the number of random patterns for LASSO.

    params_nac : dict
        parameters for NAC calculation.
        ``params_nac[apdb]``: auto_kappa.apdb.ApdbVasp
    
    scph : int
        Self-consistent phonon approach is used (1) or not (0).
    
    disp_temp : float, 500
        temperature for random displacement
    
    frac_nrandom_higher : float, 0.34
        fractional number of random displacement patterns for high-order FCs

    Return
    =======
    integer :
        ``0`` when the calculation was finished properly,
        ``-1`` when negative frequencies were found,
    
    """
    ###
    ### Calculate harmonic FCs and harmonic phonon properties
    ###
    analyze_harmonic_properties(
            almcalc, calc_force, negative_freq=negative_freq, 
            params_nac=params_nac)
    
    ### Check negative frequency
    if almcalc.minimum_frequency < negative_freq:
        
        ### If negative frequencies were found,
        log = AkLog(material_name)
        log.write_yaml()
        ak_log.negative_frequency(almcalc.minimum_frequency)
        
        ### If SCPH is not used and negative frequencies were found, return -1.
        if scph == 0:
            return -1
    
    if harmonic_only:
        msg = "\n"
        msg += " Harmonic properties have been calculated.\n"
        logger.info(msg)
        return 0
    
    ###
    ### Calculate cubic FCs
    ###
    calculate_cubic_force_constants(
            almcalc, calc_force,
            nmax_suggest=nmax_suggest, 
            frac_nrandom=frac_nrandom, 
            neglect_log=neglect_log)
    
    ###
    ### Calculate high-order FCs using LASSO
    ###
    if scph:
        calculate_high_order_force_constants(
                almcalc, calc_force,
                frac_nrandom=frac_nrandom_higher,
                disp_temp=disp_temp,
                )

        ### Perform SCP calculation
        print(" Conduct SCP calculation!")
        print(" ^^^^^^^^^^^^^^^^^^^^^^^^^")
        sys.exit()
        
    ###
    ### Calculate kappa with different k-mesh densities
    ###
    if calc_kappa:
        if scph == 0:
            calculate_thermal_conductivities(almcalc, 
                kdensities=kdensities_for_kappa,
                neglect_log=neglect_log,
                temperatures_for_spectral="300:500"
                )
        else:
            #calculate_thermal_conductivities_scph(almcalc,
            #        kdensities=kdensities_for_kappa[1],
            #        temperatures=[300])
            pass
    
    ### output log.yaml 
    log = AkLog(material_name)
    log.write_yaml()
    
    return 0

def analyze_harmonic_with_larger_supercells(
        almcalc_orig, base_dir=None,
        #
        max_natoms_init=None,
        delta_max_natoms=50,
        max_loop=1,
        #
        k_length=20,
        negative_freq=-1e-3,
        neglect_log=False,
        #
        restart=1,
        ):
    """ Analyze harmonic properties with larger supercell(s) """
    from auto_kappa.structure.supercell import estimate_supercell_matrix
    from auto_kappa.cui.ak_log import start_larger_supercell

    count = 0
    isc = 0
    sc_mat_prev = almcalc_orig.scell_matrix
    almcalc_new = None
    while isc <= max_loop:
        
        count += 1
        max_natoms = max_natoms_init + delta_max_natoms * count
        
        ### estimate a supercell from the max number of atoms
        sc_mat = estimate_supercell_matrix(
                almcalc_orig.unitcell, max_num_atoms=max_natoms
                )
        
        ### check if the supercell is changed
        diff = np.amax(np.asarray(sc_mat) - np.asarray(sc_mat_prev))
        
        if diff < 0.1:
            continue
        
        ### label
        sc_label = "%dx%dx%d" % (sc_mat[0,0], sc_mat[1,1], sc_mat[2,2])
        added_dir = "sc-" + sc_label
        
        ### Set AlamodeCalc obj
        almcalc_new = AlamodeCalc(
                almcalc_orig.primitive,
                base_directory=almcalc_orig.base_directory,
                additional_directory=added_dir,
                restart=restart,
                primitive_matrix=almcalc_orig.primitive_matrix,
                scell_matrix=sc_mat,
                cutoff2=almcalc_orig.cutoff2,
                cutoff3=almcalc_orig.cutoff3,
                magnitude=almcalc_orig.magnitude,
                magnitude2=almcalc_orig.magnitude2,
                nac=almcalc_orig.nac,
                commands=almcalc_orig.commands,
                verbosity=almcalc_orig.verbosity,
                yamlfile_for_outdir=almcalc_orig.yamlfile_for_outdir,
                )
        
        start_larger_supercell(almcalc_new)
        
        ### Set kmesh for VASP calculation
        kpts = klength2mesh(k_length, almcalc_new.supercell.cell.array)
        
        ### ASE calculator for forces
        apdb = ApdbVasp(
                almcalc_orig.unitcell,
                primitive_matrix=almcalc_orig.primitive_matrix,
                scell_matrix=sc_mat,
                command=almcalc_orig.commands["vasp"],
                )
        calc_force = apdb.get_calculator('force', kpts=kpts)
        
        ### analyze phonon properties with a supercell
        analyze_harmonic_properties(
                almcalc_new, calc_force,
                negative_freq=negative_freq)
        
        ### Finish
        if almcalc_new.minimum_frequency < negative_freq:
            ak_log.negative_frequency(almcalc_new.minimum_frequency)
            return almcalc_new
        
        isc += 1
        if isc == max_loop:
            break
        
        ### prepare for the next step
        sc_mat_prev = almcalc_new.scell_matrix
    
    return almcalc_new

def analyze_harmonic_properties(
        almcalc, calculator, 
        neglect_log=False, max_num_corrections=5,
        deltak=0.01, reciprocal_density=1500,
        negative_freq=-1e-3,
        params_nac={'apdb': None, 'kpts': None}
        ):
    """ Analyze harmonic FCs and phonon properties.

    Args
    =====
    almcalc : AlamodeCalc obj

    calculator : ASE calculator for VASP
    
    deltak : float
        resolution of phonon dispersion

    reciprocal_density : float
        resolution of DOS
    
    """
    from auto_kappa.alamode.log_parser import get_minimum_frequency_from_logfile
    
    #logger = logging.getLogger(__name__)
    
    ##### suggest and creat structures for harmonic FCs
    almcalc.write_alamode_input(propt='suggest', order=1)
    almcalc.run_alamode(propt='suggest', order=1)
    
    ### calculate forces
    almcalc.calc_forces(order=1, calculator=calculator)
    
    _ncores_orig = almcalc.commands['alamode']['ncores']
    _para_orig = almcalc.commands['alamode']['anphon_para']
    
    nac_orig = almcalc.nac
    for propt in ["fc2", "band", "dos"]:
        
        ### ver.1
        #_analyze_each_harm_propt(
        #        propt, almcalc, 
        #        neglect_log=neglect_log,
        #        max_num_corrections=max_num_corrections,
        #        deltak=deltak, reciprocal_density=reciprocal_density)
        
        ### ver.2: modified
        almcalc.analyze_harmonic_property(
                propt, 
                max_num_corrections=max_num_corrections,
                deltak=deltak, 
                reciprocal_density=reciprocal_density
                )
        
    ### optimize NAC option
    fmin = almcalc.get_minimum_frequency(which="both")
    
    ### calculate Born effective charge
    if fmin < negative_freq and almcalc.nac == 0:
        
        try:
            ## set NONANALYTICAL option for Alamode
            almcalc.nac = 2
            mode = 'nac'
            params_nac["apdb"].run_vasp(
                    mode,
                    almcalc.out_dirs[mode],
                    params_nac["kpts"],
                    print_params=True
                    )
        except Exception:
            msg = "\n Warning: cannot get parameters for NAC."
            logger.warning(msg)
         
    ### get optimal NONANALYTICAL option for Alamode
    if fmin < negative_freq and almcalc.nac != 0:
        
        nac_new = almcalc.get_optimal_nac(
                tol_neg_frac=0.03,
                max_num_corrections=max_num_corrections,
                deltak=deltak,
                reciprocal_density=reciprocal_density,
                negative_freq=negative_freq,
                )
        
        ### If a different NAC can remove negative frequencies,
        ### harmonic properties are calculated again.
        if nac_new is not None:
            
            ### message
            nac_orig = almcalc.nac
            msg = "\n Modify NONANALYTIC option from %d to %d, " % (
                    nac_orig, nac_new)
            msg += "which can remove negative frequencies."
            logger.info(msg)
            
            ### update NAC
            almcalc.nac = nac_new
            
            ### copy the calculated results with the optimal NAC option
            try:
                dir_bandos = almcalc.out_dirs["harm"]["bandos"]
                dir_nac = dir_bandos + "/nac_%d" % almcalc.nac
                msg = "\n >>> Move files in %s to %s." % (
                        almcalc.get_relative_path(dir_nac),
                        almcalc.get_relative_path(dir_bandos))
                logger.info(msg)
                if os.path.exists(dir_nac):
                    fns = glob.glob(dir_nac + "/*")
                    for ff in fns:
                        shutil.copy(ff, dir_bandos)
                neg_log = 0
            except Exception:
                neg_log = 1
            
            ### overwrite harmonic properties
            for propt in ["band", "dos"]:
                almcalc.analyze_harmonic_property(
                        propt,
                        max_num_corrections=max_num_corrections,
                        deltak=deltak, 
                        reciprocal_density=reciprocal_density,
                        neglect_log=neg_log,
                        )
    
    ###
    almcalc.commands['alamode']['ncores'] = _ncores_orig
    almcalc.commands['alamode']['anphon_para'] = _para_orig
    
    ### plot band and DOS
    almcalc.plot_bandos()
    
    ### eigenvalues at commensurate points
    almcalc.write_alamode_input(propt='evec_commensurate')
    almcalc.run_alamode(propt='evec_commensurate', neglect_log=neglect_log)
    
def calculate_cubic_force_constants(
        almcalc, calculator,
        nmax_suggest=None, frac_nrandom=None, neglect_log=False
        ):
    """ Calculate cubic force constants
    """
    ### calculate forces for cubic FCs
    almcalc.write_alamode_input(propt='suggest', order=2)
    almcalc.run_alamode(propt='suggest', order=2)
    almcalc.calc_forces(
            order=2, calculator=calculator,
            nmax_suggest=nmax_suggest,
            frac_nrandom=frac_nrandom,
            output_dfset=2,
            )
    
    ### calculate anharmonic force constants
    if almcalc.fc3_type == 'lasso':
        for propt in ['cv', 'lasso']:
            almcalc.write_alamode_input(propt=propt, order=2)
            almcalc.run_alamode(propt, order=2, neglect_log=neglect_log)
    else: 
        ## ver.1 : with ALM library
        ##almcalc.calc_anharm_force_constants()
        ## ver.2: with alm command
        almcalc.write_alamode_input(propt='fc3')
        almcalc.run_alamode(propt='fc3', neglect_log=neglect_log)

def calculate_high_order_force_constants(
        almcalc, calculator, order=5, frac_nrandom=None, 
        disp_temp=500,):
    """ Calculate high-order (up to 6th-order) force constants finally to
    calculate 4th order FCs for phonon renormalization.
    """
    ### calculate forces for high-order FCs
    almcalc.write_alamode_input(propt='suggest', order=order)
    almcalc.run_alamode(propt='suggest', order=order)
    almcalc.calc_forces(
            order=order, calculator=calculator,
            frac_nrandom=frac_nrandom,
            temperature=disp_temp,
            output_dfset=2,
            )
    
    ### calculate anharmonic force constants
    for propt in ['cv', 'lasso']:
        almcalc.write_alamode_input(propt=propt, order=order)
        almcalc.run_alamode(propt, order=order, neglect_log=False)


#def calculate_thermal_conductivities_scph(
#        almcalc, kdensities=[1000], temperatures=[300],
#        **kwargs,):
#    """ Calculate thermal conductivities with SCP theory """
#    pass

    
def calculate_thermal_conductivities(
        almcalc, 
        kdensities=[500, 1000, 1500], 
        neglect_log=False,
        temperatures_for_spectral="300:500",
        **kwargs,
        ):
    """ Calculate thermal conductivities and plot spectral info
    
    Args
    ======
    almcalc : AlamodeCalc obj

    kdensities : array of float

    neglect_log : bool

    temperatures_for_spectral : string of floats seperated by ":"

    """
    for kdensity in kdensities:
        
        kpts = get_automatic_kmesh(
                almcalc.primitive, reciprocal_density=kdensity)
            
        outdir = (
                almcalc.out_dirs['cube']['kappa_%s' % almcalc.fc3_type] + 
                "_%dx%dx%d" % (int(kpts[0]), int(kpts[1]), int(kpts[2]))
                )
        
        ###
        almcalc.write_alamode_input(
                propt='kappa', order=2, kpts=kpts, outdir=outdir, **kwargs)
        
        almcalc.run_alamode(
                propt='kappa', neglect_log=neglect_log, outdir=outdir)
        
        ### check output file for kappa
        kappa_log = outdir + "/kappa.log"
        flag = should_rerun_alamode(kappa_log) 
        
        if flag and almcalc.commands['alamode']['anphon_para'] == "mpi":
            ##
            ## if thermal conductivity was not calculated with MPI, run the 
            ## calculation again with OpenMP.
            ##
            ak_log.rerun_with_omp()
            almcalc.commands['alamode']['anphon_para'] = "omp"
            almcalc.run_alamode(
                    propt='kappa', neglect_log=neglect_log, 
                    outdir=outdir, **kwargs)
        
    ### analyze phonons
    msg = "\n"
    msg += " Plot anharmonic properties:\n"
    msg += " ---------------------------"
    logger.info(msg)
    
    out = almcalc.plot_kappa()
    
    #if out < 0:
    #    msg = "\n"
    #    msg += " Warning: thermal conductivity may be too small."
    #    logger.error(msg)
    
    ###
    try:
        almcalc.plot_lifetime(temperatures=temperatures_for_spectral)
    except Exception:
        logger.warning("\n Warning: "\
                "the figure of lifetime was not created properly.")
    
    try:
        almcalc.plot_scattering_rates(temperature=300., grain_size=1000.)
    except Exception:
        logger.warning("\n Warning: the figure of "\
                "scattering rate was not created properly.")
    
    try:
        almcalc.plot_cumulative_kappa(
                temperatures=temperatures_for_spectral, 
                wrt='frequency', xscale='linear')
    except Exception:
        logger.warning("\n Warning: the figure of "\
                "cumulative thermal conductivity was not created properly.")
    
    try:
        almcalc.plot_cumulative_kappa(
                temperatures=temperatures_for_spectral, 
                wrt='mfp', xscale='log')
    except Exception:
        logger.warning("\n Warning: the figure of "\
                "cummulative TCs was not created properly.")
    

