#
# ak_script.py
#
# This script is akrun command user interface.
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
import datetime
import time

from auto_kappa.apdb import ApdbVasp
from auto_kappa import output_directories
from auto_kappa.alamode.almcalc import AlamodeCalc
from auto_kappa.alamode.pes import calculate_pes
from auto_kappa.io.files import write_output_yaml
from auto_kappa.calculators.alamode import (
        analyze_phonon_properties,
        analyze_harmonic_with_larger_supercells,
        calculate_cubic_force_constants,
        calculate_thermal_conductivities,
        )
from auto_kappa.cui.ak_parser import get_parser
from auto_kappa.cui import ak_log
from auto_kappa.cui.initialization import (
        use_omp_for_anphon,
        get_previous_nac,
        get_required_parameters,
        get_base_directory_name,
        read_phonondb,
        )

import logging
logging.getLogger('matplotlib').setLevel(logging.WARNING)
logging.getLogger('custodian').setLevel(logging.WARNING)
logger = logging.getLogger(__name__)

import warnings
warnings.filterwarnings('ignore')

def _set_outdirs(base_dir):
    """ Prepare and return names of output directories """
    out_dirs = {}
    for k1 in output_directories.keys():
        values1 = output_directories[k1]
        if type(values1) == str:
            out_dirs[k1] = base_dir + '/' + values1
        else:
            out_dirs[k1] = {}
            for k2 in values1.keys():
                values2 = values1[k2]
                out_dirs[k1][k2] = base_dir + '/' + values2
    return out_dirs

def _stop_symmetry_error(out):
    """ Stop the calculation because of the symmetry change """
    if out == -1:
        msg = "\n Error: crystal symmetry was changed during the "\
                "relaxation calculation."
    elif out == -2:
        msg = "\n Error: too many errors for the relaxation calculation."
    else:
        msg = ""
    
    msg += "\n"
    msg += "\n STOP THE CALCULATION"
    time = datetime.datetime.now()
    msg += "\n at " + time.strftime("%m/%d/%Y %H:%M:%S")
    msg += "\n"
    logger.error(msg)
    sys.exit()
    
def main():
    
    ### Parse given parameters
    options = get_parser()
    ak_params = eval(str(options))
    
    ### Get the name of the base directory
    base_dir = get_base_directory_name(
            ak_params['material_name'], restart=ak_params['restart']
            )
    os.makedirs(base_dir, exist_ok=True)
    
    ### set logger
    logfile = base_dir + "/ak.log"
    ak_log.set_logging(filename=logfile, level=logging.DEBUG, format="%(message)s")
    
    ### Start auto-kappa
    ak_log.start_autokappa()
    
    #ak_log.print_machine_info()
    
    ### Set output directories
    out_dirs = _set_outdirs(base_dir)
    
    ### memory check
    if use_omp_for_anphon(base_dir):
        if ak_params["anphon_para"] != "omp":
            msg = "\n Change anphon_para option to \"omp\".\n"
            logger.info(msg)
            ak_params["anphon_para"] = "omp"
    
    ######################
    ### Adjust options ###
    ######################
    ### max_natoms3: maximun limit of the number of atoms for FC3
    #if options.max_natoms3 is None:
    #    options.max_natoms3 = options.max_natoms
    
    ### relaxed_cell
    if ak_params['relaxed_cell'] is not None:
        if (ak_params['relaxed_cell'].lower()[0] == "u" or 
                ak_params['relaxed_cell'].lower()[0] == "c"):
            ak_params['relaxed_cell'] = "unitcell"
        elif ak_params['relaxed_cell'].lower()[0] == "p":
            ak_params['relaxed_cell'] = "primitive"
     
    ### Get required parameters for the calculation!
    cell_types, structures, trans_matrices, kpts_used, nac = (
            get_required_parameters(
                base_directory=base_dir,
                dir_phdb=ak_params['directory'], 
                file_structure=ak_params['file_structure'],
                max_natoms=ak_params['max_natoms'], 
                #max_natoms3=ak_params['max_natoms3'],
                k_length=ak_params['k_length'],
                celltype_relax_given=ak_params['relaxed_cell'],
                )
            )
    
    ### NONANALYTIC (primitive)
    if nac != 0:
        if ak_params['nonanalytic'] is not None:
            nac = ak_params['nonanalytic']
        
        ### check previously-used NONANALYTIC parameter
        try:
            dir0 = base_dir.replace(os.getcwd(), ".")
        except Exception:
            dir0 = base_dir
        
        prev_nac = get_previous_nac(dir0)
        if prev_nac is not None and nac != prev_nac:
            msg = "\n NONANALYTIC was modified to %s" % prev_nac
            logger.info(msg)
            nac = ak_params["nonanalytic"] = prev_nac
    
    ### print parameters
    ak_log.print_options(ak_params)
    
    ak_log.print_conditions(
            cell_types=cell_types, 
            trans_matrices=trans_matrices,
            kpts_all=kpts_used,
            )
    
    ### write file
    os.makedirs(out_dirs["result"], exist_ok=True)
    filename = out_dirs["result"] + "/parameters.yaml"
    ak_log.write_parameters(
            filename,
            structures["unitcell"], 
            cell_types, trans_matrices, kpts_used, nac
            )
    
    try:
        fn_print = filename.replace(os.getcwd(), ".")
    except Exception:
        fn_print = filename
    msg = "\n Output %s" % fn_print
    logger.info(msg)
    
    ### output yaml file
    yaml_outdir = base_dir + "/output_directories.yaml"
    info = {"directory": out_dirs["result"].replace(base_dir, "."),
            "kind": "others",
            "note": "results"}
    write_output_yaml(yaml_outdir, "result", info, overwrite=False)
    
    ### For materials with negative frequencies, options.max_natoms3 may be
    ### different from options.max_natoms.
    #if options.max_natoms != options.max_natoms3:
    #    msg = " Error: max_natoms != max_natoms3 is not supported yet."
    #    logger.error(msg)
    #    exti()
    
    ### command to run VASP jobs
    command_vasp = {
            'mpirun': ak_params['mpirun'], 
            'nprocs': ak_params['ncores'], 
            'nthreads': 1, 
            'vasp': ak_params['command_vasp'],
            }
    
    ### Set ApdbVasp object
    apdb = ApdbVasp(
            structures["unitcell"],
            primitive_matrix=trans_matrices["primitive"],
            scell_matrix=trans_matrices["supercell"],
            command=command_vasp,
            #yamlfile_for_outdir=yaml_outdir
            )
    
    ### Relaxation calculation
    out = apdb.run_relaxation(
            out_dirs["relax"],
            kpts_used["relax"],
            volume_relaxation=ak_params['volume_relaxation'],
            cell_type=cell_types["relax"],
            max_error=ak_params["max_relax_error"],
            nsw_params=ak_params["nsw_params"],
            )
    
    ### Stop the calculation because of the symmetry error
    if out < 0:
        _stop_symmetry_error(out)
    
    ### output yaml file
    info = {
            "directory": out_dirs["relax"].replace(base_dir, "."), 
            "kind": "others",
            "note": "structure optimization",
            }
    write_output_yaml(yaml_outdir, "relax", info)
    
    ##############################
    ### Get relaxed structures
    ### This part was omitted in the beggining of ver.0.2.
    structures_relax = apdb.structures.copy()
    
    ### Calculate Born effective charge
    if nac:
        mode = 'nac'
        apdb.run_vasp(mode, out_dirs[mode], kpts_used["nac"], print_params=True)
        
        ### output yaml file
        info = {
                "directory": out_dirs[mode].replace(base_dir, "."), 
                "kind": "VASP",
                "note": "Born effective charge",
                }
        write_output_yaml(yaml_outdir, mode, info)
    
    ### command for ALAMODE
    command_alamode = {
        'mpirun': ak_params['mpirun'], 
        'anphon_para': ak_params['anphon_para'], 
        'ncores': ak_params['ncores'], 
        'anphon': ak_params['command_anphon'],
        'alm': ak_params['command_alm'],
        }

    ### Set AlmCalc
    almcalc = AlamodeCalc(
            structures_relax['prim'],
            base_directory=base_dir,
            restart=ak_params['restart'],
            primitive_matrix=trans_matrices['primitive'],
            scell_matrix=trans_matrices['supercell'],
            cutoff2=-1,
            cutoff3=ak_params['cutoff_cubic'],
            magnitude=ak_params['mag_harm'],
            magnitude2=ak_params['mag_cubic'],
            ##mag_high=ak_params['mag_high'],
            nac=nac,
            commands={'alamode': command_alamode, 'vasp': command_vasp},
            verbosity=ak_params['verbosity'],
            yamlfile_for_outdir=yaml_outdir
            )
    
    ### Prepare an ase.calculators.vasp.vasp.Vasp obj for force calculation
    calc_force = apdb.get_calculator('force', kpts=kpts_used['harm'])
    
    ### Analyze phonon properties
    neglect_log = ak_params['neglect_log']
    out = analyze_phonon_properties(
            almcalc,
            calc_force=calc_force,
            negative_freq=ak_params['negative_freq'],
            material_name=ak_params['material_name'],
            neglect_log=neglect_log,
            harmonic_only=ak_params['harmonic_only'],
            #
            nmax_suggest=ak_params['nmax_suggest'],
            frac_nrandom=ak_params['frac_nrandom'],
            #
            params_nac={'apdb': apdb, 'kpts': kpts_used['nac']},
            #
            scph=ak_params['scph'],
            disp_temp=ak_params['random_disp_temperature'],
            frac_nrandom_higher=ak_params['frac_nrandom_higher'],
            )
    
    ### Calculate PES
    if (almcalc.minimum_frequency < ak_params['negative_freq'] and 
            ak_params['pes'] == 2):
        almcalc.calculate_pes(negative_freq=ak_params['negative_freq'])
    
    ########################
    ##  Larger supercell  ##
    ########################
    ### calculate harmonic properteis with larger supercells
    if (almcalc.minimum_frequency < ak_params['negative_freq'] and 
            ak_params["analyze_with_largersc"] == 1):
        
        almcalc_large = analyze_harmonic_with_larger_supercells(
                almcalc,
                base_dir=base_dir,
                #
                max_natoms_init=ak_params['max_natoms'],
                delta_max_natoms=options.delta_max_natoms,
                max_loop=options.max_loop_for_largesc,
                #
                k_length=ak_params['k_length'],
                negative_freq=ak_params['negative_freq'],
                neglect_log=neglect_log,
                #
                restart=ak_params['restart'],
                )

        ### plot band and DOS
        from auto_kappa.plot.bandos import plot_bandos_for_different_sizes
        figname = almcalc.out_dirs["result"] + "/fig_bandos.png"
        plot_bandos_for_different_sizes(
                almcalc_large, almcalc, figname=figname)
        
        ### If negative frequencies could be removed, calculate cubic FCs with
        ### the supercell of initial size
        if (almcalc_large.minimum_frequency > ak_params["negative_freq"] and 
                ak_params['harmonic_only'] == 0):
            
            calculate_cubic_force_constants(
                    almcalc, calc_force,
                    nmax_suggest=ak_params['nmax_suggest'], 
                    frac_nrandom=ak_params['frac_nrandom'], 
                    neglect_log=neglect_log
                    )
            
            ### calculate kappa
            kdensities = [500, 1000, 1500]
            calculate_thermal_conductivities(
                    almcalc_large, 
                    kdensities=kdensities,
                    neglect_log=neglect_log,
                    temperatures_for_spectral="300:500",
                    fc2xml = almcalc_large.fc2xml,
                    fcsxml = almcalc.fc3xml
                    )
        
        else:
            
            ak_log.negative_frequency(almcalc_large.minimum_frequency)
            
            ### calculate PES
            if ak_params['pes'] > 0:
                almcalc_large.calculate_pes(
                        negative_freq=ak_params['negative_freq'])
                sys.exit()
    
    ### plot and print calculation durations
    from auto_kappa.io.times import get_times
    times, labels = get_times(base_dir)
    
    from auto_kappa.plot.pltalm import plot_times_with_pie
    figname = base_dir + "/result/fig_times.png"
    plot_times_with_pie(
            times, labels, 
            figname=almcalc.get_relative_path(figname))
    
    ak_log.print_times(times, labels)
        
    ### END of calculations
    ak_log.end_autokappa()

