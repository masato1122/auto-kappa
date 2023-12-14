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
import sys
import os
import os.path
import numpy as np
import datetime
import yaml
import time
import glob
import shutil

from auto_kappa.io.phonondb import Phonondb
from auto_kappa.apdb import ApdbVasp
from auto_kappa import output_directories
from auto_kappa.cui.ak_parser import get_parser
from auto_kappa.alamode.almcalc import AlamodeCalc, should_rerun_alamode
from auto_kappa.alamode.log_parser import AkLog
from auto_kappa.alamode.pes import calculate_pes
from auto_kappa.structure.crystal import get_automatic_kmesh
from auto_kappa.cui.suggest import klength2mesh
from auto_kappa.io.files import write_output_yaml

from auto_kappa.cui import ak_log

import logging
logging.getLogger('matplotlib').setLevel(logging.WARNING)
logging.getLogger('custodian').setLevel(logging.WARNING)
logger = logging.getLogger(__name__)

import warnings
warnings.filterwarnings('ignore')

def start_autokappa():
    """ Print the logo.
    Font: stick letters
    """
    from auto_kappa.version import __version__
    msg = "\n"
    msg += "              ___  __                __   __       \n"
    msg += "     /\  |  |  |  /  \ __ |__/  /\  |__) |__)  /\  \n"
    msg += "    /~~\ \__/  |  \__/    |  \ /~~\ |    |    /~~\ \n"
    msg += "    ver. %s\n" % __version__
    msg += "\n"
    time = datetime.datetime.now()
    msg += " Start at " + time.strftime("%m/%d/%Y %H:%M:%S") + "\n"
    
    #logger = logging.getLogger(__name__)
    logger.info(msg)

def print_times(times, labels):
    
    #logger = logging.getLogger(__name__)

    msg = "\n"
    msg += " Calculation times:\n"
    msg += " ==================\n"
    msg += "\n"
    logger.info(msg)
    
    nchar = 0
    ttot = np.sum(np.asarray(times))
    for i, lab in enumerate(labels):
        tt = float(times[i])         ### sec
        percentage = 100.*tt/ttot
        line = "%25s (sec): %13.2f (%.1f%%)" % (lab, tt, percentage)
        nchar = max(nchar, len(line))
        
        msg = " " + line
        logger.info(msg)
    
    msg = " " + "-" * (nchar + 5) + "\n"
    msg += " %25s (sec): %13.2f \n" % ("total", ttot)
    logger.info(msg)

def end_autokappa():
    
    time = datetime.datetime.now()
    msg = "\n"
    msg += "     ___       __  \n"
    msg += "    |__  |\ | |  \ \n"
    msg += "    |___ | \| |__/ \n"
    msg += "    \n"
    msg += " at " + time.strftime("%m/%d/%Y %H:%M:%S") + " \n"
    
    #logger = logging.getLogger(__name__)
    logger.info(msg)

def print_options(params):

    msg = "\n"
    msg += " Input parameters:\n"
    msg += " =================\n"
    for key in params.keys():
        if params[key] is not None:
            msg += " " + key.ljust(25) + " : " + str(params[key]) + "\n"
    msg += "\n"
    
    #logger = logging.getLogger(__name__)
    logger.info(msg)

def print_conditions(cell_types=None, trans_matrices=None, kpts_all=None):
    
    msg = ""
    if cell_types is not None:
        msg += " Cell type\n"
        msg += " ---------\n"
        for cal_type in cell_types:
            msg += " %6s : %10s\n" % (cal_type, cell_types[cal_type])
        msg += "\n"
     
    ###
    if trans_matrices is not None:
        msg += " Transformation matrix\n"
        msg += " ---------------------\n"
        for cell_type in trans_matrices:
            msg += " %10s : " % cell_type
            for i in range(3):
                vec = trans_matrices[cell_type][i]
                if cell_type == "primitive":
                    msg += "%.3f " * 3 % tuple(vec)
                else:
                    msg += "%d " * 3 % tuple(vec)
            msg += "\n"
        msg += "\n"
    
    ###
    if kpts_all is not None:
        msg += " k-mesh\n"
        msg += " ------\n"
        for cal_type in kpts_all:
            msg += (" %6s :" % (cal_type) + 
                    " %d" * 3 % tuple(kpts_all[cal_type]) + 
                    "\n")
     
    #logger = logging.getLogger(__name__)
    logger.info(msg)

def _get_celltype4relaxation(ctype_input, base_dir, natoms_prim=None):
    """ Return the cell type (primitive or unit (conventional) cell) used for 
    the relaxation.
    
    Algorithm
    ---------
    1. If the relaxation calculation has been already started, read the 
    previously used cell type.
    
    2. If it is the initial calculation,
    
    2-1. If the --relaxed_cell (``ctype_suggested``) is not given, the cell type
    for the relaxation is "unitcell".

    2-2. If the --relaxed_cell is given, the cell type is deteremined by this
    option.

    Parameters
    -----------
    ctype_input : string
        cell type given with the options (--relaxed_cell)

    base_dir : string
        base output directory for the automation calculation

    natoms_prim : int
        number of atoms in the primitive cell
    """
    import os.path
    import ase.io
    #from auto_kappa.io.vasp import wasfinished
    from auto_kappa.structure.crystal import get_standardized_structure_spglib
    
    #logger = logging.getLogger(__name__)

    cell_type = None
    
    ### check previous calculations
    filenames = [
            base_dir + "/relax/vasprun-xml",
            base_dir + "/relax/freeze-1/vasprun-xml"]
    
    for fn in filenames:
        if os.path.exists(fn):
            atoms = ase.io.read(fn, format='vasp-xml')
            if len(atoms) == natoms_prim:
                cell_type = "primitive"
            elif len(atoms) > natoms_prim:
                cell_type = "unitcell"
            else:
                msg = " Error: cell type cannot be deteremined."
                logger.error(msg)
                sys.exit()
    
    ###
    if cell_type is None:
        if ctype_input is None:

            cell_type = "unitcell"
        
        else:
            
            cell_type = ctype_input
    
    return cell_type

def read_phonondb(directory):
    """ Read data in Phonondb. Phonondb is used to obtain structures 
    (primitive, unit, and super cells) and k-points.
    
    Note that the definition of transformation matrix in Phononpy (Phonondb)
    is different from that in ASE and Pymatgen.
    With Phonopy definition    : P = mat_u2p_pp   @ U
    With ASE or pmg definition : P = mat_u2p_pp.T @ U
    _pp in ``mat_u2p_pp`` and ``mat_u2s_pp`` denotes Phonopy.
    
    To do
    --------
    If the material is not contained in Phonondb, the trasnformation
    matrices need to be obtained with the same manner as spglib. Pymatgen
    may not be able to be used because of different definitions of the
    transformation matrix in Phonopy and Pymatgen.
    
    >>> from auto_kappa.structure.cells import get_mat_u2p_spglib
    >>> unitcell, mat_u2p = get_mat_u2p_spglib(structure)
    """
    phdb = Phonondb(directory)
    
    unitcell = phdb.get_unitcell(format='ase')
    
    if phdb.nac == 1:
        kpts_for_nac = phdb.get_kpoints(mode='nac').kpts[0]
    else:
        kpts_for_nac = None
    
    ##
    matrices = {
            "primitive": phdb.primitive_matrix,
            "supercell": phdb.scell_matrix,
            }

    kpts_all = {
            "relax": phdb.get_kpoints(mode='relax').kpts[0],
            "harm": phdb.get_kpoints(mode='force').kpts[0],
            "nac": kpts_for_nac
            }

    return unitcell, matrices, kpts_all, phdb.nac
            
def _get_base_directory_name(label, restart=True):
    """ Return the base directory name with an absolute path
    """
    outdir = os.getcwd() + "/"
    if restart:
        """ Case 1: when the job can be restarted """
        outdir += label
    else:
        """ Case 2: when the job cannot be restarted """
        if os.path.exists(label) == False:
            outdir += label
        else:
            for i in range(2, 100):
                outdir = label + "-%d" % i
                if os.path.exists(outdir) == False:
                    break
    return outdir

def _get_previously_used_parameters(outdir, lattice_unit, cell_types=None):
    """ Read previously used parameters: transformation matrices and k-meshes.
    
    Parameters
    -----------
    outdir : string
        base output directory
    
    lattice_unit : ndarray, shape=(3,3)
        lattice vectors of unitcell

    Returns
    -------
    params_prev : dict
        keys1 = "trans_matrix" and "kpts"
        keys2 = "relax", "nac", "harm", and "cube"
        "trans_matrix" : transformation matrix
        "kpts" : k-mesh for different calculations
    
    """
    from pymatgen.io.vasp import Poscar, Kpoints
    
    ### prepare dictionaries
    cal_types = ["relax", "nac", "harm", "cube"]
    param_types = ["trans_matrix", "kpts"]
    params_prev = {}
    for k1 in cal_types:
        params_prev[k1] = {}
        for k2 in param_types:
            params_prev[k1][k2] = None
    
    ### kpoints for the relaxation
    for i in range(4):
        
        if i == 0:
            """ for relaxation
            """
            label = "relax"
            dirs = [outdir + "/relax/full-1", outdir + "/relax"]
        
        elif i == 1:
            """ for nac
            """
            label = "nac"
            dirs = [outdir + "/nac"]
        
        elif i == 2:
            """ for harmonic FCs
            """
            label = "harm"
            dirs = [outdir + "/harm/force/prist"]
        
        elif i == 3:
            """ for cubic FCs
            """
            label = "cube"
            dirs = [
                    outdir + "/cube/force_fd/prist",
                    outdir + "/cube/force_lasso/prist",
                    outdir + "/cube/force/prist",
                    ]
        
        ##
        params_prev[label]["kpts"] = None
        params_prev[label]["trans_matrix"] = None
        for dd in dirs:
            fn_kpts = dd + "/KPOINTS"
            fn_structure = dd + "/POSCAR"
            if os.path.exists(fn_kpts) and os.path.exists(fn_structure):
                
                ### k-mesh
                kpts = Kpoints.from_file(fn_kpts).kpts[0]
                params_prev[label]["kpts"] = kpts

                ### structure
                structure = Poscar.from_file(
                        fn_structure, check_for_POTCAR=False).structure
                lattice = structure.lattice.matrix
                
                ### transformation matrix
                ##
                ## cell_trans.T = cell_orig.T @ Mtrans
                ## => Mtrans = (cell_orig.T)^-1 @ cell_trans.T
                ##
                #### modified at April 13, 2023
                mat = np.dot(np.linalg.inv(lattice_unit.T), lattice.T)
                
                ###
                if cell_types[label] in ["primitive"]:
                    params_prev[label]["trans_matrix"] = mat
                else:
                    params_prev[label]["trans_matrix"] = np.rint(mat).astype(int)
                
    return params_prev

def _make_structures(unitcell, 
        primitive_matrix=None, 
        supercell_matrix=None, 
        #supercell_matrix3=None
        ):
    """
    Parameters
    ------------
    unitcell : Structure obj

    Return
    -------
    structures : dict of Structure obj

    """
    from auto_kappa.structure.crystal import change_structure_format
    from phonopy.structure.cells import get_supercell, get_primitive

    unit_pp = change_structure_format(unitcell, format="phonopy")
    
    structures = {}

    structures["unitcell"] = unitcell
    
    if primitive_matrix is not None:
        structures["primitive"] = change_structure_format(
                get_primitive(unit_pp, primitive_matrix), format='ase')
    
    if supercell_matrix is not None:
        structures["supercell"] = change_structure_format(
                get_supercell(unit_pp, supercell_matrix), format='ase')
    
    #if supercell_matrix3 is not None:
    #    structures["supercell3"] = change_structure_format(
    #            get_supercell(unit_pp, supercell_matrix3), format='ase')
    
    return structures

def _determine_kpoints_for_all(
        params_suggest, params_prev,
        prioritize_previous_kpts=True
        ):
    """ Return k-mesh for all the calculation
    
    Algorithm
    ----------
    If ``prioritize_previous_kpts`` is True:

    mat1 = params_suggest[cal_type]["trans_matrix"]
    mat2 = params_prev[cal_type]["trans_matrix"]
    kpts1 = params_suggest[cal_type]["kpts"]
    kpts2 = params_prev[cal_type]["kpts"]
    
    - If ``mat1`` is equal to ``mat2``, use ``kpts2``. Note that ``kpts2`` may
      be None.
    
    - If ``mat1`` is not equal to ``mat2``, use kpts1.

    Parameters
    ----------
    params_*** : dict
        keys = "relax", "nac", "harm", and "cube"
    
    params_***[key] : dict of array
        keys = "trans_matrix" and "kpts"
        transformation matrix and kmesh

    params_suggest :
        parameters suggested or used in Phonondb

    params_prev :
        parameters used in the previous calculation
    
    prioritize_previous_kpts : bool
    
    """
    kpts_all = {"relax": None, "nac": None, "harm": None, "cube": None}
    for cal_type in kpts_all:

        mat1 = params_suggest[cal_type]["trans_matrix"]
        mat2 = params_prev[cal_type]["trans_matrix"]
        
        kpts1 = params_suggest[cal_type]["kpts"]
        kpts2 = params_prev[cal_type]["kpts"]
        
        if mat2 is not None:
            if np.allclose(mat1, mat2):
                if prioritize_previous_kpts and kpts2 is not None:
                    kpts_all[cal_type] = kpts2
        
        if kpts_all[cal_type] is None:
            kpts_all[cal_type] = kpts1
        
        ### check
        if kpts_all[cal_type] is None:
            msg = "\n Error: cannot obtain k-mesh info properly."
            msg += "\n k-mesh for %s is None." % cal_type
            logger.error(msg)
            sys.exit()
    
    return kpts_all

def _get_required_parameters(
        base_directory=None,
        dir_phdb=None, file_structure=None,
        max_natoms=None, 
        k_length=None,
        celltype_relax_given=None,
        ):
    """ Return the required parameters for the automation calculation: 
    structures, transformation matrices, and k-meshes.
    
    Parameters
    -----------
    base_directory : string
        base output directory
    
    dir_phdb : string
        Directory of Phonondb

    file_structure : string
        structure file name. ``dir_phdb`` or ``file_structure`` must be given
        while, if both of them are given, ``dir_phdb`` is used.
    
    max_natoms(3) : int
        maximum limit of the number of atoms in the supercell for FC2(3)
    
    k_length : float
    celltype_relax_given : string

    Return
    --------
    structures : dict of Structure obj.
        keys=["primitive", "unitcell", "supercell", "supercell3"]
    
    trans_matrices : dict of arrays
        keys=["primitive", "supercell", "supercell3"]
        transformation matrices for the corresponding structures

    kpts_suggested : dict of arrays
        keys=["primitive", "unitcell", "supercell", supercell3]
        kpoints suggested based on the given structures

    kpts_prev : dict of arrays
        keys=["relax", "nac", "harm", "cube"]
        kpoints which have been used for previous calculations
    
    """
    structures = None
    trans_matrices = None
    kpts_all = None
    cell_types = None
    if dir_phdb is not None:
        """ Case 1: Phonondb directory is given.
        Read data in the given Phonondb directory and suggest parameters for
        FC3.
        
        Note
        -----
        For Phonondb, the conventional cell is used for both of the
        relaxation and NAC calculations. Auto-kappa, however, can accept
        different cell types while the default types are the conventional and 
        primitive cells for the relaxation and NAC, respectively. Therefore,
        the k-meshes used for Phonondb basically will not be used for the
        automation calculation.
        
        Example of the contents in Phonondb directory
        ----------------------------------------------
        >>> BORN       FORCE_SETS   INCAR-nac    KPOINTS-force  KPOINTS-relax
        >>> phonon.yaml   POSCAR-unitcell        disp.yaml  INCAR-force
        >>> INCAR-relax  KPOINTS-nac    PAW_dataset.txt  phonopy.conf
        >>> POSCAR-unitcell.yaml

        """
        ### Read Phonondb directory
        unitcell, trans_matrices, _, nac = read_phonondb(dir_phdb)
        
        trans_matrices["unitcell"] = np.identity(3).astype(int)
        
        ### Suggest the structure for FC3
        #if max_natoms3 == max_natoms:
        #    trans_matrices["supercell3"] = trans_matrices["supercell"].copy()
        #
        #else:
        #    from auto_kappa.structure.supercell import estimate_supercell_matrix
        #    trans_matrices["supercell3"] = estimate_supercell_matrix(
        #            unitcell, max_num_atoms=max_natoms3)
        
        ### Set structures
        structures = _make_structures(
                unitcell,
                primitive_matrix=trans_matrices['primitive'],
                supercell_matrix=trans_matrices['supercell'],
                #supercell_matrix3=trans_matrices['supercell3'],
                )
        
        ### get suggested k-mesh
        kpts_suggested = {
                "primitive": klength2mesh(
                    k_length, structures["primitive"].cell.array),
                "unitcell": klength2mesh(
                    k_length, structures["unitcell"].cell.array),
                "supercell": klength2mesh(
                    k_length, structures["supercell"].cell.array),
                #"supercell3": klength2mesh(
                #    k_length, structures["supercell3"].cell.array),
                }
    
    elif file_structure is not None:
        """ Case 2: only a structure is given.
        Every required parameters are suggested.

        """
        from auto_kappa.cui.suggest import suggest_structures_and_kmeshes
        
        structures, trans_matrices, kpts_suggested = (
                suggest_structures_and_kmeshes(
                        filename=file_structure,
                        max_natoms=max_natoms,
                        #max_natoms3=max_natoms3,
                        k_length=k_length,
                        )
                    )
        
        ### This part can be modified. So far, NAC is considered for materials
        ### which is not included in Phonondb.
        nac = 2
    
    else:
        """ Case 3: error
        """
        msg = "\n Error: --directory or --file_structure must be given."
        logger.error(msg)
        sys.exit()
    
    ### Cell type used for the relaxation
    ### primitive or unitcell (conventional)
    ### All the elements of ``cell_types`` must be given.
    cell_type_relax = _get_celltype4relaxation(
            celltype_relax_given, base_directory,
            natoms_prim=len(structures['primitive'])
            )
    
    cell_types = {
            "relax": cell_type_relax,
            "nac": "primitive",
            "harm": "supercell",
            "cube": "supercell",
            }
    
    ### get previously used transformation matrices and kpoints
    params_prev = _get_previously_used_parameters(
            base_directory, structures["unitcell"].cell.array,
            cell_types=cell_types)
    
    ### prepare dictionary to compare with ``params_prev``
    params_suggest = {}
    for cal_type in params_prev:
        cell_type = cell_types[cal_type]
        params_suggest[cal_type] = {}
        params_suggest[cal_type]['trans_matrix'] = trans_matrices[cell_type]
        params_suggest[cal_type]['kpts'] = kpts_suggested[cell_type]
    
    ### determine k-mesh info for all calculations
    prioritize_previous_kpts = True
    kpts_all = _determine_kpoints_for_all(
            params_suggest, params_prev,
            prioritize_previous_kpts=prioritize_previous_kpts
            )
    
    return cell_types, structures, trans_matrices, kpts_all, nac
    
def _write_parameters(
        outfile, unitcell, cell_types, trans_matrices, kpts_used, nac):
    """ Output parameters.yaml
    """
    #logger = logging.getLogger(__name__)
    
    params = {}
    
    ### unitcell
    params["lattice"] = unitcell.cell.array.tolist()
    params["symbols"] = unitcell.get_chemical_symbols()
    params["positions"] = unitcell.get_scaled_positions().tolist()
    
    ### cell types
    params["cell_types"] = cell_types
    
    ### transformation matrix
    params["trans_matrices"] = {}
    for key in trans_matrices:
        params["trans_matrices"][key] = np.asarray(trans_matrices[key]).tolist()
    
    ### k-mesh
    params["kpts"] = {}
    for key in kpts_used:
        params["kpts"][key] = np.asarray(kpts_used[key]).tolist()
    
    ### prepare the directory
    with open(outfile, "w") as f:
        yaml.dump(params, f)
    
    #msg = " Output %s" % outfile
    #logger.info(msg)

def _start_larger_supercell(almcalc):
        
    ### print
    line = "Analyze with a larger supercell (%dx%dx%d)" % (
            almcalc.scell_matrix[0][0],
            almcalc.scell_matrix[1][1],
            almcalc.scell_matrix[2][2])
    msg = "\n"
    msg += " ###" + "#" * (len(line)) + "###\n"
    msg += " ## " + " " * (len(line)) + " ##\n"
    msg += " ## " + line              + " ##\n"
    msg += " ## " + " " * (len(line)) + " ##\n"
    msg += " ###" + "#" * (len(line)) + "###"
    
    #logger = logging.getLogger(__name__)
    logger.info(msg)
    
    msg = "\n"
    msg += " Number of atoms : %d\n" % len(almcalc.supercell)
    msg += " Cell size : %.2f, %.2f, %.2f\n" % (
            np.linalg.norm(almcalc.supercell.cell[0]),
            np.linalg.norm(almcalc.supercell.cell[1]),
            np.linalg.norm(almcalc.supercell.cell[2]),
            )
    logger.info(msg)

def analyze_harmonic_with_larger_supercells(
        almcalc_orig, 
        base_dir=None,
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

    from auto_kappa.structure.supercell import estimate_supercell_matrix
    
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
        
        _start_larger_supercell(almcalc_new)
        
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
        if almcalc_new.minimum_frequency > negative_freq:
            ak_log.negative_frequency(almcalc_new.minimum_frequency)
            return almcalc_new
        
        isc += 1
        if isc == max_loop:
            break
        
        ### prepare for the next step
        sc_mat_prev = almcalc_new.scell_matrix
    
    return almcalc_new
    
#def _get_optimal_nac(almcalc, neclect_log=None, **args):
#    """ Analyze harmonic properties with diffrent methods for NAC when negative 
#    frequencies were obtained. """
#    
#    nac_orig = almcalc.nac
#    for nac in [2, 1, 3]:
#        if nac == nac_orig:
#            continue
#        almcalc.nac = nac
#        for propt in ["band", "dos"]:
#            
#            outdir = almcalc.out_dirs["harm"]["bandos"] + "/nac_%d" % nac
#            fcsxml = "../../../result" + almcalc.outfiles["harm_xml"]
#            
#            _analyze_each_harm_propt(
#                    propt, almcalc, outdir=outdir, fcsxml=fcsxml, **args)
#            
#            print(nac)
#            exit()
#    print(" HEEEEEEEEEEEEEEEEEEEEEERE")

#def _analyze_each_harm_propt(
#        propt, almcalc, outdir=None, fcsxml=None,
#        neglect_log=False, max_num_corrections=None,
#        deltak=None, reciprocal_density=None):
#    """ Analyze each harmonic property 
#    
#    Args
#    ------
#    propt : string
#        "fc2", "band", or "dos"
#    
#    """
#    ### make an input script for the ALAMODE job
#    if propt == "band":
#        almcalc.write_alamode_input(propt=propt, deltak=deltak, outdir=outdir)
#    
#    elif propt == "dos":
#        nkts = get_automatic_kmesh(
#                almcalc.primitive, 
#                reciprocal_density=reciprocal_density)
#        almcalc.write_alamode_input(propt=propt, nkts=nkts, outdir=outdir)
#    
#    elif propt == "fc2":
#        almcalc.write_alamode_input(propt=propt, outdir=outdir)
#    
#    else:
#        msg = "\n Error: %s is not supported." % propt
#        logger.error(msg)
#    
#    ### get log file name
#    if propt == "fc2":
#        logfile = almcalc.out_dirs['harm']['force'] + "/%s.log" % propt
#    else:
#        logfile = almcalc.out_dirs['harm']['bandos'] + "/%s.log" % propt
#    
#    count = 0
#    ncores = almcalc.commands['alamode']['ncores']
#    has_error = False
#    while True:
#    
#        if count == 0:
#            neglect_log = 0
#        else:
#            neglect_log = 1
#        
#        ### calculate phonon property
#        almcalc.run_alamode(propt=propt, neglect_log=neglect_log, outdir=outdir)
#        count += 1 
#        
#        ### check log file
#        flag = should_rerun_alamode(logfile) 
#        
#        ### Check phonon band
#        ## This part solves a bug in the old calculation.
#        ## In old version, eigenvalues has sometimes contained ``nan`` values.
#        if propt == "band":
#            file_band = (
#                    almcalc.out_dirs['harm']['bandos'] + 
#                    "/%s.bands" % almcalc.prefix)
#            flag = _should_rerun_band(file_band)
#        
#        ### the property has been calculated properly.
#        if flag == False:
#            break
#        
#        ### modify the MPI and OpenMPI conditions or finish the job
#        ### Note: ALAMODE job sometimes exceeds the memory.
#        has_error = False
#        if count == 1:
#            
#            ### Change to OpenMP parallelization
#            ak_log.rerun_with_omp()
#            almcalc.commands['alamode']['anphon_para'] == "omp"
#        
#        elif count > 1 and count < max_num_corrections:
#            
#            ### modify the number of threads (even number may be better?)
#            ncores /= 2
#            if ncores > 1:
#                ncores += int(ncores % 2)
#            
#            almcalc.commands['alamode']['ncores'] = ncores
#            msg = "\n Rerun the ALAMODE job with "\
#                    "OMP_NUM_THREADS/SLURM_CPUS_PER_TASK = %d" % ncores
#            logger.error(msg)
#            
#            if ncores == 1:
#                count = max_num_corrections
#        
#        else:
#            has_error = True
#            break
#    
#    if has_error:
#        msg = "\n Error: ALAMODE job for %s has not been finished properly." % propt
#        msg += "\n Abort the job."
#        logger.error(msg)
#        sys.exit()

def analyze_harmonic_properties(
        almcalc, calculator, 
        neglect_log=False, max_num_corrections=5,
        deltak=0.01, reciprocal_density=1500,
        negative_freq=-1e-3
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
    #logger = logging.getLogger(__name__)
    
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
        from auto_kappa.io.vasp import get_dfset
        for propt in ['cv', 'lasso']:
            almcalc.write_alamode_input(propt=propt, order=2)
            almcalc.run_alamode(propt, order=2, neglect_log=neglect_log)
    else: 
        ## ver.1 : with ALM library
        ##almcalc.calc_anharm_force_constants()
        ## ver.2: with alm command
        almcalc.write_alamode_input(propt='fc3')
        almcalc.run_alamode(propt='fc3', neglect_log=neglect_log)

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
                propt='kappa', neglect_log=neglect_log, outdir=outdir, **kwargs)
        
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
    
    #logger = logging.getLogger(__name__)
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
    
def analyze_phonon_properties(
        almcalc,
        calc_force=None,
        negative_freq=-1e-3,
        material_name=None,
        neglect_log=False,
        #
        harmonic_only=False,
        calc_kappa=True,
        #
        nmax_suggest=None,
        frac_nrandom=None,
        ):
    """ Analyze phonon properties
    """
    #logger = logging.getLogger(__name__)
    
    ###
    ### 1. Calculate harmonic FCs and harmonic phonon properties
    ###
    analyze_harmonic_properties(almcalc, calc_force, negative_freq=negative_freq)
    
    ### Check negative frequency
    if almcalc.minimum_frequency < negative_freq:
        
        log = AkLog(material_name)
        log.write_yaml()
        
        ak_log.negative_frequency(almcalc.minimum_frequency)
        
        return -1
    
    if harmonic_only:
        msg = "\n"
        msg += " Harmonic properties have been calculated.\n"
        logger.info(msg)
        return 1
    
    ###
    ### 2. Calculate cubic FCs
    ###
    calculate_cubic_force_constants(
        almcalc, calc_force,
        nmax_suggest=nmax_suggest, 
        frac_nrandom=frac_nrandom, 
        neglect_log=neglect_log)
    
    ###
    ### 3. Calculate kappa with different k-mesh densities
    ###
    if calc_kappa == False:
        return 2
    
    calculate_thermal_conductivities(almcalc, 
        kdensities=[500, 1000, 1500], 
        neglect_log=neglect_log,
        temperatures_for_spectral="300:500"
        )
    
    ### output log.yaml 
    log = AkLog(material_name)
    log.write_yaml()

def _plot_bandos_for_different_sizes(
        almcalc1, almcalc2, figname="fig_bandos.png"
        ):
    
    ### Plot band and DOS
    from auto_kappa.plot.bandos import plot_bandos
    lab1 = ""
    lab2 = ""
    for j in range(3):
        lab1 += "%d" % (almcalc1.scell_matrix[j][j])
        lab2 += "%d" % (almcalc2.scell_matrix[j][j])
        if j != 2:
            lab1 += "x"
            lab2 += "x"
    
    plot_bandos(
            directory=almcalc1.out_dirs["harm"]["bandos"], 
            prefix=almcalc1.prefix,
            directory2=almcalc2.out_dirs["harm"]["bandos"], 
            prefix2=almcalc2.prefix,
            fig_labels=[lab1, lab2],
            figname=figname,
            )

def _use_omp_for_anphon(base_dir):
    from auto_kappa.alamode.log_parser import exceed_memory
    line1 = base_dir + "/harm/bandos/dos.log"
    line2 = base_dir + "/*/harm/bandos/dos.log"
    line3 = base_dir + "/cube/kappa*/kappa.log"
    for line in [line1, line2, line3]:
        fns = glob.glob(line)
        if len(fns) > 0:
            for fn in fns:
                if exceed_memory(fn):
                    return True
    return False
    
def main():
    
    options = get_parser()
    ak_params = eval(str(options))
    
    ### Get the name of the base directory
    base_dir = _get_base_directory_name(
            ak_params['material_name'], restart=ak_params['restart']
            )
    os.makedirs(base_dir, exist_ok=True)
    
    ### set logger
    logfile = base_dir + "/ak.log"
    ak_log.set_logging(filename=logfile, level=logging.DEBUG, format="%(message)s")
    #logger = logging.getLogger(__name__)
    
    ### Start auto-kappa
    start_autokappa()
    
    #msg = "\n Current directory: %s" % os.getcwd()
    #logger.info(msg)

    ### Set output directories
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
    
    ### memory check
    if _use_omp_for_anphon(base_dir):
        if ak_params["anphon_para"] != "omp":
            msg = "\n"
            msg += " Change anphon_para option to \"omp\".\n"
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
            _get_required_parameters(
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
    
    ### print parameters
    print_options(ak_params)
    
    print_conditions(
            cell_types=cell_types, 
            trans_matrices=trans_matrices,
            kpts_all=kpts_used,
            )
    
    ### write file
    os.makedirs(out_dirs["result"], exist_ok=True)
    filename = out_dirs["result"] + "/parameters.yaml"
    _write_parameters(
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
            )
    
    ### Stop because of the change of symmetry during the structure opt.
    if out < 0:
        if out == -1:
            msg = "\n Error: crystal symmetry was changed during the "\
                    "relaxation calculation."
        elif out == -2:
            msg = "\n Error: too many number of errors for the relaxation "\
                    "calculation."
        ##
        msg += "\n"
        msg += "\n STOP THE CALCULATION"
        time = datetime.datetime.now()
        msg += "\n at " + time.strftime("%m/%d/%Y %H:%M:%S")
        msg += "\n"
        logger.error(msg)
        sys.exit()
        
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

    ### Born effective charge
    if nac:
        mode = 'nac'
        apdb.run_vasp(mode,
                out_dirs[mode],
                kpts_used["nac"],
                print_params=True
                )
        
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
            cutoff3=ak_params['cutoff3'],
            magnitude=0.01,
            magnitude2=ak_params['magnitude2'],
            nac=nac,
            commands={'alamode': command_alamode, 'vasp': command_vasp},
            verbosity=ak_params['verbosity'],
            yamlfile_for_outdir=yaml_outdir
            )
    
    ### ASE calculator for forces
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
            )
    
    ### Calculate PES
    #if (almcalc.minimum_frequency < ak_params['negative_freq'] 
    #        and ak_params['pes'] == 2):
    if (almcalc.minimum_frequency < 10 and ak_params['pes'] == 2):
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
        figname = almcalc.out_dirs["result"] + "/fig_bandos.png"
        _plot_bandos_for_different_sizes(
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
    
    print_times(times, labels)
        
    ### END of calculations
    end_autokappa()

