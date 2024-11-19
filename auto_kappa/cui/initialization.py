#
# initialization.py
#
# This script helps to initialize the automation calculation based on the given
# parameters as well as parameters used for the previous calculation.
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

from auto_kappa.io.phonondb import Phonondb
from auto_kappa import output_directories
from auto_kappa.cui.suggest import klength2mesh

import logging
logger = logging.getLogger(__name__)

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
    from auto_kappa.structure.crystal import get_standardized_structure_spglib
    
    cell_type = None
    
    ### check previous calculations
    dir_relax = base_dir + "/" + output_directories["relax"]
    filenames = [dir_relax + "/vapsun.xml", dir_relax + "/freeze-1/vasprun.xml"]
    
    for fn in filenames:
        if os.path.exists(fn):
            atoms = ase.io.read(fn, format='vasp-xml')
            if len(atoms) == natoms_prim:
                cell_type = "primitive"
            elif len(atoms) > natoms_prim:
                cell_type = "unitcell"
            else:
                msg = " Error: cell type cannot be deteremined."
                msg += " Stop the calculation."
                logger.error(msg)
                sys.exit()
    
    ### determine ``cell_type`` ("unitcell" or "primitive")
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

def get_base_directory_name(label, restart=True):
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

def get_required_parameters(
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
                )
        
        ### get suggested k-mesh
        kpts_suggested = {
                "primitive": klength2mesh(
                    k_length, structures["primitive"].cell.array),
                "unitcell": klength2mesh(
                    k_length, structures["unitcell"].cell.array),
                "supercell": klength2mesh(
                    k_length, structures["supercell"].cell.array),
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

def get_previous_nac(base_dir):
    """ Get previously-used NAC parameter. If it cannot be found, return None.
    """
    ### Check bug during the VASP job for NAC
    dir_nac = base_dir + "/" + output_directories["nac"]
    file_err = dir_nac + "/std_err.txt"
    if os.path.exists(file_err):
        lines = open(file_err, 'r').readlines()
        for line in lines:
            ### If the previous VASP job was stopped due to a bug,
            if "Please submit a bug report." in line:
                msg = ("\n The previous VASP calculation in %s was aborted "
                        "due to a bug." % (dir_nac))
                msg += f'\n See {file_err} for details.'
                logger.info(msg)
                return 0
    
    ### Check the optimal NAC in previous calculations.
    dir_bandos = base_dir + "/" + output_directories["harm"]["bandos"]
    for mode in ["band", "dos"]:
        logfile = dir_bandos + "/%s.log" % mode
        if os.path.exists(logfile) == False:
            continue
        try:
            lines = open(logfile, 'r').readlines()
            for line in lines:
                if "NONANALYTIC =" in line:
                    data = line.split()
                    prev_nac = int(data[2])
                    return prev_nac
        except Exception:
            pass
    ###
    return None

def use_omp_for_anphon(base_dir):
    """ Read log files for previous calculations and check the memory.
    """
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

