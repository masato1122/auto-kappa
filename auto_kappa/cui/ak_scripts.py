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
import os.path
import numpy as np
import datetime
import warnings
import yaml

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
    msg = "\n"
    msg += "              ___  __                __   __       \n"
    msg += "     /\  |  |  |  /  \ __ |__/  /\  |__) |__)  /\  \n"
    msg += "    /~~\ \__/  |  \__/    |  \ /~~\ |    |    /~~\ \n"
    msg += "    ver. %s\n" % __version__
    msg += "\n"
    time = datetime.datetime.now()
    msg += " Start at " + time.strftime("%m/%d/%Y %H:%M:%S") + "\n"
    msg += "\n"
    print(msg)

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

def print_conditions(cell_types=None, trans_matrices=None, kpts_all=None):
    
    if cell_types is not None:
        print(" Cell type")
        print(" ---------")
        for cal_type in cell_types:
            print(" %6s : %10s" % (cal_type, cell_types[cal_type]))
        print("")
    
    ###
    if trans_matrices is not None:
        print(" Transformation matrix")
        print(" ---------------------")
        for cell_type in trans_matrices:
            print(" %10s : " % cell_type, end="")
            for i in range(3):
                vec = trans_matrices[cell_type][i]
                if cell_type == "primitive":
                    print("%.3f " * 3 % tuple(vec), end=" ")
                else:
                    print("%d " * 3 % tuple(vec), end=" ")
            print("")
        print("")
            
    ###
    if kpts_all is not None:
        print(" k-mesh")
        print(" ------")
        for cal_type in kpts_all:
            print(" %6s :" % (cal_type), " %d" * 3 % tuple(kpts_all[cal_type]))
        print("")    
    print("")

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
    from auto_kappa.io.vasp import wasfinished
    from auto_kappa.structure.crystal import get_standardized_structure_spglib
    
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
                warnings.warn(" Error: cell type cannot be deteremined.")
                exit()
    
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
        """ Case 1: when the job can be restarted
        """
        outdir += label
    
    else:
        """ Case 2: when the job cannot be restarted
        """
        if os.path.exists(label) == False:
            
            outdir += label
        
        else:
            for i in range(2, 100):
                outdir = label + "-%d" % i
                if os.path.exists(outdir) == False:
                    break
    
    return outdir

def _get_previously_used_parameters(outdir, lattice_unit):
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
                structure = Poscar.from_file(fn_structure).structure
                lattice = structure.lattice.matrix

                params_prev[label]["trans_matrix"] = \
                        np.rint(np.dot(
                            lattice, 
                            np.linalg.inv(lattice_unit)).T
                            ).astype(int)
    
    return params_prev

def _make_structures(unitcell, 
        primitive_matrix=None, supercell_matrix=None, supercell_matrix3=None):
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
    
    if supercell_matrix3 is not None:
        structures["supercell3"] = change_structure_format(
                get_supercell(unit_pp, supercell_matrix3), format='ase')
    
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
            print(" Error: cannot obtain k-mesh info properly.")
            print(" k-mesh for %s is None." % cal_type)
            exit()
    
    return kpts_all

def _get_required_parameters(
        base_directory=None,
        dir_phdb=None, file_structure=None,
        max_natoms=None, max_natoms3=None,
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
        if max_natoms3 == max_natoms:
            trans_matrices["supercell3"] = trans_matrices["supercell"]
        
        else:
            from auto_kappa.structure.supercell import estimate_supercell_matrix
            trans_matrices["supercell3"] = estimate_supercell_matrix(
                    unitcell, max_num_atoms=max_natoms3)
        
        ### Set structures
        structures = _make_structures(
                unitcell,
                primitive_matrix=trans_matrices['primitive'],
                supercell_matrix=trans_matrices['supercell'],
                supercell_matrix3=trans_matrices['supercell3'],
                )
        
        ### get suggested k-mesh
        from auto_kappa.cui.suggest import klength2mesh
        kpts_suggested = {
                "primitive": klength2mesh(
                    k_length, structures["primitive"].cell.array),
                "unitcell": klength2mesh(
                    k_length, structures["unitcell"].cell.array),
                "supercell": klength2mesh(
                    k_length, structures["supercell"].cell.array),
                "supercell3": klength2mesh(
                    k_length, structures["supercell3"].cell.array),
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
                        max_natoms3=max_natoms3,
                        k_length=k_length,
                        )
                    )
        
        ### This part can be modified. So far, NAC is performed for materials
        ### which is not included in Phonondb.
        nac = 1
    
    else:
        """ Case 3: error
        """
        warnings.warn(" Error: --directory or --file_structure must be given.")
        exit()
    
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
            "cube": "supercell3",
            }
    
    ### get previously used transformation matrices and kpoints
    params_prev = _get_previously_used_parameters(
            base_directory, structures["unitcell"].cell.array)
    
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
    
    print(" Output", outfile)


def main():
    
    times = {}
    
    options = get_parser()
    
    start_autokappa()
    
    ### Get the name of the base directory
    if options.restart == 1:
        restart = True
    else:
        restart = False
    
    base_dir = _get_base_directory_name(options.material_name, restart=restart)
    
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
    
    ######################
    ### Adjust options ###
    ######################
    ### max_natoms3: maximun limit of the number of atoms for FC3
    if options.max_natoms3 is None:
        options.max_natoms3 = options.max_natoms
    
    ### relaxed_cell
    if options.relaxed_cell is not None:
        if (options.relaxed_cell.lower()[0] == "u" or 
                options.relaxed_cell.lower()[0] == "c"):
            options.relaxed_cell = "unitcell"
        elif options.relaxed_cell.lower()[0] == "p":
            options.relaxed_cell = "primitive"
    
    ### Get required parameters for the calculation!!!!!
    (
            cell_types, structures, 
            trans_matrices, kpts_used, nac
            ) = _get_required_parameters(
                    base_directory=base_dir,
                    dir_phdb=options.directory, 
                    file_structure=options.file_structure,
                    max_natoms=options.max_natoms, 
                    max_natoms3=options.max_natoms3,
                    k_length=options.k_length,
                    celltype_relax_given=options.relaxed_cell,
                    )
    
    ### print parameters
    print_options(options)
    
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
    
    ### For materials with negative frequencies, options.max_natoms3 can be
    ### different from options.max_natoms.
    if options.max_natoms != options.max_natoms3:
        warnings.warn(" Error: max_natoms != max_natoms3 is not supported yet.")
        exti()
    
    ### command to run VASP jobs
    command_vasp = {
            'mpirun': options.mpirun, 
            'nprocs': options.ncores, 
            'nthreads': 1, 
            'vasp': options.command_vasp,
            }
    
    ### Set ApdbVasp object
    apdb = ApdbVasp(
            structures["unitcell"],
            primitive_matrix=trans_matrices["primitive"],
            scell_matrix=trans_matrices["supercell"],
            command=command_vasp,
            )
    
    ### Relaxation calculation
    mode = 'relax'
    if options.volume_relaxation == 0:
        flag = False
    else:
        flag = True
    
    apdb.run_relaxation(
            out_dirs[mode],
            kpts_used["relax"],
            volume_relaxation=flag,
            cell_type=cell_types["relax"]
            )

    ### Born effective charge
    if nac == 1:
        mode = 'nac'
        apdb.run_vasp(mode,
                out_dirs[mode],
                kpts_used["nac"],
                print_params=True
                )
    
    ### Set AlmCalc
    commands = {
            'alamode': {
                'mpirun': options.mpirun, 
                'anphon_para': options.anphon_para, 
                'ncores': options.ncores, 
                'anphon': options.command_anphon,
                'alm': options.command_alm,
                },
            'vasp': command_vasp,
            }
    
    almcalc = AlamodeCalc(
            structures["primitive"],
            base_directory=base_dir,
            restart=options.restart,
            primitive_matrix=trans_matrices["primitive"],
            scell_matrix=trans_matrices["supercell"],
            cutoff2=-1, cutoff3=options.cutoff3, 
            magnitude=0.01,
            magnitude2=options.magnitude2,
            nac=nac,
            commands=commands,
            verbosity=options.verbosity
            )
    
    ### ASE calculator for forces
    mode = 'force'
    calc_force = apdb.get_calculator(
            mode,
            kpts=kpts_used["harm"],
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
    almcalc.run_alamode(propt='evec_commensurate', neglect_log=1)
    
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
    
    ### calculate kappa with different k-mesh densities
    for kdensity in [500, 1000, 1500]:
        kpts = get_automatic_kmesh(
                almcalc.primitive, reciprocal_density=kdensity)
        outdir = (
                out_dirs['cube']['kappa_%s' % almcalc.fc3_type] + 
                "_%dx%dx%d" % (int(kpts[0]), int(kpts[1]), int(kpts[2]))
                )
        almcalc.write_alamode_input(
                propt='kappa', order=2, kpts=kpts, outdir=outdir)
        almcalc.run_alamode(
                propt='kappa', neglect_log=neglect_log, outdir=outdir)
    
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

