# -*- coding: utf-8 -*-
#
# apdb.py
#
# ApdbVasp class treats structures and run a few calculations such as relaxation
# calculation and Born effective charge calculation for anharmonic phonon
# database (APDB).
#
# Copyright (c) 2022 Masato Ohnishi
#
# This file is distributed under the terms of the MIT license.
# Please see the file 'LICENCE.txt' in the root directory
# or http://opensource.org/licenses/mit-license.php for information.
#
import os.path
import os
import sys
import numpy as np
import glob

import ase.io
from ase.calculators.vasp import Vasp
from ase.calculators.vasp.create_input import GenerateVaspInput

from phonopy.structure.cells import get_supercell

from auto_kappa.structure.crystal import (
        get_standardized_structure_spglib,
        change_structure_format,
        )
from auto_kappa.calculator import run_vasp, backup_vasp
from auto_kappa.io.vasp import print_vasp_params, wasfinished
from auto_kappa.structure.crystal import get_spg_number
from auto_kappa.cui import ak_log
#from auto_kappa.io.files import write_output_yaml
from auto_kappa.vasp.params import get_amin_parameter

import logging
logger = logging.getLogger(__name__)

class ApdbVasp():
    
    def __init__(
            self, unitcell, 
            primitive_matrix=None,
            scell_matrix=None,
            encut_scale_factor=1.3,
            command={'mpirun': 'mpirun', 'nprocs': 2, 'vasp': 'vasp'},
            amin_params = {},
            #yamlfile_for_outdir=None,
            ):
        """
        Args
        -------
        unitcell : structure object
            original unitcell (conventional) structure
        
        primitive_matrix : float, shape=(3,3)
            transformation matrix from the unitcell to the primitive cell with 
            the definition in Phonopy. Note that the definition in Phonopy is 
            not same as that in Pymatgen and ASE
            
        scell_matrix : float, shape=(3,3)
            transformation matrix from the unitcell to the supercell with 
            the definition in Phonopy, which is not same as Pymatgen and ASE

        Note
        ------
        Translational vectors of the primitive cell can be calculated as
        $$
        pcell = primitive_matrix.T @ unitcell.cell
        pcell = np.dot(primitive_matrix.T, unitcell.cell)
        $$
        
        """
        ### Transformation matrices
        ### The definition is the same as that for Phonopy.
        ### Please see the tutorial of Phonopy in detail.
        ### https://phonopy.github.io/phonopy/setting-tags.html#basic-tags
        self._mat_u2p = primitive_matrix
        self._mat_u2s = scell_matrix
        
        if primitive_matrix is None:
            msg = " Error: primitive_matrix must be given."
            logger(msg)
        
        if scell_matrix is None:
            msg = " Error: scell_matrix must be given."
            logger(msg)
        
        ### set structure variables
        ### Every structures will be stored in ``self._trajectory``.
        ### For example, self.trajectory[0] is the initial structures and 
        ### self.trajectory[1] is the latest structures.
        self._structures = None
        self._trajectory = []
        self.update_structures(unitcell)
        
        ### VASP command
        self._command = command
        
        ### parameters
        self.encut_factor = encut_scale_factor

        #self._yamlfile_for_outdir = yamlfile_for_outdir
        
        ### AMIN parameters
        from auto_kappa import default_amin_parameters
        self.amin_params = {}
        for key in default_amin_parameters:
            if key in amin_params.keys():
                if amin_params[key] is not None:
                    self.amin_params[key] = amin_params[key]
            if key not in self.amin_params.keys():
                self.amin_params[key] = default_amin_parameters[key]
        
    @property
    def primitive_matrix(self):
        return self._mat_u2p
    
    @property
    def scell_matrix(self):
        return self._mat_u2s
    
    @property
    def command(self):
        return self._command
    
    #@property
    #def yamlfile_for_outdir(self):
    #    return self._yamlfile_for_outdir
    
    def update_command(self, val):
        self._command.update(val)
    
    def update_structures(self, unitcell, format='ase', standardization=True):
        """ Update unit, primitive, supercells with the given new unit cell.
        Args
        -----
        unitcell : structures obj
            unit cell structure
        """
        if standardization:
            unitcell = get_standardized_structure_spglib(
                    unitcell, to_primitive=False, format=format)
        ##
        structures = self.get_structures(unitcell, format=format)
        self._structures = structures
        self._trajectory.append(structures)
    
    def get_structures(self, unitcell, format='ase'):
        """ Get primitive and supercells with the stored unitcell and 
        transformation matrices.
        """
        try:
            from phonopy import Phonopy
            phonon = Phonopy(
                    change_structure_format(unitcell, format='phonopy'),
                    self._mat_u2s,
                    primitive_matrix=self._mat_u2p
                    )
            unit = change_structure_format(phonon.unitcell , format=format) 
            prim = change_structure_format(phonon.primitive , format=format) 
            sc   = change_structure_format(phonon.supercell , format=format)
        except Exception:
            from auto_kappa.structure.crystal import get_primitive_structure_spglib
            unit = change_structure_format(unitcell , format=format) 
            prim = get_primitive_structure_spglib(unitcell)
            prim = change_structure_format(prim, format=format)
            sc = get_supercell(
                    change_structure_format(unitcell, format='phonopy'),
                    self.scell_matrix)
            sc = change_structure_format(sc, format=format)
        
        ###
        structures = {"unit": unit, "prim": prim, "super": sc}
        return structures
    
    @property
    def structures(self):
        if self._structures is not None:
            return self._structures
        else:
            return None
    
    @property
    def trajectory(self):
        return self._trajectory
    
    @property
    def primitive(self):
        return self._structures['prim']
        ##
        #if self.relaxed_structure['prim'] is not None:
        #    return self.relaxed_structure['prim']
        #
        #elif self.original_structure['prim'] is not None:
        #    return self.original_structure['prim']
        #
        #else:
        #    return None
    
    @property
    def unitcell(self):
        return self._structures['unit']
        ##
        #if self.relaxed_structure['unit'] is not None:
        #    
        #    ### relaxed unitcell    
        #    return self.relaxed_structure['unit']
        #
        #elif self.relaxed_structure['prim'] is not None:
        #    
        #    ### ver. old: may contain a critical error
        #    mat_p2u = np.rint(np.linalg.inv(self._mat_u2p)).astype(int)
        #    #
        #    ### ver. corrected
        #    #mat_p2u = np.linalg.inv(self.primitive_matrix).T
        #
        #    ### relaxed prim => relaxed unit
        #    self.relaxed_structure['unit'] = make_supercell(
        #            self.relaxed_structure['prim'], mat_p2u)
        #    
        #    return self.relaxed_structure['unit']
        #
        #elif self.original_structure['unit'] is not None:
        #    return self.original_structure['unit']
        # 
        #else:
        #    return None
    
    @property
    def supercell(self):
        return self._structures['super']
        ##
        #if self.relaxed_structure['scell'] is not None:
        #    
        #    ### relaxed supercell
        #    return self.relaxed_structure['scell']
        #
        #elif self.relaxed_structure['prim'] is not None:
        #    
        #    ### relaxed unit => relaxed scell
        #    self.relaxed_structure['scell'] = make_supercell(
        #            self.unitcell, self._mat_u2s)
        #    
        #    return self.relaxed_structure['scell']
        #
        #elif self.original_structure['scell'] is not None:
        #    return self.original_structure['scell']
        #
        #else:
        #    return None
    
    def get_calculator(self, mode, directory=None, kpts=None, **args):
        """ Return VASP calculator created by ASE
        Args
        ------
        mode : string
            'relax', 'force', 'nac', or 'md'
        
        directory : string
            output directory
        
        kpts : list of float, shpae=(3,)
        
        """
        from auto_kappa.calculator import get_vasp_calculator
        
        ### get structure (Atoms obj)
        if 'relax' in mode.lower() or mode.lower() == 'nac':
            structure = self.primitive
        elif 'force' in mode.lower() or mode.lower() == 'md':
            structure = self.supercell
        
        calc = get_vasp_calculator(mode, 
                directory=directory, 
                atoms=structure,
                kpts=kpts,
                encut_scale_factor=self.encut_factor,
                **args,
                )
        calc.command = '%s -n %d %s' % (
                self.command['mpirun'], 
                self.command['nprocs'],
                self.command['vasp'])
        return calc
    
    def run_relaxation(
            self, directory: './out', kpts: None,
            standardize_each_time=True,
            volume_relaxation=0,
            cell_type='p',
            force=False, num_full=2, verbosity=1,
            max_error=None, nsw_params=None, 
            **args
            ):
        """ Perform relaxation calculation, including full relaxation 
        calculations (ISIF=3) with "num_full" times and a relaxation of atomic
        positions (ISIF=2). See descriptions for self.run_vasp for details.
        
        Args
        =======

        directory : string
            working directory for VASP

        kpts : array, shape=(3)
            k-mesh for VASP

        standardize_each_time : bool
        
        volume_relaxation : int

        cell_type : string

        force : bool

        num_full : int,
            Number of relxation calculation w/o any restriction [default: 2]

        verbosity : int
        
        max_error : int
            Max number of retry the calculation. If error.{max_error}.tar(.gz) 
            exists, stop the calculation.
        
        args : dictionary
            input parameters for VASP
        
        Return
        ========
        
        integer :
            If negative value, stop the job.
            -1 : symmetry error
            -2 : too many errors
        
        """
        ### relaxation cell type
        if cell_type[0].lower() == 'p':
            cell_type = 'primitive'
            to_primitive = True
        elif cell_type[0].lower() == 'c' or cell_type[0].lower() == 'u':
            cell_type = 'conventional'
            to_primitive = False
        else:
            msg = " Error"
            logger.info(msg)
            sys.exit()
        
        ### message
        if verbosity != 0:
            #line = "Cell type for structure optimization: %s" % cell_type
            line = "Structure optimization"
            msg = "\n\n " + line
            msg += "\n " + "=" * (len(line))
            msg += "\n\n Cell type : %s" % cell_type
            logger.info(msg)
         
        ### Get the relaxed structure obtained with the old version
        ### For the old version, the xml file is located under ``directory``.
        if wasfinished(directory, filename='vasprun.xml'):
            filename = directory + "/vasprun.xml"
            prim = ase.io.read(filename, format='vasp-xml')
            
            mat_p2u = np.linalg.inv(self.primitive_matrix)
            mat_p2u = np.array(np.sign(mat_p2u) * 0.5 + mat_p2u, dtype="intc")
            unitcell = get_supercell(
                    change_structure_format(prim, format='phonopy'),
                    mat_p2u)
            unitcell = change_structure_format(unitcell) 
            self.update_structures(unitcell)
            msg = "\n Already finised with the old version (single full relaxation)"
            return 0
        
        ### NSW parameters
        out = _parse_nsw_params(nsw_params)
        nsw_init = out[0]
        nsw_diff = out[1]
        nsw_min = out[2]
        
        ### symmetry
        spg_before = get_spg_number(self.unitcell)
        
        ### perform relaxation calculations
        count = 0
        count_err = 0
        max_sym_err = 2
        while True:
            
            ### set working directory and mode
            if count < num_full:
                ## full relxation
                num = count + 1
                dir_cur = directory + "/full-%d" % num
                mode = 'relax-full'
            
            else:
                ## relaxation of atomic positions
                num = count - num_full + 1
                dir_cur = directory + "/freeze-%d" % num
                mode = 'relax-freeze'
            
            #### check the number of errors
            #if max_error is not None:
            #    if too_many_errors(dir_cur, max_error=max_error):
            #        return -2
            
            #### determine NSW parameter based on the number of errors
            args['nsw'] = _get_nsw_parameter(
                    dir_cur, nsw_init=nsw_init, 
                    nsw_diff=nsw_diff, nsw_min=nsw_min)
            
            ### print message
            if verbosity != 0:
                line = "%s (%d)" % (mode, num)
                msg = "\n " + line
                msg += "\n " + "-" * len(line)
                logger.info(msg)
            
            ##
            if count == 0:
                if count_err == 0:
                    print_params = True
            else:
                fn = dir_pre + "/CONTCAR"
                if os.path.exists(fn) == False:
                    msg = "\n Error: %s does not exist." % fn
                    logger.error(msg)
                    sys.exit()
                
                print_params = False
            
            ### get the structure used for the analysis
            if to_primitive:
                structure = self.primitive
            else:
                structure = self.unitcell
            
            ### run a relaxation calculation
            ### out == -1 : symmetry was changed
            out = self.run_vasp(
                    mode, dir_cur, kpts, 
                    structure=structure, force=force, 
                    print_params=print_params,
                    cell_type=cell_type,
                    verbosity=0,
                    standardization=standardize_each_time,
                    **args
                    )
            
            if out == -1:
                
                ### backup failed result
                backup_vasp(dir_cur, delete_files=True)
                
                ### set ISYM = 2 explicitly
                #args["isym"] = 2
                
                count_err += 1
                if max_sym_err == count_err:
                    msg =  "\n The calculation was failed %d times." % (count_err)
                    msg += "\n Abort the relaxation calculation."
                    logger.info(msg)
                    return -1
                else:
                    logger.info("\n Retry the relaxation calculation.")
                    continue
            
            ### update
            dir_pre = dir_cur
            
            count += 1
            count_err = 0
            if count == num_full + 1:
                break
        
        ### update structures
        self.update_structures(self.unitcell, standardization=True)
        
        ### strict relaxation with Birch-Murnaghan EOS
        if volume_relaxation:
            from auto_kappa.vasp.relax import StrictRelaxation
            outdir = directory + "/volume" 
            
            if to_primitive:
                structure = self.primitive
            else:
                structure = self.unitcell
            
            init_struct = change_structure_format(structure, format='pmg')
            relax = StrictRelaxation(init_struct, outdir=outdir)
            Vs, Es = relax.with_different_volumes(
                    kpts=kpts, 
                    command=self.command,
                    )
            figname = outdir + '/fig_bm.png'
            relax.plot_bm(figname=figname)
            relax.print_results()
            struct_opt = relax.get_optimal_structure()
            outfile = outdir + "/POSCAR.opt"
            struct_opt.to(filename=outfile)
            struct_ase = change_structure_format(struct_opt, format='ase') 
            
            ###
            if to_primitive:
                _mat_p2u = np.linalg.inv(self.primitive_matrix)
                _mat_p2u = np.array(np.sign(_mat_p2u) * 0.5 + _mat_p2u, dtype="intc")
                unitcell = get_supercell(struct_ase, _mat_p2u)
            else:
                unitcell = struct_ase.copy()
            
            self.update_structures(unitcell)
        
        ### Check the crystal symmetry before and after the relaxation
        spg_after = get_spg_number(self.primitive)
        
        self._write_relax_yaml({
            'directory': directory,
            'cell_type': cell_type,
            'structure': self.unitcell,
            'spg': [spg_before, spg_after],
            'volume_relaxation': volume_relaxation,
            })
        
        if spg_before != spg_after:
            ak_log.symmetry_error(spg_before, spg_after)
            return -1
         
        return 0
    
    def _write_relax_yaml(self, params):
        import yaml
        outfile = params['directory'] + '/relax.yaml'
        structure = change_structure_format(params['structure'], format='pymatgen') 
        
        ### lattice vectors
        lattice = []
        for v1 in structure.lattice.matrix:
            lattice.append([])
            for val in v1:
                lattice[-1].append(float(val))
        
        ### fractional coords
        frac_coord = []
        for pos in structure.frac_coords:
            frac_coord.append([])
            for j in range(3):
                frac_coord[-1].append(float(pos[j]))
        
        ### species
        species = [el.name for el in structure.species]
        
        dict_data = {
                'directory': params['directory'],
                'cell_type_for_relaxation': params['cell_type'],
                'spg_before': params['spg'][0],
                'spg_after': params['spg'][1],
                'lattice': lattice,
                'positions': frac_coord,
                'species': species,
                'volume_relaxation': params['volume_relaxation'],
                }

        with open(outfile, 'w') as f:
            yaml.dump(dict_data, f)
            
    
    def run_vasp(self, mode: None, directory: './out', kpts: None, 
            structure=None, cell_type=None,
            method='custodian', force=False, print_params=False, 
            standardization=True, verbosity=1, **args
            ):
        """ Run relaxation and born effective charge calculation
        
        Args
        -------
        mode : string
            "relax-full", "relax-freeze", "force", "nac", or "md"
        
        directory : string
            output directory
        
        kpts : array of float, shape=(3,)

        structure : structure obj

        cell_tyep : string
            cell type of ``structure``: primitive or conventional
            This is used only for ``mode = relax-***``
        
        method : string
            "custodian" or "ase"

        force : bool, default=False
            If it's True, the calculation will be done forcelly even if it had
            been already finished.
        
        args : dict
            input parameters for VASP
        
        Return
        --------
        integer :
            0. w/o error
            1. symmetry was changed during the relaxation calculation
        
        """
        if verbosity != 0:
            line = "VASP calculation (%s)" % (mode)
            msg = "\n\n " + line
            msg += "\n " + "=" * (len(line))
            logger.info(msg)
        
        ### set OpenMP
        omp_keys = ["OMP_NUM_THREADS", "SLURM_CPUS_PER_TASK"]
        for key in omp_keys:
            os.environ[key] = str(self.command['nthreads'])
        
        ### perform the calculation
        if wasfinished(directory, filename='vasprun.xml') and force == False:
            msg = "\n The calculation has been already done."
            logger.info(msg)
        
        else:
            
            #### ver.1 relax with one shot
            calc = self.get_calculator(
                    mode.lower(), directory=directory, kpts=kpts, **args
                    )
            ##
            if print_params:
                print_vasp_params(calc.asdict()['inputs'])
            
            if structure is None:
                structure = self.primitive
            
            ### get AMIN
            amin = get_amin_parameter(
                    calc.directory, structure.cell.array, **self.amin_params)
            if amin is not None:
                calc.set(amin=amin)

            ### run a VASP job
            run_vasp(calc, structure, method=method)
        
        ### set back OpenMP 
        for key in omp_keys:
            os.environ[key] = "1"
         
        ### Read the relaxed structure
        if 'relax' in mode.lower():
            
            c0 = cell_type.lower()[0]
            
            vasprun = directory + "/vasprun.xml"
            
            if c0 == 'c' or c0 == 'u':
                
                new_unitcell = ase.io.read(vasprun, format='vasp-xml')

            else:
                
                ### read primitive and transform it to the unit cell
                new_prim = ase.io.read(vasprun, format='vasp-xml')
                _mat_p2u = np.linalg.inv(self.primitive_matrix)
                _mat_p2u = np.array(np.sign(_mat_p2u) * 0.5 + _mat_p2u, dtype="intc")
                new_unitcell = get_supercell(
                        change_structure_format(new_prim, format='phonopy'),
                        _mat_p2u)
            
            ##
            num_init = get_spg_number(structure)
            num_mod = get_spg_number(new_unitcell)
            if num_init != num_mod:
                ak_log.symmetry_error(num_init, num_mod)
                return -1
            
            self.update_structures(new_unitcell, standardization=standardization)
        
        return 0

def too_many_errors(directory, max_error=100):
    """ check the number of errors in ``directory`` """    
    for file_err in glob.glob(directory+"/error.*"):
        try:
            num = int(file_err.split("/")[-1].split(".")[1])
        except Exception:
            continue
        if num >= max_error:
            return True
    return False

#def _get_number_of_errors(directory):
#    """ Get and return the number of errors in the given directory. """
#    num_errors = 0
#    ### number of errors
#    for suffix in ["tar", "tar.gz"]:
#        line = directory + "/error.*." + suffix
#        fns = glob.glob(line)
#        num_errors += len(fns)
#    ####
#    #line = directory + "/INCAR"
#    #fns = glob.glob(line)
#    #num_errors += len(fns)
#    return num_errors

def _parse_nsw_params(line, params_default=[200, 10, 20]):
    """ Return NSW params with an array 
    Args
    ======
    
    line : string, "**:**:**"

    Return
    =======
    
    array, shape=(3)
        initial, interval, and minimum NSW
    """
    data = line.split(":")
    params = []
    for j in range(3):
        try:
            params.append(int(data[j]))
        except Exception:
            params.append(int(params_default[j]))
    return params

def _get_nsw_parameter(directory, nsw_init=200, nsw_diff=10, nsw_min=20):
    """ Determine the number of NSW based on the number of errors """
    from auto_kappa.vasp.params import get_number_of_errors
    num_errors = get_number_of_errors(directory)
    nsw = max(nsw_min, nsw_init - nsw_diff * num_errors)
    return nsw

#def _get_amin_parameter(directory, lattice, len_tol=50, amin_set=0.01):
#    """ Get and return AMIN """
#    
#    num_errors = _get_number_of_errors(directory)
#    if num_errors == 0:
#        return None
#    
#    ### if error exists,
#    amin = None
#    for j in range(3):
#        length = np.linalg.norm(lattice[j])
#        if length > len_tol:
#            return amin_set
#    return None

