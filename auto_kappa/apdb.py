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

import ase.io
from ase.calculators.vasp import Vasp
from ase.calculators.vasp.create_input import GenerateVaspInput

from phonopy import Phonopy
from phonopy.structure.cells import get_supercell

from auto_kappa.structure.crystal import (
        get_standardized_structure_spglib,
        change_structure_format,
        )
from auto_kappa.calculator import run_vasp
from auto_kappa.io.vasp import print_vasp_params, wasfinished
from auto_kappa.structure.crystal import get_spg_number
from auto_kappa.cui import ak_log
#from auto_kappa.io.files import write_output_yaml

import logging
logger = logging.getLogger(__name__)

class ApdbVasp():
    
    def __init__(
            self, unitcell, 
            primitive_matrix=None,
            scell_matrix=None,
            encut_scale_factor=1.3,
            command={'mpirun': 'mpirun', 'nprocs': 2, 'vasp': 'vasp'},
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
                    unitcell, to_primitive=False, format=format
                    )
        ##
        structures = self.get_structures(unitcell, format=format)
        self._structures = structures
        self._trajectory.append(structures)
    
    def get_structures(self, unitcell, format='ase'):
        """ Get primitive and supercells with the stored unitcell and 
        transformation matrices.
        """
        phonon = Phonopy(
                change_structure_format(unitcell, format='phonopy'),
                self._mat_u2s,
                primitive_matrix=self._mat_u2p
                )
        
        unit = change_structure_format(phonon.unitcell , format=format) 
        prim = change_structure_format(phonon.primitive , format=format) 
        sc   = change_structure_format(phonon.supercell , format=format) 
        
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
    
    def get_calculator(self, mode, directory=None, kpts=None):
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
                )
        calc.command = '%s -n %d %s' % (
                self.command['mpirun'], 
                self.command['nprocs'],
                self.command['vasp'])
        return calc
    
    def run_relaxation(self, directory: './out', kpts: None,
            standardize_each_time=True,
            volume_relaxation=0,
            cell_type='p',
            force=False, num_full=2, verbosity=1): 
        """ Perform relaxation calculation, including full relaxation 
        calculations (ISIF=3) with "num_full" times and a relaxation of atomic
        positions (ISIF=2). See descriptions for self.run_vasp for details.
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
            line = "Structure optimization with %s cell" % cell_type
            msg = "\n\n " + line + "\n"
            msg += " " + "=" * (len(line))
            logger.info(msg)
         
        #### Get the relaxed structure obtained with the old version
        #### For the old version, the xml file is located under ``directory``.
        if wasfinished(directory, filename='vasprun.xml'):
            filename = directory + "/vasprun.xml"
            prim = ase.io.read(filename, format='vasp-xml')
            
            mat_p2u = np.linalg.inv(self.primitive_matrix)
            unitcell = get_supercell(
                    change_structure_format(prim, format='phonopy'),
                    mat_p2u)
            unitcell = change_structure_format(unitcell) 
            self.update_structures(unitcell)
            msg = "\n Already finised with the old version (single full relaxation)"
            return 0
        
        spg_before = get_spg_number(self.unitcell)
        #print(spg_before)
        #exit()
         
        ### perform relaxation calculations
        for ical in range(num_full+1):
            
            ### set working directory and mode
            if ical < num_full:
                num = ical + 1
                dir_cur = directory + "/full-%d" % num
                mode = 'relax-full'
            
            else:
                num = ical - num_full + 1
                dir_cur = directory + "/freeze-%d" % num
                mode = 'relax-freeze'
            
            if verbosity != 0:
                line = "%s (%d)" % (mode, num)
                msg = "\n " + line + "\n"
                msg += " " + "-" * len(line)
                
                logger.info(msg)
            
            ###
            if ical == 0:
                print_params = True
            
            else:
                fn = dir_pre + "/CONTCAR"
                if os.path.exists(fn) == False:
                    msg = "\n Error: %s does not exist." % fn
                    logger.error(msg)
                    sys.exit()
                
                #msg = ""
                #msg += " Update the primitive structure: " + fn
                #structure = ase.io.read(fn, format='vasp')
                print_params = False
            
            ### get the structure used for the analysis
            if to_primitive:
                structure = self.primitive

            else:
                structure = self.unitcell
            
            #### standardization
            #if standardize_each_time:
            #    
            #    structure = get_standardized_structure_spglib(
            #            structure, to_primitive=to_primitive, format='ase'
            #            )
            #
            
            ### run a relaxation calculation
            self.run_vasp(
                    mode, dir_cur, kpts, 
                    structure=structure, force=force, 
                    print_params=print_params,
                    cell_type=cell_type,
                    verbosity=0,
                    standardization=standardize_each_time
                    )
            
            ### output to yaml
            #name = "%s-%d" % (mode, ical)
            #info = {"directory": dir_cur.replace(,
            #        "kind": "",
            #        "note": ""}
            #write_output_yaml(self.yamlfile_for_outdir, name, info)
            #exit()
             
            ### update
            dir_pre = dir_cur

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
                unitcell = get_supercell(
                        struct_ase, 
                        np.linalg.inv(self.primitive_matrix)
                        )
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
            sys.exit()

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
            standardization=True, verbosity=1
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
        """
        if verbosity != 0:
            line = "VASP calculation (%s)" % (mode)
            msg = "\n\n " + line + "\n"
            msg += " " + "=" * (len(line)) + "\n"
            logger.info(msg)
        
        os.environ['OMP_NUM_THREADS'] = str(self.command['nthreads'])
        
        if wasfinished(directory, filename='vasprun.xml') and force == False:
            msg = "\n The calculation has been already finished."
            logger.info(msg)
        
        else:
        
            #### ver.1 relax with one shot
            calc = self.get_calculator(
                    mode.lower(), directory=directory, kpts=kpts
                    )
            ##
            if print_params:
                print_vasp_params(calc.asdict()['inputs'])
            
            if structure is None:
                structure = self.primitive
            
            run_vasp(calc, structure, method=method)
        
        os.environ["OMP_NUM_THREADS"] = "1"
        
        ### Read the relaxed structure
        if 'relax' in mode.lower():
            
            c0 = cell_type.lower()[0]
            
            vasprun = directory + "/vasprun.xml"
            
            if c0 == 'c' or c0 == 'u':
                
                new_unitcell = ase.io.read(vasprun, format='vasp-xml')

            else:
                
                ### read primitive and transform it to the unit cell
                new_prim = ase.io.read(vasprun, format='vasp-xml')
                mat_p2u = np.linalg.inv(self.primitive_matrix)
                new_unitcell = get_supercell(
                        change_structure_format(new_prim, format='phonopy'),
                        mat_p2u)
            
            ##
            num_init = get_spg_number(structure)
            num_mod = get_spg_number(new_unitcell)
            if num_init != num_mod:
                
                ak_log.symmetry_error(num_init, num_mod)
                sys.exit()
            
            self.update_structures(new_unitcell, standardization=standardization)
     

