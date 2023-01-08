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
import warnings
import numpy as np

import ase.io
from ase.calculators.vasp import Vasp
from ase.calculators.vasp.create_input import GenerateVaspInput

from ase.build import make_supercell
from auto_kappa.structure.crystal import (
        get_primitive_matrix, get_primitive_structure, 
        get_standardized_structure,
        change_structure_format,
        )
from auto_kappa.calculator import run_vasp
from auto_kappa.io.vasp import print_vasp_params, wasfinished
from auto_kappa.structure.crystal import get_spg_number

class ApdbVasp():
    
    def __init__(self, unitcell, 
            primitive_matrix=None,
            scell_matrix=None,
            encut_scale_factor=1.3,
            #auto_lreal_scell_size=False,
            command={'mpirun': 'mpirun', 'nprocs': 2, 'vasp': 'vasp'},
            ):
        """
        Args
        -------
        primitive, unitcell, supercell : ASE Atoms objects
        
        """
        ### matrices
        self._primitive_matrix = primitive_matrix
        self.scell_matrix = scell_matrix

        ### prepare dictionaries for structures
        self.original_structure = {'prim': None, 'unit': None, 'scell': None}
        
        self.original_structure['prim'] = \
                get_primitive_structure(unitcell, self.primitive_matrix)
        
        self.original_structure['unit'] = unitcell.copy()
        
        if scell_matrix is not None:
            self.original_structure['scell'] = \
                    make_supercell(unitcell, scell_matrix)
        
        ###
        self.relaxed_structure = {
                'prim': None,
                'unit': None,
                'scell': None,
                }
        
        ### VASP command
        self._command = command
        
        ### parameters
        self.encut_factor = encut_scale_factor
        #self.lreal_size = auto_lreal_scell_size
    
    @property
    def primitive_matrix(self):
        if self._primitive_matrix is None:
            prim_mat = get_primitive_matrix(self.unitcell)
            self._primitive_matrix = prim_mat
        return self._primitive_matrix
    
    @property
    def command(self):
        return self._command

    def update_command(self, val):
        self._command.update(val)

    @property
    def primitive(self):
        """ Relaxed structures are returned in prior to original structures.
        
        """
        if self.relaxed_structure['prim'] is not None:
            #print(" Return relaxed structure.")
            return self.relaxed_structure['prim']
        
        elif self.original_structure['prim'] is not None:
            print("")
            print(" Return UNrelaxed structure.")
            return self.original_structure['prim']
        
        else:
            return None
    
    @property
    def unitcell(self):
        
        if self.relaxed_structure['unit'] is not None:
            
            ### relaxed unitcell    
            return self.relaxed_structure['unit']
        
        elif self.relaxed_structure['prim'] is not None:
            
            ### relaxed prim => relaxed unit
            self.relaxed_structure['unit'] = make_supercell(
                    self.relaxed_structure['prim'],
                    np.linalg.inv(self.primitive_matrix)
                    )
            
            return self.relaxed_structure['unit']
        
        elif self.original_structure['unit'] is not None:
            print(" Return UNrelaxed structure.")
            return self.original_structure['unit']
         
        else:
            return None
    
    @property
    def supercell(self):
        
        if self.relaxed_structure['scell'] is not None:
            
            ### relaxed supercell
            return self.relaxed_structure['scell']
        
        elif self.relaxed_structure['prim'] is not None:
            
            ### relaxed unit => relaxed scell
            self.relaxed_structure['scell'] = make_supercell(
                    self.unitcell, self.scell_matrix)
            
            return self.relaxed_structure['scell']
        
        elif self.original_structure['scell'] is not None:
            print(" Return UNrelaxed structure.")
            return self.original_structure['scell']
        
        else:
            return None
    
    def get_calculator(self, mode, directory=None, kpts=None):
        """
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
                #auto_lreal_scell_size=self.lreal_size,
                )
        calc.command = '%s -n %d %s' % (
                self.command['mpirun'], 
                self.command['nprocs'],
                self.command['vasp'])
        return calc
    
    def run_relaxation(self, directory: './out', kpts: None,
            standardize_each_time=True,
            volume_relaxation=False,
            cell_type='p',
            force=False, num_full=2, verbosity=1): 
        """ Perform relaxation calculation, including full relaxation 
        calculations (ISIF=3) with "num_full" times and a relaxation of atomic
        positions (ISIF=2). See descriptions for self.run_vasp for details.
        
        """
        spg_before = get_spg_number(self.primitive)
        
        ### relaxation cell type
        if cell_type[0].lower() == 'p':
            cell_type = 'primitive'
            to_primitive = True
        elif cell_type[0].lower() == 'c' or cell_type[0].lower() == 'u':
            cell_type = 'conventional'
            to_primitive = False
        else:
            print(" Error")
            exit()
        
        ### message
        if verbosity != 0:
            line = "Structure optimization with %s cell" % cell_type
            msg = "\n\n " + line + "\n"
            msg += " " + "=" * (len(line))
            print(msg)
        
        ### Get the relaxed structure obtained with the old version
        ### For the old version, the xml file is located under ``directory``.
        if wasfinished(directory, filename='vasprun.xml'):
            filename = directory + "/vasprun.xml"
            prim = ase.io.read(filename, format='vasp-xml')
            prim_stand = get_standardized_structure(prim, to_primitive=True)
            self.update_structures(prim_stand, cell_type='prim')
            msg = "\n Already finised with the old version (single full relaxation)"
            print(msg)
            return 0
        
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
                print(msg)
             
            if ical == 0:
                structure = self.primitive
                print_params = True
            else:
                fn = dir_pre + "/CONTCAR"
                if os.path.exists(fn) == False:
                    warnings.warn(" Error: %s does not exist." % fn)
                    exit()
                
                #print("")
                #print(" Update the primitive structure:", fn)
                structure = ase.io.read(fn, format='vasp')
                print_params = False
            
            ### standardization
            if standardize_each_time:
                structure = get_standardized_structure(
                        structure, to_primitive=to_primitive,
                        )
            
            ### run a relaxation calculation
            self.run_vasp(
                    mode, dir_cur, kpts, 
                    structure=structure, force=force, 
                    print_params=print_params,
                    cell_type=cell_type,
                    verbosity=0
                    )
            
            dir_pre = dir_cur
        
        ### get and set the standardized primitive structure
        ## old
        #prim_stand = get_standardized_structure(
        #        self.primitive, to_primitive=to_primitive,
        #        )
        #self.update_structures(prim_stand, cell_type=cell_type)
        ## modified
        structure_stand = get_standardized_structure(
                self.primitive, to_primitive=to_primitive,
                )
        self.update_structures(structure_stand, cell_type=cell_type)
        
        ### strict relaxation with Birch-Murnaghan EOS
        if volume_relaxation:
            from auto_kappa.vasp.relax import StrictRelaxation
            outdir = directory + "/volume" 
            #init_struct = change_structure_format(prim_stand, format='pmg')
            init_struct = change_structure_format(structure_stand, format='pmg')
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
            self.update_structures(struct_ase, cell_type=cell_type)
        
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
            print("")
            warnings.warn(" WARRNING: The crystal symmetry was changed due to the structure relaxation")
            exit()
    
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
            method='custodian', force=False, print_params=False, verbosity=1):
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
            print(msg)
        
        os.environ['OMP_NUM_THREADS'] = str(self.command['nthreads'])
        
        if wasfinished(directory, filename='vasprun.xml') and force == False:
            print("")
            print(" The calculation had been already finished.")
         
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
            self.set_relaxed_structures(directory, cell_type=cell_type)
     
    
    def set_relaxed_structures(self, directory, cell_type=None):
        """ Set self.relaxed_structure
        directory : string
            Directory for "relax" calculation
        """
        filename = directory + '/vasprun.xml'
        
        if os.path.exists(filename) == False:
            print("")
            msg = " WARNING: %s does not exist."
            warnings.warn(msg)
            print("")
            return None
        
        else:
            print("")
            print(" Update the primitive structure:", filename)
            prim = ase.io.read(filename, format='vasp-xml', index=-1)
            self.update_structures(prim, cell_type=cell_type)
            
    
    def update_structures(self, structure, cell_type='prim'):
        """ Update stored structures
        Args
        -----
        structure : structure obj
        cell_type : string
            type of "structure": primitive or convenctional
        """
        if cell_type[0].lower() == 'p':
            primitive = structure.copy()
        elif cell_type[0].lower() == 'c' or cell_type[0].lower() == 'u':
            primitive = get_primitive_structure(structure, self.primitive_matrix)
        else:
            print(" Error")
            exit()

        self.relaxed_structure['prim'] = primitive
        
        self.relaxed_structure['unit'] = make_supercell(
                primitive, 
                np.linalg.inv(self.primitive_matrix)
                )
        
        self.relaxed_structure['scell'] = make_supercell(
                self.unitcell, 
                self.scell_matrix
                )
        

