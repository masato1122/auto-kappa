# -*- coding: utf-8 -*-
import os.path
import os
import warnings
import numpy as np

import ase.io
from ase.calculators.vasp import Vasp
from ase.calculators.vasp.create_input import GenerateVaspInput

from ase.build import make_supercell
from .structure.crystal import (
        get_primitive_matrix, get_primitive_structure, change_structure_format)
from .calculator import run_vasp
from .io.vasp import print_vasp_params

class ApdbVasp():
    
    def __init__(self, unitcell, 
            primitive_matrix=None,
            scell_matrix=None,
            encut_scale_factor=1.3,
            auto_lreal_scell_size=65,
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
        self.lreal_size = auto_lreal_scell_size
    
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
        from .calculator import get_vasp_parameters
        
        calc = Vasp(setups='recommended', xc='pbesol')
        
        ### get structure (Atoms obj)
        if mode.lower() == 'relax' or mode.lower() == 'nac':
            structure = self.primitive
        elif mode.lower() == 'force' or mode.lower() == 'md':
            structure = self.supercell
        
        calc = get_vasp_parameters(mode, 
                directory=directory, 
                atoms=structure,
                kpts=kpts,
                encut_scale_factor=self.encut_factor,
                auto_lreal_scell_size=self.lreal_size,
                )
        calc.command = '%s -n %d %s' % (
                self.command['mpirun'], 
                self.command['nprocs'],
                self.command['vasp'])
        return calc
    
    def run_vasp(self, mode: None, directory: './out', kpts: None, 
            method='custodian', force=False, print_params=False):
        """ Run relaxation calculation
        
        Args
        -------
        mode : string
            "relax", "force", "nac", or "md"
        
        directory : string
            output directory
        
        kpts : array of float, shape=(3,)
        
        method : string
            "custodian" or "ase"

        force : bool, default=False
            If it's True, the calculation will be done forcelly even if it had
            been already finished.
            
        """
        from .io.vasp import wasfinished
        print("")
        msg = "VASP calculation (%s)" % (mode)
        border = "=" * (len(msg) + 2)
        print("", msg)
        print(border)


        if wasfinished(directory, filename='vasprun.xml') and force == False:
            print("")
            print(" The calculation had been already finished.")
            print("")
         
        else:
            calc = self.get_calculator(
                    mode.lower(), directory=directory, kpts=kpts
                    )
            
            ##
            if print_params:
                print_vasp_params(calc.asdict()['inputs'])
            
            os.environ['OMP_NUM_THREADS'] = str(self.command['nthreads'])
            run_vasp(calc, self.primitive, method=method)
            os.environ.pop("OMP_NUM_THREADS", None)
        
        ### Read the relaxed structure
        if mode.lower() == 'relax':
            
            self.set_relaxed_structures(directory)

    def set_relaxed_structures(self, directory):
        """ Set self.relaxed_structure
        directory : string
            Directory for "relax" calculation
        """
        filename = directory + '/vasprun.xml'
        
        if os.path.exists(filename) == False:
            print("")
            msg = " WARNING: %s calculation might not work properly." % mode
            warnings.warning(msg)
            print("")
            return None
        
        else:
            print(" Set relaxed structures.")
            prim = ase.io.read(filename, format='vasp-xml', index=-1)
            self.relaxed_structure['prim'] = prim
            self.relaxed_structure['unit'] = self.unitcell
            self.relaxed_structure['scell'] = self.supercell
        
