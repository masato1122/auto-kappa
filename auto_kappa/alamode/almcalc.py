#
# almcalc.py
#
# Interface for ALAMODE
#
# Copyright (c) 2022 Masato Ohnishi
#
# This file is distributed under the terms of the MIT license.
# Please see the file 'LICENCE.txt' in the root directory
# or http://opensource.org/licenses/mit-license.php for information.
#
"""
To Do
=============
2. Use different (larger) supercell size for 3rd-order force constants

"""
# -*- coding: utf-8 -*-
import os
import os.path
import numpy as np
import warnings
import logging

import ase
import subprocess
import shutil
import pandas as pd
import f90nml
import glob

from phonopy import Phonopy
from phonopy.structure.cells import get_primitive
from phonopy.structure.cells import get_supercell

from auto_kappa import output_directories, output_files
from auto_kappa.units import AToBohr, BohrToA
from auto_kappa.structure.crystal import change_structure_format, get_formula
from auto_kappa.io.vasp import wasfinished, get_dfset, read_dfset, print_vasp_params
from auto_kappa.calculator import run_vasp

from auto_kappa.io.vasp import write_born_info
from auto_kappa.io.alm import AlmInput, AnphonInput

class AlamodeCalc():
    
    ### k1: alamode type, k2: mode, k3: propt
    propts_mode_type = {
            'anphon': {
                'phonons': ['band', 'dos', 'evec_commensurate'],
                'rta': ['kappa'],
                'scph': ['scph'],
                },
            'alm':{
                'suggest': ['suggest'],
                'optimize': ['fc2', 'fc3', 'cv', 'lasso'],
                }
            }
    
    def __init__(self, 
            prim_given, 
            base_directory='apdb_output',
            primitive_matrix=None,
            scell_matrix=None, 
            scell_matrix3=None,
            restart=1,
            cutoffs=None, nbody=[2,3,3,2,2], magnitude=0.01, magnitude2=0.03,
            cutoff2=-1, cutoff3=4.3, 
            order_lasso=5,
            nac=None,
            commands=None,
            verbosity=0,
            ):
        """ Calculations with ALM and Anphon are managed with this class.
        
        Args
        =======
        prim_given : primitive structure
            Different formats such as pymatgen and ASE are accepted while the
            format is changed to ase-Atoms in this module.

        base_directory : string
            The top directory of the output directories. For example, material
            ID of Materials Project or composition [apdb_output]
        
        restart : int, default 1
            The calculation will restart (1) or NOT restart (0) when the
            directory exists.

        scell_matrix : array of float, shape=(3,3)
            transformation matrix from the unitcell to a supercell
        
        primitive_matrix : array of float, shape=(3,3)
            transformation matrix from the unitcell to the primitive cell
        
        scell_matrix3 : array of float, shape=(3,3)
            transformation matrix from the unitcell to a supercell used to
            calculate cubic FCs

        cutoffs : shape=(order, num_elems, num_elems), unit=[Ang]
            cutoff radii for force constants. If it's not given, it will be set
            automatically with cutoff2 and cutoff3. Default cutoff radii are 
            "-1" and "4.3" Angstrom for harmonic and cubid FCs, respectively.
        
        ndoby : list of int
            nbody[i]-body interaction is considered for {i+2}-order FCs.
            If it is not given, default value is [2,3,3,2,2,...].
            See the tutorial of Alamode for details.
         
        cutoff2, cutoff3 : float, unit=[Ang]
            cutoff radii.
            If cutoffs is None, cutoff2 and cutoff3 are used.
            For lasso calculation, cuttoff3 is used for higher order FCs.
        
        nac : int
            If nac=0, nac is not considered while, if nac != 0, nac is considered.
    
        verbosity : int
            If verbosity=1, details are output while, if 0, not.

        commands : dict
            Number of processes and threads are given.
            default={
                'mpirun': 'mpirun', "anphon_para": "omp", 
                "ncores": 1, "anphon": "anphon", "alm": "alm",
                }
        
        Example
        ----------
        >>> almcalc = AlamodeCalc(
        >>>     primitive,
        >>>     base_directory='mpid-149',
        >>>     primitive_matrix=pmat,
        >>>     scell_matrix=smat,
        >>>     cutoff2=-1, cutoff3=4.3,
        >>>     nbody=[2,3,3,2], magnitude=0.01,
        >>>     nac=1,
        >>>     commands={'mpirun':'mpirun', 'anphon_para':'omp', 'ncores':1,'anphon':'anphon'},
        >>>     verbosity=0
        >>>     )
        >>> calc_force = apdb.get_calculator('force', 
        >>>     "./mp-149/harm/force", [4,4,4],
        >>>     )
        >>> almcalc.calc_forces(1, calc_force)
        >>> almcalc.calc_harmonic_force_constants()
        >>> almcalc.write_alamode_input(propt='dos')
        >>> almcalc.run_alamode(propt='dos')
        """
        ### set name of the output directory
        self.set_basedir_name(base_directory, restart)
        
        ###
        self._commands = {
                'vasp':{
                    'mpirun': 'mpirun', "nprocs": 1, "nthreads": 1, 
                    'vasp': 'vasp'
                    },
                'alamode':{
                    'mpirun': 'mpirun', 'anphon_para': 1, 'ncores': 2,
                    'alm': 'alm', 'anphon': 'anphon',
                    }
                }
        self._commands.update(commands)
        
        ### output directories
        self.out_dirs = {}
        for k1 in output_directories.keys():
            values1 = output_directories[k1]
            if type(values1) == str:
                self.out_dirs[k1] = self.base_directory + '/' + values1
            else:
                self.out_dirs[k1] = {}
                for k2 in values1.keys():
                    values2 = values1[k2]
                    self.out_dirs[k1][k2] = self.base_directory + '/' + values2
        
        ###
        self.outfiles = output_files
        
        ### matrices
        if primitive_matrix is None:
            warnings.warn(" Error: primitive_matrix must be given.")
        
        if scell_matrix is None:
            warnings.warn(" Error: scell_matrix must be given.")
        
        if scell_matrix3 is None:
            scell_matrix3 = scell_matrix.copy()
        
        self._primitive_matrix = primitive_matrix
        self._scell_matrix = scell_matrix
        self._scell_matrix3 = scell_matrix3
        
        ### make the unitcell with ``get_supercell`` in Phonopy
        unit_pp = get_supercell(
                change_structure_format(prim_given, format='phonopy'),
                np.linalg.inv(primitive_matrix)
                )
        
        ### set unitcell, primitive cell, and two kinds of supercells
        self._primitive = None
        self._unitcell = None
        self._supercell = None
        self._supercell3 = None
        
        self._set_structures(unit_pp)
        
        self._num_elems = None
        
        ### nbody and cutoffs
        ### set cutoff radii: Ang => Bohr
        nbody = None
        if nbody is None:
            self.set_nbody_automatically()
        else:
            self.nbody = nbody
        
        self._cutoff2 = cutoff2
        self._cutoff3 = cutoff3
        self._cutoffs = cutoffs
        
        ##
        self._magnitude = magnitude       ### Ang
        self._magnitude2 = magnitude2     ### Ang
        
        self.force_calculators = {}
        
        self.nac = nac
        
        self.lasso = False
        
        self._fc3_type = None            ### 'df' or 'lasso'
        
        self._nmax_suggest = None

        self._commands = commands

        #self._order_lasso = order_lasso
        
        self._file_result = None
        self._file_isotope = None
        self._scat = None
        
        ###
        self._prefix = get_formula(self.primitive)  ## prefix for input files
        self._frequency_range = None
        #self._minimum_frequency = None
        
        self.verbosity = verbosity
    
    def _set_structures(self, unit_given, format='ase'):
        """ Set every structures with the given format: primitive, unit, and two
        kinds of supercells. If supercell can be read from vasprun.xml, xml file
        will be used.
        """
        ### unit cell
        unit_pp = change_structure_format(unit_given, format="phonopy")
        self._unitcell = change_structure_format(unit_pp, format=format)
        
        ### primitive cell
        prim_pp = get_primitive(unit_pp, self.primitive_matrix)
        self._primitive = change_structure_format(prim_pp, format=format)
        
        ### supercell for harmonic FCs
        xml = self.out_dirs['harm']['force'] + '/prist/vasprun.xml'
        try:
            sc = ase.io.read(xml, format='vasp-xml')
        except Exception:
            sc = get_supercell(unit_pp, self.scell_matrix)
        
        self._supercell = change_structure_format(sc, format=format)
         
        ### supercell for cubic FCs
        xml_fd = self.out_dirs['cube']['force_fd'] + '/psist/vasprun.xml'
        xml_lasso = self.out_dirs['cube']['force_lasso'] + '/psist/vasprun.xml'
        try:
            sc3 = ase.io.read(xml_fd, format='vasp-xml')
        except Exception:
            try:
                sc3 = ase.io.read(xml_lasso, format='vasp-xml')
            except Exception:
                sc3 = get_supercell(unit_pp, self.scell_matrix3)
        
        self._supercell3 = change_structure_format(sc3, format=format)
        
    @property
    def magnitude(self):
        return self._magnitude

    @property
    def magnitude2(self):
        return self._magnitude2
    
    @property
    def fc3_type(self):
        return self._fc3_type
    
    @property
    def base_directory(self):
        return self._base_directory
    
    def set_basedir_name(self, dir0, restart):
        """
        Note 
        -----
        If restart is off, the directory name needs to be determined carefully.
        """
        if restart == 1:
            self._base_directory = dir0
        else:
            if os.path.exists(dir0):
                for i in range(2, 100):
                    dir1 = dir0 + '-%d' % i
                    if os.path.exists(dir1) == False:
                        self._base_directory = dir1
                        break
            else:
                self._base_directory = dir0
    
    @property
    def commands(self):
        return self._commands
    
    def update_commands(self, val):
        """ Command is update """
        self._commands.update(val)
    
    @property
    def prefix(self):
        if self._prefix is None:
            self._prefix = get_formula(self.primitive)
        return self._prefix

    @property
    def nbody(self):
        return self._nbody

    @nbody.setter
    def nbody(self, nbody):
        self._nbody = nbody
        
    def set_nbody_automatically(self, order=5):
        nbody = []
        for i in range(order):
            if i <= 1:
                nbody.append(i+2)
            elif i == 2:
                nbody.append(3)
            else:
                nbody.append(2)
        ## update
        self._nbody = nbody
     
    @property
    def cutoff2(self):
        return self._cutoff2
    @cutoff2.setter
    def cutoffs(self, c2):
        self._cutoff2 = c2

    @property
    def cutoff3(self):
        return self._cutoff3
    @cutoff3.setter
    def cutoffs(self, c3):
        self._cutoff3 = c3
    
    @property
    def cutoffs(self):
        if self._cutoffs is not None:
            return self._cutoffs
        else:
            ### Bohr => Angstrom
            self._cutoffs = get_cutoffs_automatically(
                    cutoff2=self.cutoff2,
                    cutoff3=self.cutoff3,
                    num_elems=self.num_elems,
                    order=len(self.nbody)
                    ) * AToBohr
        return self._cutoffs
    
    @cutoffs.setter
    def cutoffs(self, matrix):
        self._cutoffs = matrix

    @property
    def nmax_suggest(self):
        return self._nmax_suggest

    @nmax_suggest.setter
    def nmax_suggest(self, val):
        self._nmax_suggest = val

    @property
    def num_elems(self):
        if self._num_elems is None:
            self._num_elems = len(list(
                set(self.primitive.get_chemical_symbols())
                ))
        return self._num_elems

    @property
    def primitive_matrix(self):
        return self._primitive_matrix
    
    @property
    def scell_matrix(self):
        return self._scell_matrix
    
    @property
    def scell_matrix3(self):
        return self._scell_matrix3
    
    @property
    def primitive(self):
        return self._primitive

    @property
    def unitcell(self):
        if self._unitcell is None:
            warnings.warn(" Error: unitcell is not yet defined.")
        else:
            return self._unitcell
    
    @property
    def supercell(self):
        """
        """
        if self._supercell is None:
            warnings.warn(" Error: supercell is not yet defined.")
        else:
            return self._supercell

    @property
    def supercell3(self):
        """
        """
        if self._supercell3 is None:
            warnings.warn(" Error: supercell for cubic FCs is not yet defined.")
        else:
            return self._supercell3
        
    @property
    def frequency_range(self):
        if self._frequency_range is None:
            fn = self.out_dirs['harm']['bandos'] + '/' + self.prefix + '.bands'
            if os.path.exists(fn) == False:
                print(" Phonon dispersion needs to be calculated.")
            else:
                fmin, fmax = _read_frequency_range(fn, format='anphon')
                self._frequency_range = [fmin, fmax]
        return self._frequency_range
    
    @property
    def minimum_frequency(self):
        """ Read minimum frequencies from log files for band (band.log) and DOS
        (dos.log) calculations and return the value
        """
        from auto_kappa.alamode.log_parser import (
                get_minimum_frequency_from_logfile)
        
        log_band = self.out_dirs['harm']['bandos'] + "/band.log"
        log_dos = self.out_dirs['harm']['bandos'] + "/dos.log"
        
        fmin = 100000.
        try:
            out_band = get_minimum_frequency_from_logfile(log_band)
            fmin = min(fmin, out_band['minimum_frequency'])
        except Exception:
            out_band = None
        
        try:
            out_dos = get_minimum_frequency_from_logfile(log_dos)
            fmin = min(fmin, out_dos['minimum_frequency'])
        except Exception:
            out_dos = None
        
        if fmin > 1e3:
            return None
        else:
            return fmin
    
    def _get_file_pattern(self, order):
        
        if order == 1:
            filename = (self.out_dirs['harm']['suggest'] + 
                    '/%s.pattern_HARMONIC' % (self.prefix))
        elif order == 2:
            filename = (self.out_dirs['cube']['suggest'] + 
                    '/%s.pattern_ANHARM%d' % (self.prefix, order+1))
        else:
            warnings.warn(" Error: order %d is not supported." % order)
            exit()
        return filename
    
    def _get_logfile_suggest(self, order):
        
        if order == 1:
            filename = self.out_dirs['harm']['suggest'] + '/suggest.log'
        elif order == 2:
            filename = self.out_dirs['cube']['suggest'] + '/suggest.log'
        else:
            warnings.warn(" Error: order %d is not supported." % order)
            exit()
        return filename
    
    def _get_number_of_suggested_structures(self, order: None):
        """ Return the number of the suggested structures for FD method """
        
        file_pattern = self._get_file_pattern(order)
        
        ##
        lines = open(file_pattern, 'r').readlines()
        count = 0
        for ll in lines:
            data = ll.split()
            if 'Basis' in ll:
                continue
            if ':' in ll and len(data) == 2:
                count += 1
        return count

    def _get_number_of_free_fcs(self, order: None):
        """ Return the number of the FCs for FD method """
        
        file_log = self._get_logfile_suggest(order)
        
        ##
        nfc_free = {}
        lines = open(file_log, 'r').readlines()
        count = 0
        for ll in lines:
            data = ll.split()
            if len(data) == 0:
                continue
            if 'Number of free' in ll:
                if "HARMONIC" in ll:
                    nfc_free[1] = int(data[-1])
                elif "ANHARM" in ll:
                    order = int(data[3][-1]) - 1
                    nfc_free[order] = int(data[-1])
        return nfc_free
    
    def get_suggested_structures(self, order: None, nrandom=0, disp_mode='fd',
            temperature=None, number_of_displacements=1, classical=False,
            ):
        """ Return structures in which atoms are displaced with the given method,
        ``disp_mode``.
        
        order : int 
            1 (harmonic) or 2 (cubic)
        
        disp_mode : string
            "fd" for finite displacement or
            "random_normalcoordinate" for random displacement in normal 
            coordinates
         
        nrandom : int
            Number of patterns generated with the random displacement
        
        """
        if order == 1:
            structure = self.supercell
        else:
            structure = self.supercell3
            #exit()
        
        ### get displacements
        print("")
        print(" Displacement mode :", disp_mode)
        print("")
        if disp_mode == "fd":
            
            ### ver.1 with ALM library
            #all_disps = get_finite_displacements(
            #        structure,
            #        order=order, 
            #        cutoffs=self.cutoffs[:order],  ## Bohr
            #        nbody=self.nbody,
            #        mag=self.mag                   ## Ang
            #        )
            
            ### ver.2 with alm command and generated files
            file_pattern = self._get_file_pattern(order)
            
            all_disps = self._get_displacements("fd", file_pattern=file_pattern)
         
        elif disp_mode == "random_normalcoordinate":
            
            ### 
            logfile = self.out_dirs['harm']['evec'] + '/evec_commensurate.log'
            if _alamode_finished(logfile) == False:
                propt = 'evec_commensurate'
                self.write_alamode_input(propt=propt)
                self.run_alamode(propt=propt)
            
            ###
            fevec = self.out_dirs['harm']['evec'] + '/' + self.prefix + '.evec'
            all_disps = self._get_displacements(
                    disp_mode,
                    file_evec=fevec,
                    temperature=temperature, 
                    number_of_displacements=number_of_displacements,
                    classical=classical)
         
        elif disp_mode == "random":
            
            all_disps = self._get_displacements(
                    disp_mode,
                    number_of_displacements=number_of_displacements,
                    )
         
        else:
            warnings.warn(" ERROR: displacement mode %s is not supported" %
                    (disp_mode))
            exit()
        
        ### error
        if all_disps is None:
            msg = "\n Error: Failed to obtain displacement patterns.\n"
            print(msg)
            exit()
        
        ## get pristine structure
        structures = {}
        structures['prist'] = structure.copy()
        
        ## set calculator for each pattern
        for i, displacements in enumerate(all_disps):
            
            scell = structure.copy()
            
            ## displace atoms
            scell.translate(displacements)
            
            structures[i+1] = scell
        
        return structures
    
    def _get_displacements(
            self, displacement_mode: None,
            file_pattern=None, 
            #
            file_evec=None, number_of_displacements=1, temperature=500., 
            classical=False,
            ):
        """ Return displacements
        
        Args
        -----
        file_pattern : string
            **_pattern.HARMONIC, **_pattern.ANHARM3, ...
        
        """
        from auto_kappa.alamode.tools.VASP import VaspParser
        from auto_kappa.alamode.tools.GenDisplacement import AlamodeDisplace
        
        codeobj = VaspParser()
        codeobj.set_initial_structure(self.supercell)
        
        ###########################################################
        ##
        ## This part may need to be fixed.
        ##
        almdisp = AlamodeDisplace(
                displacement_mode, codeobj,
                file_evec=file_evec,
                primitive=self.primitive,
                verbosity=self.verbosity
                )
        
        if almdisp is None:
            warnings.warn(" Error: Couldn't obtain AlamodeDisplace object properly.")
            return None
        
        msg = "\n"
        msg += " Generate displacements with an Alamode tool\n"
        msg += " Displacement mode : %s\n" % displacement_mode
        msg += " %d patterns will be generated.\n" % (number_of_displacements)
        msg += "\n"
        print(msg)
        
        if displacement_mode == 'fd':
            
            if file_pattern is None:
                warnings.warn(" Error: file_pattern must be given.")
                exit()

            header_list, disp_list = almdisp.generate(
                    file_pattern=[file_pattern], magnitude=self.magnitude,
                    )
            
        elif displacement_mode == 'random_normalcoordinate':
            
            header_list, disp_list = almdisp.generate(
                    temperature=temperature,
                    number_of_displacements=number_of_displacements,
                    classical=classical
                    )

        elif displacement_mode == 'random':
            
            header_list, disp_list = almdisp.generate(
                    number_of_displacements=number_of_displacements,
                    magnitude=self.magnitude2
                    )

        else:
            msg = " Error: %s is not supported." % (displacement_mode)
            warnings.warn(msg)
            exit()
        
        all_disps = np.zeros_like(disp_list)
        for i, each in enumerate(disp_list):
            all_disps[i] = np.dot(each, self.supercell.cell)
        
        return all_disps 
    
    
    def calc_forces(self, order: None, calculator=None, 
            nmax_suggest=100, frac_nrandom=10., 
            temperature=500., classical=False,
            calculate_forces=True, output_dfset=1,
            ):
        ##nmax_suggest=100, frac_nrandom=0.02, 
        """ Calculate forces for harmonic or cubic IFCs and make a DFSET file in
        out_dirs['result'] directory.
        VASP output will be stored in self.out_dirs['harm/cube']['force'].
        When the number of suggested patterns with ALM, nsuggest, is larger than
        ``nmax_suggest``, LASSO will be used and atoms are displaced with a random 
        displacement with normal coordinate at the given temperature. This
        method is implemanted in displace.py, an Alamode tool. The number of 
        the generated patterns for LASSO is 
        max(int(``nsuggest`` * ``frac_nrandom`` + 0.5), ``nmax_suggest``).
        
        Args
        ======
        order : int
            1 (harmonic) or 2 (cubic)
        
        calculator : ase.calculators.vasp.Vasp
        
        calculate_forces : bool, [False]
            If False, structures will be created, but will be not calculated.
            If True, generated structures will be calculated.
            
        nmax_suggest : int, default 200
            If the number of suggested structures exceeds 
        
        frac_nrandom : float, default 0.02
            Ratio of the number of patterns generated with a random displacement
            (``nrandom``) to that of the suggested patterns with ALM
            (``nsuggest``). \\
            ``nrandom``:math:`= \\max(` 
                ``frac_nrandom`` :math:`\\times` ``nsuggest``,
                ``nsuggest`` :math:`)`
        
        output_dfset : int or bool, default=1
            If False or 0, DFSET file is not written.
            If 1, output DFSET file only when it doesn't exist.
            If 2, output DFSET file even when it exists.
        
        """
        if calculate_forces:
            line = "Force calculation (order: %s)" % order
        else:
            line = "Generate displacement structures (order: %s)" % order
        
        msg = "\n\n " + line + "\n"
        msg += " " + "=" * (len(line)) + "\n"
        print(msg)
        
        self._nmax_suggest = nmax_suggest
        
        ### get suggsted structures with ALM
        ### structures : dict of structure objects
        if order <= 2:
            
            nsuggest = self._get_number_of_suggested_structures(order)
            
            nfcs = self._get_number_of_free_fcs(order)[order]
            
            if order == 1:
                natoms = len(self.supercell)
            else:
                natoms = len(self.supercell3)
            
            msg = "\n Number of the suggested structures with ALM : %d\n" % (nsuggest)
            print(msg)
            
        else:
            
            nsuggest = nmax_suggest + 1
            
            nfcs = self._get_number_of_free_fcs(order)
                
            natoms = len(self.supercell3)
            
            print("")
            print(" This part is still under development.")
            print("")
            exit()
        
        ### output directory
        if order == 1:
            outdir0 = self.out_dirs['harm']['force']
        elif order == 2:
            
            if nsuggest <= nmax_suggest:
                
                self._fc3_type = 'fd'
                outdir0 = self.out_dirs['cube']['force_fd']
            
            else:

                self._fc3_type = 'lasso'
                
                ### ver.2
                outdir0 = self.out_dirs['cube']['force_lasso']
        
        elif order > 2:
            
            outdir0 = self.out_dirs['higher']['force']
        
        else:
            warnings.warning(" WARNING: given order (%d) is not supported yet." % order)
            exit()
        
        ### Compressive sensing, LASSO, is used when the number of the 
        ### suggested structure (nsuggest) exceeds nmax_suggest.
        ### When nsuggest > nmax_suggest, max(nsuggest * frac_nrandom,
        ### nmax_suggest) structures will be generated with random displacement
        ### based on the normal coordinates.
        if order != 1 and nsuggest > nmax_suggest:
        
            ###
            #self.lasso = True
            
            ### ver.1
            #nrandom = int(frac_nrandom * nsuggest + 0.5)
            #
            ### ver.2
            nrandom = int(frac_nrandom * nfcs / natoms)
             
            ngenerated = max(5, nrandom)
            
            msg = "\n"
            msg += " Maximum limit of the number of suggested patterns : %d\n" % (nmax_suggest)
            msg += "\n"
            msg += " Number of suggested patterns exceeds the maximum limit.\n"
            msg += "\n"
            msg += " Fractional number of the random patterns : %.3f\n" % (frac_nrandom)
            print(msg)
            
            if order == 2:
                ## FC3 is obtained with random-displacement method
                ## with a fixed displacement magnitude
                structures = self.get_suggested_structures(
                        order, 
                        disp_mode='random',
                        number_of_displacements=ngenerated,
                        )

            elif order > 2:

                print("")
                print("")
                print(" NEED TO BE CHECKED")
                print("")
                print("")
                exit()
                ## High order FCs are obtained with
                ## a random-displacment based on normal coordinate
                structures = self.get_suggested_structures(
                        order, 
                        disp_mode='random_normalcoordinate',
                        number_of_displacements=ngenerated,
                        temperature=temperature,
                        classical=classical
                        )
        
        else:
             
            structures = self.get_suggested_structures(order, disp_mode='fd')
        
        num_done = 0
        nsuggest =  len(structures)
        for ii, key in enumerate(list(structures.keys())):
            
            outdir = outdir0 + '/' + str(key)
            if wasfinished(outdir):
                print(" %s: skipped" % outdir)
                num_done += 1
                continue
            
            ## set output directory
            calculator.directory = outdir
            
            if ii == 0:
                print_vasp_params(calculator.asdict()['inputs'])
            
            structure = structures[key].copy()
            
            ### ver.1: with ASE
            #structure.calc = calculator
            #ene = structure.get_potential_energy()
            #print(" %s: %.3f eV" % (outdir, ene))
            #
            ### ver.2: with Custodian
            if os.path.exists(calculator.directory) == False:
                os.makedirs(calculator.directory, exist_ok=True)
                calculator.write_input(structure)
            
            if calculate_forces:
                
                run_vasp(calculator, structure, method='custodian')
                num_done += 1
            
            print(" %s" % (calculator.directory))
        
        ### output DFSET
        if num_done == len(structures) and output_dfset:
             
            if order == 1:
                
                fn0 = self.outfiles['harm_dfset']
            
            elif order == 2:
                
                fn0 = self.outfiles['cube_%s_dfset' % self.fc3_type]
            
            elif order > 2:

                fn0 = self.outfiles['lasso_dfset']

            else:
                warnings.warn(" Error: order=", order, " is not supported yet.")
                exit()
            
            os.makedirs(self.out_dirs['result'], exist_ok=True)
            offset_xml = outdir0 + '/prist/vasprun.xml'
            outfile = self.out_dirs['result'] + "/" + fn0
            if os.path.exists(outfile) == False or output_dfset == 2:
                out = get_dfset(
                        outdir0, offset_xml=offset_xml, outfile=outfile,
                        nset=nsuggest-1,
                        )
            else:
                msg = "\n"
                msg += " %s already exists\n" % (outfile)
                print(msg)
        print("")
    
    def write_alamode_input(
            self, propt=None, order=None, kpts=[15,15,15],
            outdir=None,
            **kwargs
            ):
        """ Write Anphon Input

        Args
        -----
        propt : string
            "band", "phonons", "kappa", 'evec_commensurate'
        
        kpts : list of int, shape=(3)
            k-points
        """
        ### get alamode_type and mode
        out = self._get_alamodetype_mode(propt)
        if out is None:
            print("")
            warnings.warn(" Error: %s is not supported yet." % propt)
            exit()
        alamode_type = out[0]
        mode = out[1]
        
        fc2xml = None
        fc3xml = None
        
        ### prepare filenames
        if propt == 'band' or propt == 'dos':
            
            dir_work = self.out_dirs['harm']['bandos']
            fcsxml = '../../result/' + self.outfiles['harm_xml']
            
        elif propt == 'evec_commensurate':
            
            dir_work = self.out_dirs['harm']['evec']
            fcsxml = '../../result/' + self.outfiles['harm_xml']
            
        elif propt == 'kappa':
            
            dir_work = self.out_dirs['cube']['kappa_%s' % self.fc3_type]
            fcsxml = '../../result/' + self.outfiles['cube_%s_xml' % self.fc3_type]
            
        elif propt == 'lasso' or propt == 'cv':
            
            if order == 2:
                dir_work = self.out_dirs['cube'][propt]
                dfset = ('../../result/' + 
                        self.outfiles['cube_%s_dfset' % self.fc3_type])
            
            else:
                dir_work = self.out_dirs['higher'][propt]
                dfset = '../../result/' + self.outfiles['lasso_dfset']
                
                ## get file name of FC3 (fc3xml)
                fc3xml = ('../../result/' + 
                        self.outfiles['cube_%s_xml' % self.fc3_type])
                fn_check = (
                        self.out_dirs['result'] + '/' + 
                        self.outfiles['cube_%s_xml' % self.fc3_type]
                        )
                if os.path.exists(fn_check) == False:
                    warnings.warn(" Error: FC3 has not been calculated.")
                    exit()
            
            ## get fc2xml
            fc2xml = '../../result/' + self.outfiles['harm_xml']
            fn_check = (
                    self.out_dirs['result'] + '/' + self.outfiles['harm_xml']
                    )
            if os.path.exists(fn_check) == False:
                warnings.warn(" Error: FC2 has not been calculated.")
                exit()
                
        elif propt == 'fc2':
            
            dir_work = self.out_dirs['harm']['force']
            dfset = '../../result/' + self.outfiles['harm_dfset']
            fc2xml = None
        
        elif propt == 'fc3':
            
            dir_work = self.out_dirs['cube']['force_%s' % self.fc3_type]
            dfset = ('../../result/' + 
                    self.outfiles['cube_%s_dfset' % self.fc3_type])
            fc2xml = '../../result/' + self.outfiles['harm_xml']
            
        elif propt == 'suggest':
            
            if order is None:
                print("")
                print(" Order is not given. Set order = 1.")
                order = 1
            
            dfset = None
            fc2xml = None
            
            if order == 1:
                dir_work = self.out_dirs['harm']['suggest']
            elif order == 2:
                dir_work = self.out_dirs['cube']['suggest']
            elif order > 2:
                dir_work = self.out_dirs['higher']['suggest']
            else:
                warnings.warn(" Error: order=", order, " is not supported.")
                exit()
        else:
            print("")
            warnings.warn(" Error: %s is not supported yet." % (propt))
            exit()

        ##
        born_xml = self.out_dirs['nac'] + '/vasprun.xml'
        
        ## output directory
        if outdir is not None:
            dir_work = outdir
        
        ## prepare directory
        os.makedirs(dir_work, exist_ok=True)
        
        ## non-analytical term
        if self.nac != 0 and alamode_type == 'anphon':
            
            borninfo = 'BORNINFO'
            
            ## make BORNINFO file
            outfile = dir_work + '/BORNINFO'
            write_born_info(born_xml, outfile=outfile) 
        else:
            borninfo = None
        
        ## set kpmode
        kpmode = None
        if propt == 'band':
            kpmode = 1
        elif propt == 'dos' or propt == 'kappa':
            kpmode = 2
        elif propt == 'evec_commensurate':
            kpmode = 0
        
        if alamode_type == 'anphon':
            
            ### set input file for anphon
            inp = AnphonInput.from_structure(
                    self.primitive,
                    mode=mode,
                    kpmode=kpmode,
                    fcsxml=fcsxml,
                    nonanalytic=self.nac, borninfo=borninfo
                )
        
            ### set primitive cell with Pymatgen-structure
            inp.set_primitive(
                    change_structure_format(
                        self.primitive, format="pymatgen-structure"
                        )
                    )
        
        elif alamode_type == 'alm':
            
            ### set input file for alm
            inp = AlmInput.from_structure(
                    self.supercell,
                    mode=mode,
                    dfset=dfset,
                    fc2xml=fc2xml,
                    fc3xml=fc3xml,
                    nonanalytic=self.nac, borninfo=borninfo
                )
        
        ###
        self._prefix = inp['prefix']
        
        ## set kpoints and write a file
        filename = dir_work + '/' + propt + '.in'
        if propt == 'band':
            
            inp.set_kpoint(deltak=0.01)
        
        elif propt == 'dos':
            
            if self.frequency_range is not None:
                df = self.frequency_range[1] - self.frequency_range[0]
                fmin = self.frequency_range[0] - df * 0.05
                fmax = self.frequency_range[1] + df * 0.05
            
            npoints = 301
            inp.update({'emin': fmin})
            inp.update({'emax': fmax})
            inp.update({'delta_e': (fmax-fmin)/(npoints-1)})
            inp.update({'kpts': kpts})
            inp.update({'pdos': 1})
        
        elif propt == 'evec_commensurate':
            
            ###### supercell matrix wrt primitive cell
            #### ver. old
            #### critical error!!!
            ##mat_p2s_tmp_old = np.dot(self.scell_matrix, np.linalg.inv(self.primitive_matrix))
            ##
            #### ver. corrected
            ##Uinv = np.linalg.inv(self.unitcell.cell)
            ##mat_p2s_tmp = np.dot(
            ##        np.dot(
            ##            Uinv,
            ##            np.linalg.inv(
            ##                self.primitive_matrix.T
            ##                )
            ##            ),
            ##        np.dot(
            ##            self.unitcell.cell.T,
            ##            self.scell_matrix
            ##            )
            ##        )
            ### ver. corrected 2
            mat_p2s_tmp = np.dot(
                    np.linalg.inv(self.primitive_matrix),
                    self.scell_matrix
                    )
            
            mat_p2s = np.rint(mat_p2s_tmp).astype(int)
            diff_max = np.amax(abs(mat_p2s - mat_p2s_tmp))
            if diff_max > 1e-3:
                print("")
                print(" CAUTION: please check the cell size of primitive and "\
                        "supercell")
                print(diff_max)
            
            ### commensurate points
            from auto_kappa.structure.crystal import get_commensurate_points
            comm_pts = get_commensurate_points(mat_p2s)
            inp.update({'printevec': 1})
            inp.set_kpoint(kpoints=comm_pts)
        
        elif propt == 'kappa':
            
            inp.update({'kpts': kpts})
            inp.update({'nac':  self.nac})
            inp.update({'isotope': 2})
            inp.update({'kappa_coherent': 1})
            inp.update({'tmin': 50})
            inp.update({'tmax': 1000})
            inp.update({'dt': 50})
        
        elif propt in ['cv', 'lasso', 'fc2', 'fc3', 'suggest']:
            """
            Comments
            ---------
            * maxalpha is automatically set by the code.
            * minalpha = maxalpha * 1e-6 

            * set nbody and cutoffs
            
            """
            if propt == 'fc2':
                order = 1
            elif propt == 'fc3':
                order = 2
            elif propt == 'cv' or propt == 'lasso':
                if order is None:
                    #order = self.order_lasso
                    warnings.warn(" Error: order must be given.")
                    exit()
            elif propt == 'suggest':
                if order is None:
                    warnings.warn(" ERROR: order must be given.")
                    exit()
            else:
                print(" Error")
                exit()
            
            if len(self.nbody) < order:
                self.set_nbody_automatically()
            inp.update({'norder': order})
            inp.update({'nbody': [self.nbody[i] for i in range(order)]})
            
            ### set cutoffs for alamode
            cutoffs_alamode = {}
            for i1 in range(self.num_elems):
                for i2 in range(self.num_elems):
                    lab = "%s-%s" % (
                            inp.as_dict()['kd'][i1],
                            inp.as_dict()['kd'][i2]
                            )
                    cutoffs_alamode[lab] = np.where(
                            self.cutoffs[:order,i1,i2]<0., 
                            None, self.cutoffs[:order,i1,i2]
                            )
            
            inp.update({'cutoff': cutoffs_alamode})
            
            if propt in ['cv', 'lasso']:
                inp.update({'lmodel': 'enet'})
                inp.update({'l1_ratio': 1.0})
                inp.update({'maxiter': 10000})   ## use the original default, 10000
                inp.update({'conv_tol': 1e-10})  ## strincter than the default, 1e-8
                if propt == 'cv':
                    inp.update({'cv': 5})
                    ### Alamode package set automatically.
                    #inp.update({'cv_maxalpha': 1e-2})
                    #inp.update({'cv_minalpha': 1e-8})
                    inp.update({'cv_nalpha': 50})
                elif propt == 'lasso':
                    
                    alpha = self.get_suggested_l1alpha(order=order)
                    inp.update({'cv': 0})
                    inp.update({'l1_alpha': alpha})      ### read l1_alpha
        
        else:
            warnings.warn(" Error: %s is not supported." % propt)
            exit()
        
        inp.update(kwargs)
        inp.to_file(filename=filename)
    
    def get_suggested_l1alpha(self, order=None):
        
        if order == 2:
            fn = self.out_dirs['cube']['cv']+'/'+self.prefix+'.cvscore'
        else:
            fn = self.out_dirs['higher']['cv']+'/'+self.prefix+'.cvscore'
        
        try:
            lines = open(fn, 'r').readlines()
            for ll in lines:
                if "Minimum CVSCORE" in ll:
                    alpha = float(ll.split()[-1])
                    return alpha
            return None
        except Exception:
            warnings.warn(" Warning: cannot find %s" % fn)
            return None
    
    def _get_alamodetype_mode(self, propt):
        
        ### get alamode_type and mode 
        flag = False
        for k1 in self.propts_mode_type:
            for k2 in self.propts_mode_type[k1]:
                if propt in self.propts_mode_type[k1][k2]:
                    return k1, k2
        if flag == False:
            return None
     
    def run_alamode(self, propt=None, order=None, neglect_log=False, 
            outdir=None, **args):
        """ Run anphon
        
        Args
        ---------
        propt : string
            "band", "dos", "evec_commensurate", or "kappa"

        neglect_log : bool
            If False, anphon will not be run if the same calculation had been
            conducted while, if True, anphon will be run forecely.
        
        """
        ### get alamode_type and mode
        out = self._get_alamodetype_mode(propt)
        if out is None:
            warnings.warn(" Error: %s is not supported yet." % propt)
            exit()
        alamode_type = out[0]
        mode = out[1]

        ## change directory
        if propt == 'band' or propt == 'dos':
            
            workdir = self.out_dirs['harm']['bandos']

            ### store output file names
            if propt == 'band':
                ext = 'bands'
            else:
                ext = propt
        
        elif propt == 'kappa':
            
            #if newversion == False:
            #    if self.lasso == False:
            #        workdir = self.out_dirs['cube']['kappa_fd']
            #    else:
            #        workdir = self.out_dirs['higher']['kappa']
            #else:
            ##    
            workdir = self.out_dirs['cube']['kappa_%s' % self.fc3_type]
        
        elif propt == 'evec_commensurate':
            
            workdir = self.out_dirs['harm']['evec']
        
        elif propt == 'cv' or propt == 'lasso':
            
            #if newversion == False:
            #    workdir = self.out_dirs['higher'][propt]
            #else:
            ##
            if order == 2:
                workdir = self.out_dirs['cube'][propt]
            else:
                workdir = self.out_dirs['higher'][propt]
            
        elif propt == 'fc2':

            workdir = self.out_dirs['harm']['force']
        
        elif propt == 'fc3':

            workdir = self.out_dirs['cube']['force_%s' % self.fc3_type]
        
        elif propt == 'suggest':
            
            if order is None:
                order = 1
            
            if order == 1:
                workdir = self.out_dirs['harm']['suggest']
            elif order == 2:
                workdir = self.out_dirs['cube']['suggest']
            elif order > 2:
                workdir = self.out_dirs['higher']['suggest']
            else:
                print("")
                warnings.warn(" Error: order must be gien properly.")
                exit()

        else:
            print("")
            warnings.warn(" WARNNING: %s property is not supported." % (propt))
            exit()
        
        ### modify work directory
        if outdir is not None:
            workdir = outdir
        
        ### print title
        msg = "\n"
        line = "Run %s for %s:" % (alamode_type, propt)
        msg += " " + line + "\n"
        msg += " " + "-" * (len(line)) + "\n\n"
        msg += " Working directory : " + workdir + "\n"
        print(msg)
            
        filename = "%s.in" % propt
        logfile = propt + '.log'
        
        ## prepare command and environment
        if alamode_type == 'anphon' and self.commands['alamode']['anphon_para'] == 'mpi':
            _nprocs = self.commands['alamode']['ncores']
            _nthreads = 1
        else:
            _nprocs = 1
            _nthreads = self.commands['alamode']['ncores']
            
        val = run_alamode(filename, logfile, workdir=workdir, 
                neglect_log=neglect_log,
                mpirun=self.commands['alamode']['mpirun'], 
                nprocs=_nprocs,
                nthreads=_nthreads,
                command=self.commands['alamode'][alamode_type],
                )
        
        if val == 2:
            print(" %s has been already calculated." % propt)
        
        if mode == 'optimize' and propt in ['lasso', 'fc2', 'fc3']:
            
            ### copy the generated FCs file to "result" directory
            self._copy_generated_fcsfiles(
                    propt=propt, order=order,
                    )
            
            if propt in ['lasso']:
                self._plot_cvsets(order=order)

    def _copy_generated_fcsfiles(self, propt=None, order=None):
        """ Copyr a FCs file into the "result" directory.
        """
        if propt == 'lasso':
            
            #if newversion == False:
            #    fn1 = self.out_dirs['higher']['lasso']+'/'+self.prefix+'.xml'
            #    fn2 = self.out_dirs['result']+'/'+self.outfiles['lasso_xml']
            #else:
            ##
            if order == 2:
                fn1 = self.out_dirs['cube']['lasso']+'/'+self.prefix+'.xml'
                fn2 = self.out_dirs['result']+'/'+self.outfiles['cube_lasso_xml']
            else:
                fn1 = self.out_dirs['higher']['lasso']+'/'+self.prefix+'.xml'
                fn2 = self.out_dirs['result']+'/'+self.outfiles['lasso_xml']
        
        elif propt == 'fc2':
            fn1 = self.out_dirs['harm']['force']+'/'+self.prefix+'.xml'
            fn2 = self.out_dirs['result']+'/'+self.outfiles['harm_xml']
         
        elif propt == 'fc3':
            fn1 = (self.out_dirs['cube']['force_%s' % self.fc3_type] + 
                    '/' + self.prefix+'.xml')
            
            #if newversion == False:
            #    fn2 = self.out_dirs['result']+'/'+self.outfiles['cube_xml']
            #else:
            ##
            fn2 = (
                    self.out_dirs['result'] + '/' + 
                    self.outfiles['cube_%s_xml' % self.fc3_type]
                    )
        
        else:
            print("")
            warnings.warn(" Error: %s is not supported yet." % propt)
            exit()
        
        ##
        if os.path.exists(fn1) == False:
            
            warnings.warn(' %s does not exist.' % fn1)
        
        else:
            if os.path.exists(fn2):
                msg = "\n"
                msg += " %s was overwritten.\n" % (fn2)
            else:
                msg = "\n"
                msg += " %s was created.\n" % fn2
            print(msg)
            shutil.copy(fn1, fn2)
    
    #def write_lifetime(self, temperature=300, outfile=None):
    #    
    #    from auto_kappa.alamode.tools.analyze_phonons import (
    #            write_lifetime_at_given_temperature)
    #    
    #    msg = "\n"
    #    msg += " ### Analyze lifetime at %.1f K" % (temperature)
    #    print(msg)
    #    
    #    ### output file name
    #    if self.lasso:
    #        dir_kappa = self.out_dirs['lasso']['kappa']
    #    else:
    #        dir_kappa = self.out_dirs['cube']['kappa']
    #    
    #    file_result = dir_kappa + '/%s.result' % (self.prefix)
    #    file_isotope = dir_kappa + '/%s.self_isotope' % (self.prefix)
    #    
    #    if outfile is None:
    #        outfile = dir_kappa + '/tau_%dK.dat' % (int(temperature))
    #    
    #    ### analyze lifetime
    #    write_lifetime_at_given_temperature(
    #            file_result, file_isotope=file_isotope, 
    #            temperature=temperature, outfile=outfile)
    #    
    #    print("")
    #    print(" Output", outfile)
    #    print("")
    

    @property
    def file_result(self):
        return self._file_result

    @property
    def file_isotope(self):
        return self._file_isotope

    @property
    def scat(self):
        return self._scat
    
    def set_scattering_info(self, grain_size=None, temperature=300):
        
        #### output file name
        #if self.lasso:
        #    dir_kappa = self.out_dirs['lasso']['kappa']
        #else:
        #    dir_kappa = self.out_dirs['cube']['kappa_%s' % self.fc3_type]
        
        dirs_kappa = self.get_kappa_directories()
        
        ## get directory name for kappa with the finest k grid
        dir_kappa = None
        lab_kpts = None
        ksum = 0
        for lab in dirs_kappa:
            data = lab.lower().split("x")
            if len(data) != 3:
                continue
            kpts = np.asarray([int(d) for d in data])
            if np.sum(kpts) > ksum:
                dir_kappa = dirs_kappa[lab]
                lab_kpts = lab
                ksum = np.sum(kpts)
        
        ### file check
        self._file_result = dir_kappa + '/%s.result' % (self.prefix)
        self._file_isotope = dir_kappa + '/%s.self_isotope' % (self.prefix)
        
        if os.path.exists(self.file_result) == False:
            warnings.warn(" Error: %s does not exist." % self.file_result)
        
        if os.path.exists(self.file_isotope) == False:
            warnings.warn(" Error: %s does not exist." % self.file_isotope)
            self._file_isotope = None
        
        ###
        from .analyzer.scattering import Scattering
        self._scat = Scattering(
                self.file_result, 
                file_isotope=self.file_isotope,
                grain_size=grain_size, 
                temperature=temperature
                )
        
    def plot_lifetime(self, temperatures="300:500", grain_size=None):
        
        data = temperatures.split(":")
        ts = [float(t) for t in data]
        
        dfs = {}
        for t in ts:
            self.set_scattering_info(temperature=t, grain_size=None)
            omegas = self.scat.result['frequencies']
            taus = self.scat.lifetime
            
            dfs[int(t)] = pd.DataFrame()

            nk = len(omegas)
            nbands = len(omegas[0])
            dfs[int(t)]['frequency'] = omegas.reshape(nk*nbands)
            dfs[int(t)]['lifetime'] = taus.reshape(nk*nbands)
        
        figname = self.out_dirs['result'] + '/fig_lifetime.png'

        ### plot a figure
        from auto_kappa.plot.pltalm import plot_lifetime
        plot_lifetime(dfs, figname=figname)
    
    
    def plot_scattering_rates(self, temperature=300., grain_size=1000.):
        
        if self.scat is None:
            self.set_scattering_info(
                    temperature=temperature,
                    grain_size=grain_size)
        
        ## set temperature 
        if abs(self.scat.temperature - temperature) > 0.1:
            self.scat.change_temperature(temperature)
        
        ## set grain size
        if self.scat.size is None:
            self.scat.change_grain_size(grain_size)
        else:
            if abs(self.size - grain_size) > 1.:
                self.scat.cahnge_grain_size(grain_size)
        
        ## get frequencies
        frequencies = self.scat.result['frequencies']
        n1 = len(frequencies)
        n2 = len(frequencies[0])
        frequencies = frequencies.reshape(n1*n2)
        
        ## get scattering rates and set labels
        labels = {}
        scat_rates = {}
        for key in self.scat.scattering_rates:
            
            scat_rates[key] = self.scat.scattering_rates[key].reshape(n1*n2)
            
            if key == 'phph':
                labels[key] = '3ph (%dK)' % int(temperature)
            elif key == 'isotope':
                labels[key] = 'isotope'
            elif key == 'boundary':
                labels[key] = 'L=%dnm' % grain_size
            else:
                labels[key] = key
        
        ## plot a figure
        figname = self.out_dirs['result'] + '/fig_scat_rates.png'
        from auto_kappa.plot.pltalm import plot_scattering_rates
        plot_scattering_rates(frequencies, scat_rates, labels, figname=figname)
        
    
    def _plot_cvsets(self, order=None):
        
        from auto_kappa.plot.pltalm import plot_cvsets
        print("")
        print(" ### Plot CV results ###")
        
        #if newversion == False:
        #    figname = self.out_dirs['result'] + '/fig_cvsets.png'
        #    plot_cvsets(
        #            directory=self.out_dirs['lasso']['cv'], 
        #            figname=figname
        #            )
        #else:
        ##
        if order == 2:
            figname = self.out_dirs['result'] + '/fig_cvsets_cube.png'
            plot_cvsets(
                    directory=self.out_dirs['cube']['cv'],
                    figname=figname)
        else:
            figname = self.out_dirs['result'] + '/fig_cvsets.png'
            plot_cvsets(
                    directory=self.out_dirs['higher']['cv'], 
                    figname=figname
                    )
        
        print("")

    def plot_bandos(self, **args):
        
        ### set figure name
        if 'figname' not in args.keys():
            figname = self.out_dirs['result'] + '/fig_bandos.png'
        else:
            figname = args['figname']
        
        from auto_kappa.plot.bandos import plot_bandos
        
        ### output title
        line = "Plot band and DOS:"
        msg = "\n " + line + "\n"
        msg += " " + "-" * len(line) + "\n"
        print(msg)
         
        fig = plot_bandos(
                directory=self.out_dirs['harm']['bandos'], 
                prefix=self.prefix, 
                figname=figname,
                **args
                )
    
    def get_kappa_directories(self):
        
        ll_tmp = self.out_dirs['cube']["kappa_%s" % self.fc3_type] + "_"
        line = self.out_dirs['cube']["kappa_%s" % self.fc3_type] + "_*"
        dirs = glob.glob(line)
        dirs_done = {}
        for dd in dirs:
            fn_log = dd + "/kappa.log"
            if os.path.exists(fn_log) == False:
                continue
            if _alamode_finished(fn_log):
                label = dd.split(ll_tmp)[1]
                dirs_done[label] = dd
        return dirs_done

    def plot_kappa(self):
        
        #if self.lasso == False:
        #    dir_kappa = self.out_dirs['cube']["kappa_%s" % self.fc3_type]
        #else:
        #    dir_kappa = self.out_dirs['lasso']["kappa"]
        
        
        dirs_kappa = self.get_kappa_directories()
        dfs = {}
        for lab in dirs_kappa:
            dfs[lab] = _read_kappa(dirs_kappa[lab], self.prefix)
        
        from auto_kappa.plot.pltalm import plot_kappa
        figname = self.out_dirs['result'] + '/fig_kappa.png'
        plot_kappa(dfs, figname=figname)
    
    
    def plot_cumulative_kappa(self, temperatures="100:300:500", wrt='frequency',
            figname=None, xscale='linear', nbins=150):
        
        ## set grain size
        if self.scat.size is not None:
            self.scat.change_grain_size(None)
        
        ## set temperatures
        data = temperatures.split(":")
        ts = [float(t) for t in data]
        dfs = {}
        for t in ts:
            self.scat.change_temperature(t)
            dfs[int(t)] = self.scat.get_cumulative_kappa(
                    temperature=t, wrt=wrt, xscale=xscale, nbins=nbins)
        
        ##
        lab_kappa = "${\\rm \\kappa_{lat}}$"
        if 'freq' in wrt:
            xlabel = "Frequency (${\\rm cm^{-1}}$)"
            unit1 = "${\\rm Wm^{-1}K^{-1}/cm^{-1}}$"
            ylabel1 = "Spectral %s (%s)" % (lab_kappa, unit1)
        else:
            xlabel = "Mean free path (nm)"
            unit1 = "${\\rm Wm^{-1}K^{-1}/nm}$"
            ylabel1 = "Spectral %s (%s)" % (lab_kappa, unit1)
        
        if figname is None:
            figname = self.out_dirs['result'] + '/fig_cumu_%s.png' % wrt
            
        from auto_kappa.plot.pltalm import plot_cumulative_kappa
        plot_cumulative_kappa(dfs, xlabel=xlabel, figname=figname, 
                xscale=xscale, ylabel=ylabel1)


def run_alamode(filename, logfile, workdir='.', neglect_log=False,
        mpirun='mpirun', nprocs=1, nthreads=1, command='anphon'):
    """ Run alamode with a command (alm or anphon)
    
    Args
    ======
    filename : string
        input script of Alamode in workdir
    
    logfile : string
        log file name in workdir
    
    workdir : string
        work directory
    
    """
    cmd = "%s -n %d %s %s" %(mpirun, nprocs, command, filename)
    
    ### change directory
    dir_init = os.getcwd()
    os.chdir(workdir)
    
    ## If the job has been finished, the same calculation is not conducted.
    ## The job status is judged from *.log file.
    if _alamode_finished(logfile) == False or neglect_log:
        
        os.environ['OMP_NUM_THREADS'] = str(nthreads)
        
        ## run the job!!
        with open(logfile, 'w') as f:
            
            ### ver.1
            #proc = subprocess.Popen(
            #        cmd.split(), env=os.environ, stdout=f,
            #        stderr=subprocess.PIPE)
            #
            ### ver.2
            proc = subprocess.Popen(
                    cmd, shell=True, env=os.environ, stdout=f,
                    stderr=subprocess.PIPE)
            
            proc.wait()
        
        val = 1

    else:
        val = 2
    
    #### Return to the original directory
    ## sometimes this way leads to a problem.
    #os.environ.pop('OMP_NUM_THREADS', None)
    os.environ["OMP_NUM_THREADS"] = "1"
    os.chdir(dir_init)
    
    return val

def get_cutoffs_automatically(cutoff2=-1, cutoff3=4.3, num_elems=None, order=5):
    """ 
    Args
    -------
    #cutoffs : shape=(order, num_elems, num_elems), unit=[Ang]
    #    If cutoffs is given properly, cutoffs is used and cutoff2 and
    #    cutoff3 are neglected.
    #    If cutoffs is not given, cutoffs are given automatically with
    #    cutoff2 and cutoff3.
    cutoff2 : float, unit=Ang
    cutoff3 : float, unit=Ang
    """
    cutoffs = []
    n = num_elems
    for i in range(order):
        if i == 0:
            cc = cutoff2
        else:
            cc = cutoff3
        cutoffs.append(np.asarray(np.ones((n,n)) * cc))
    return np.asarray(cutoffs)

def _read_frequency_range(filename, format='anphon'):
    """ read minimum and maximum frequencies from .bands file created by anphon 
    """
    if format == 'anphon':
        data = np.genfromtxt(filename)
        fmax = np.max(data[:,-1])
        fmin = np.min(data[:,1])
        return fmin, fmax
    else:
        warnings.warn(" Error: %s is not supported yet." % format)
        exit()


def _alamode_finished(logfile):
    try:
        lines = open(logfile, 'r').readlines()
        n = len(lines)
        for i in range(10):
            line = lines[n-1-i]
            data = line.split()
            if len(data) != 0:
                if 'Job finished' in line:
                    return True
                #else:
                #    return False
    except Exception:
        return False
    
    return False

def run_alm(structure, order, cutoffs, nbody, mode=None,
        displacements=None, forces=None, outfile=None,
        fc2info=None, lasso=False, lasso_type='alm', verbosity=0
        ):
    """ get ALM object
    Note : length unit is Bohr

    Args
    ==========
    order : int
    cutoff : float, unit=Bohr
    ndoby : shape=(order)
    
    displacements : shape=(npatterns, natoms, 3)
    forces : shape=(npatterns, natoms, 3)
    outfile : string
    
    lasso : bool

    lasso_type : string
        "alm" or "scikit"
    
    Return
    ===========
    If mode == 'suggest':
    patterns : shape=(nstruct, natoms, 3)

    If mode == 'optimize':
    alm.get_fc(order), which are fc*_values, elem*_indices
    
    """
    if outfile is not None:
        if os.path.exists(outfile):
            msg = "\n"
            msg += " ALM calculation (%s) has been already finished.\n" % (mode)
            msg += " See: %s" % (outfile)
            msg += "\n"
            print(msg)
            return None
    
    from alm import ALM

    if type(structure) == "ase.atoms.Atoms":
        atoms = structure.copy()
    else:
        atoms = change_structure_format(structure, format='ase')
    
    lave = atoms.cell * AToBohr
    xcoord = atoms.get_scaled_positions()
    kd = atoms.get_atomic_numbers()
    
    #if lasso and lasso_type == 'scikit':
    #    msg = "\n"
    #    msg = " Perform LASSO with Scikit_learn...\n"
    #    msg = "\n"
    #    print(msg)
    #    from .alamode_tools.lasso import (
    #            get_training_data,
    #            run_lasso_by_scikit_learn
    #            )
    #    from . import default_lassobyscikit_parameters
    #    X, y = get_training_data(
    #            [lave, xcoord, kd], 
    #            displacements, forces,
    #            maxorder=order, cutoff=cutoffs, nbody=nbody
    #            )
    #    ###
    #    ### Default parameters may need to be adjusted.
    #    ###
    #    fc, alphas_lasso, rmse_mean, cv_mean = run_lasso_by_scikit_learn(
    #            X, y, len(atoms), len(forces), forces.ravel(),
    #            **default_lassobyscikit_parameters,
    #            )
        
    ###
    with ALM(lave, xcoord, kd) as alm:
        
        alm.define(order, cutoff_radii=cutoffs, nbody=nbody)
        alm.set_verbosity(verbosity)
        
        if mode == 'suggest':
            info = alm.suggest()
            if info == 1:
                warnings.warn(" Error during ALM calculation")
                return None
            else:
                patterns = alm.get_displacement_patterns(order)
                return patterns
            
        elif mode == 'optimize':
            
            alm.displacements = displacements
            alm.forces = forces
            
            ### for lasso
            #if lasso:
            #    alm.set_constraint(translation=True)
            
            ## freeze fc2
            if fc2info is not None:
                alm.freeze_fc(fc2info[0], fc2info[1])
            
            ## lasso
            if lasso:
                
                ### cross-validation
                optcontrol = {'linear_model': 2,
                              'cross_validation': 5,
                              'num_l1_alpha': 50}
                alm.set_optimizer_control(optcontrol)
                alm.optimize()
                
                ### prepare for optimization with lasso
                optcontrol['cross_validation'] = 0      ## change N -> 0
                optcontrol['l1_alpha'] = alm.get_cv_l1_alpha()
                alm.set_optimizer_control(optcontrol)
                
            ### optimization!!
            info = alm.optimize()
            
            if outfile is not None:
                alm.save_fc(outfile, format='alamode')
                print(" Output", outfile)
            return alm.get_fc(order)
        
        else:
            warnings.warn(" WARNING: mode %s is nto supported" % (mode))
    

def _read_kappa(dir_kappa, prefix):
    """ Read .kl and .kl_coherent files and return pandas.DataFrame object
    
    Args
    ------
    dir_kappa : string
        directory in which .kl and .kl_coherent should exist.
    
    prefix : string
    
    """
    fn_kp = dir_kappa + '/' + prefix + '.kl'
    fn_kc = dir_kappa + '/' + prefix + '.kl_coherent'
    
    ##
    if os.path.exists(fn_kp) == False:
        return None
    
    df = pd.DataFrame()
    data = np.genfromtxt(fn_kp)
    
    if os.path.exists(fn_kc) == False:
        fn_kc = None
        data2 = None
    else:
        data2 = np.genfromtxt(fn_kc)

        if np.max(abs(data[:,0] - data2[:,0])) > 0.5:
            warnings.warn(" Warning: temperatures are incompatible.")
    
    nt = len(data)
    df['temperature'] = data[:,0]
    dirs = ['x', 'y', 'z']
    for i1 in range(3):
        d1 = dirs[i1]
        
        if data2 is not None:
            dd = dirs[i1]
            lab2 = 'kc_%s%s' % (dd, dd)
            df[lab2] = data2[:,i1+1]
        
        for i2 in range(3):
            d2 = dirs[i2]
            num = i1*3 + i2 + 1
            lab = 'kp_%s%s' % (d1, d2)
            df[lab] = data[:,num]
    
    kave = (df['kp_xx'].values + df['kp_yy'].values + df['kp_zz'].values) / 3.
    df['kp_ave'] = kave
    
    if data2 is not None:
        kave = (df['kc_xx'].values + df['kc_yy'].values + df['kc_zz'].values) / 3.
        df['kc_ave'] = kave
        
        for pre in ['xx', 'yy', 'zz', 'ave']:
            key = 'ksum_%s' % pre
            df[key] = df['kp_%s'%pre].values + df['kc_%s'%pre].values
    
    return df


#def get_finite_displacements(structure, order, cutoffs=None, nbody=None, magnitude=None):
#    """ Generate displacement patterns with ALM for FCs calculation
#    Args
#    ======
#    cutoff : unit=[Bohr]
#    magnitude : unit=[Ang]
#    """
#    nbody_new = []
#    for i in range(order):
#        nbody_new.append(nbody[i])
#    
#    ## get patterns with ALM
#    patterns = run_alm(
#            structure, order, cutoffs[:order], nbody_new,
#            mode='suggest')
#    
#    if patterns is None:
#        return None
#    
#    natoms = len(structure)
#    
#    ### convert patterns to vectors
#    all_disps = []
#    for pattern in patterns:
#        all_disps.append(np.zeros((natoms, 3)))
#        for each in pattern:
#            if each[2].lower() == 'cartesian':
#                all_disps[-1][each[0]] = each[1] * magnitude
#            else:
#                print(" Error during getting patterns.")
#                exit()
#    
#    return all_disps

#def _lasso_cv(alm):
#    optcontrol = {'linear_model': 2,
#                  'cross_validation': 5,
#                  'num_l1_alpha': 50}
#    alm.set_optimizer_control(optcontrol)
#    alm.optimize()
#    return alm.get_cv_l1_alpha()
#
#def _lasso_optimize(alm, cv_l1_alpha):
#    optcontrol = {'linear_model': 2,
#                  'cross_validation': 0,  # change 2 -> 0
#                  'l1_alpha': cv_l1_alpha}
#    alm.set_optimizer_control(optcontrol)
#    alm.optimize()

