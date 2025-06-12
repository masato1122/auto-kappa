# -*- coding: utf-8 -*-
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
import sys
import os
import os.path
import numpy as np

import ase
import shutil
import pandas as pd
import glob
import itertools

import ase.io
from phonopy.structure.cells import get_primitive, get_supercell

from auto_kappa import output_directories, output_files, default_amin_parameters
from auto_kappa.structure.crystal import change_structure_format, get_formula
from auto_kappa.alamode.runjob import run_alamode
from auto_kappa.io.vasp import write_born_info
from auto_kappa.io.alm import AnphonInput
from auto_kappa.io.files import write_output_yaml

from auto_kappa.alamode.log_parser import get_minimum_frequency_from_logfile
from auto_kappa.alamode.io import wasfinished_alamode
from auto_kappa.alamode.compat import (
    check_directory_name_for_pristine,
    was_primitive_changed,
    was_tolerance_changed,
    backup_previous_results,
    adjust_keys_of_suggested_structures
)
import auto_kappa.alamode.helpers as helper
from auto_kappa.alamode.helpers import AlamodeForceCalculator, AlamodeInputWriter, NameHandler
from auto_kappa.utils import get_version
from auto_kappa.units import BohrToA, AToBohr

import logging
logger = logging.getLogger(__name__)


class AlamodeCalc(AlamodeForceCalculator, AlamodeInputWriter, NameHandler):
    
    ### k1: alamode type, k2: mode, k3: propt
    propts_mode_type = {
            'anphon': {
                'phonons': ['band', 'dos', 'evec_commensurate'],
                'rta':  ['kappa', 'kappa_scph'],
                'scph': ['scph'],
                },
            'anphon_ver2': {
                'rta': ['kappa_4ph', 'kappa_scph_4ph'],
            },
            'alm':{
                'suggest': ['suggest'],
                'optimize': ['fc2', 'fc3', 'cv', 'lasso'],
                }
            }
    
    def __init__(self, 
            prim_given, 
            base_directory='apdb_output',
            additional_directory=None,
            primitive_matrix=None,
            scell_matrix=None, 
            #scell_matrix3=None,
            restart=1,
            cutoffs=None, nbody=[2,3,3,2,2], 
            magnitude=0.01, magnitude2=0.03, 
            ##mag_high=0.03,
            cutoff2=-1, cutoff3=4.3, 
            #order_lasso=5,
            nac=None,
            commands=None,
            verbosity=0,
            yamlfile_for_outdir=None,
            ):
        """ This class helps to manage ALAMODE calculations including force
        calculations with VASP.
        
        Args
        ------

        prim_given : primitive structure
            Different formats such as pymatgen and ASE are accepted while the
            format is changed to ase-Atoms in this module.

        base_directory : string
            The top directory of the output directories. For example, material
            ID of Materials Project or composition [apdb_output]
        
        additional_directory : string
            ``base_directory`` + ``additional_directory`` will be used as the
            base directory.
        
        restart : int, default 1
            The calculation will restart (1) or NOT restart (0) when the
            directory exists.

        scell_matrix : array of float, shape=(3,3)
            transformation matrix from the unitcell to a supercell
        
        primitive_matrix : array of float, shape=(3,3)
            transformation matrix from the unitcell to the primitive cell
        
        #scell_matrix3 : array of float, shape=(3,3)
        #    transformation matrix from the unitcell to a supercell used to
        #    calculate cubic FCs

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
                "nprocs": 1, "anphon": "anphon", "alm": "alm",
                }
        
        yamlfile_for_outdir : string
            yaml file name to write output directory names
        
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
        >>>     commands={'mpirun':'mpirun', 'anphon_para':'omp', 'nprocs':1,'anphon':'anphon'},
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
        super().__init__()
        
        ### set name of the output directory
        self.set_basedir_name(base_directory, restart)
        
        ###
        self._commands = {
                'vasp':{
                    'mpirun': 'mpirun', "nprocs": 1, "nthreads": 1, 
                    'vasp': 'vasp'
                    },
                'alamode':{
                    'mpirun': 'mpirun', 'anphon_para': 1, 'nprocs': 2,
                    'alm': 'alm', 'anphon': 'anphon', 'anphon_ver2': 'anphon.2.0'
                    }
                }
        self._commands.update(commands)
        
        ### output directories
        self._additional_directory = additional_directory
        
        self.out_dirs = {}
        for k1 in output_directories.keys():
            
            base_dir = None
            if additional_directory is not None:
                if k1 not in ["relax", "nac"]:
                    base_dir = self.base_directory + "/" + additional_directory
            if base_dir is None:
                base_dir = self.base_directory
            
            ###
            values1 = output_directories[k1]
            if type(values1) == str:
                self.out_dirs[k1] = base_dir + '/' + values1
            else:
                self.out_dirs[k1] = {}
                for k2 in values1.keys():
                    values2 = values1[k2]
                    self.out_dirs[k1][k2] = base_dir + '/' + values2
        
        ###
        self.outfiles = output_files
        
        ### matrices
        if primitive_matrix is None:
            logger.warning("\n Error: primitive_matrix must be given.")
        
        if scell_matrix is None:
            logger.warning("\n Error: scell_matrix must be given.")
        
        #if scell_matrix3 is None:
        #    scell_matrix3 = scell_matrix.copy()
        
        self._primitive_matrix = primitive_matrix
        self._scell_matrix = scell_matrix
        #self._scell_matrix3 = scell_matrix3
        
        ### make the unitcell with ``get_supercell`` in Phonopy
        mat_p2u = np.linalg.inv(primitive_matrix)
        mat_p2u = np.array(np.sign(mat_p2u) * 0.5 + mat_p2u, dtype="intc")
        unit_pp = get_supercell(
                change_structure_format(prim_given, format='phonopy'), 
                mat_p2u
                )
        
        ### set unitcell, primitive cell, and two kinds of supercells
        self._primitive = None
        self._unitcell = None
        self._supercell = None
        #self._supercell3 = None
        
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
        nac_prev = self._get_previous_nac()
        if nac_prev is not None:
            if self.nac != nac_prev:
                msg = "\n NONANALYTIC option is modified from %d to %d, " % (
                        self.nac, nac_prev)
                msg += "which was used for the previous calculation."
                logger.info(msg)
            self.nac = nac_prev
        
        self.lasso = False
        self._fc3_type = ""            ### 'df', 'lasso', or ''
        self._nmax_suggest = None
        self._commands = commands
        self._scat = None
        
        ###
        self._prefix = get_formula(self.primitive)  ## prefix for input files
        self._frequency_range = None
        #self._minimum_frequency = None
        
        self._verbosity = verbosity
        self._yamlfile_for_outdir = yamlfile_for_outdir
    
        self._fc2xml = None
        self._fc3xml = None
        self._fcsxml = None
        self._higher_fcsxml = None
        self._harm_dfset = None
        self._cube_dfset = None
        self._higher_dfset = None
    
    def _get_previous_nac(self):
        """ Get previously used NAC option """
        nac_propt = {}
        count = 0
        for propt in ["band"]:
            ### alamode input file
            filename = self.out_dirs["harm"]["bandos"] + "/%s.in" % propt
            if os.path.exists(filename) == False:
                continue
            ### get NAC option in the input file
            lines = open(filename, 'r').readlines()
            for line in lines:
                if "nonanalytic" in line.lower():
                    data = line.split()
                    if len(data) >= 3:
                        nac_propt[propt] = int(data[2])
                        count += 1
        ###
        if len(nac_propt) != 0:
            return nac_propt["band"]
        else:
            return None
    
    def get_relative_path(self, abs_path):
        """ Get relative path 
        
        Args
        ------
        abs_path : string
            absolute path
        """
        try:
            dir1 = "./" + self.base_directory.split("/")[-1]
            relative_path = abs_path.replace(self.base_directory, dir1)
            return relative_path
        except Exception:
            return abs_path

    @property
    def fc2xml(self):
        if self._fc2xml is None:
            self._fc2xml = self.out_dirs['harm']['force'] + f"/{self.prefix}.xml"
        if os.path.exists(self._fc2xml) == False:
            self._fc2xml = None
        return self._fc2xml
    
    @property
    def fc3xml(self):
        if self._fc3xml is None:
            if self.fc3_type == 'fd':
                self._fc3xml = self.out_dirs['cube'][f'force_fd'] + f"/{self.prefix}.xml"
            elif self.fc3_type == 'lasso':
                self._fc3xml = self.out_dirs['cube']['lasso'] + f"/{self.prefix}.xml"
            else:
                raise ValueError(
                    f"\n Error: fc3_type must be 'fd' or 'lasso' but not '{self.fc3_type}'. "
                    "If you want to use the default, set fc3_type='df'.")
        if os.path.exists(self._fc3xml) == False:
            self._fc3xml = None
        return self._fc3xml
    
    @property
    def higher_fcsxml(self):
        if self._higher_fcsxml is None:
            self._higher_fcsxml = self.out_dirs['higher']['lasso'] + f"/{self.prefix}.xml"
        if os.path.exists(self._higher_fcsxml) == False:
            self._higher_fcsxml = None
        return self._higher_fcsxml
    
    @property
    def fcsxml(self):
        # if self._fcsxml is None:
        #     self._fcsxml = self.out_dirs["result"] + "/FCS.xml"
        if os.path.exists(self._fcsxml) == False:
            self._fcsxml = None
        return self._fcsxml
    
    @property
    def harm_dfset(self):
        if self._harm_dfset is None:
            self._harm_dfset = self.out_dirs['result'] + '/' + self.outfiles['harm_dfset']
        return self._harm_dfset
    @property
    def cube_dfset(self):
        if self._cube_dfset is None:
            self._cube_dfset = self.out_dirs['result'] + '/' + self.outfiles[f'cube_{self.fc3_type}_dfset']
        return self._cube_dfset
    @property
    def higher_dfset(self):
        if self._higher_dfset is None:
            self._higher_dfset = self.out_dirs['result'] + '/' + self.outfiles['higher_dfset']
        return self._higher_dfset
    
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
    
    @property
    def additional_directory(self):
        return self._additional_directory
    
    def set_basedir_name(self, dir0, restart):
        """
        Note 
        -----
        If restart is off, the directory name needs to be determined carefully.
        """
        if restart:
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
    def verbosity(self):
        return self._verbosity
    
    @property
    def yamlfile_for_outdir(self):
        return self._yamlfile_for_outdir
    
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
            ### Angstrom => Bohr
            self._cutoffs = helper.get_cutoffs_automatically(
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
    
    #@property
    #def scell_matrix3(self):
    #    return self._scell_matrix3
    
    @property
    def primitive(self):
        return self._primitive

    @property
    def unitcell(self):
        if self._unitcell is None:
            msg = " Error: unitcell is not yet defined."
            logger.error(msg)
        else:
            return self._unitcell
    
    @property
    def supercell(self):
        """
        """
        if self._supercell is None:
            msg = " Error: supercell is not yet defined."
            logger.error(msg)
        else:
            return self._supercell

    #@property
    #def supercell3(self):
    #    """
    #    """
    #    if self._supercell3 is None:
    #        msg = " Error: supercell for cubic FCs is not yet defined."
    #        logger.error(msg)
    #    else:
    #        return self._supercell3
        
    @property
    def frequency_range(self):
        if self._frequency_range is None:
            fn = self.out_dirs['harm']['bandos'] + '/' + self.prefix + '.bands'
            if os.path.exists(fn) == False:
                msg = " Phonon dispersion needs to be calculated."
                logger.info(msg)
            else:
                fmin, fmax = helper.read_frequency_range(fn, format='anphon')
                self._frequency_range = [fmin, fmax]
        return self._frequency_range
    
    @property
    def minimum_frequency(self):
        return self.get_minimum_frequency(which="both")
    
    def get_minimum_frequency(self, which="both"):
        """ Read minimum frequencies from log files for band (band.log) and DOS
        (dos.log) calculations and return the value
        """
        #from auto_kappa.alamode.log_parser import (
        #        get_minimum_frequency_from_logfile)
        
        if which == "both":
            properties = ["band", "dos"]
        else:
            properties = [which]
        
        fmin = 1e9
        for propt in properties:
            
            logfile = self.out_dirs["harm"]["bandos"] + "/%s.log" % propt
            
            try:
                out = get_minimum_frequency_from_logfile(logfile)
                fmin = min(fmin, out["minimum_frequency"])
            except Exception:
                msg = " Error: %s may contain error." % logfile
                logger.error(msg)
        
        if fmin > 100.:
            ### alamode log file does not show extremely large negative values
            ### properly.
            return -1e6
        else:
            return fmin
    
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
    
    def _get_number_of_free_fcs(self, order: None, verbosity=False):
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
                    ith = int(data[3][-1]) - 1
                    nfc_free[ith] = int(data[-1])
        
        ### print
        if verbosity:
            msg = "\n Number of free FCs :"
            names = {1: "harmonic", 2: "cubic"}
            for ii in nfc_free:
                if ii in names:
                    name = names[ii]
                else:
                    name = "%dth-order" % (ii+1)
                msg += "\n - %s : %s " % (
                        name.ljust(10), str(nfc_free[ii]).rjust(6))
            logger.info(msg)
        
        return nfc_free
    
    def get_suggested_structures(self, order: None, nrandom=0, disp_mode='fd',
            temperature=None, number_of_displacements=1, classical=False,
            can_return_none=False):
        """ Return structures in which atoms are displaced with the given method,
        ``disp_mode``.
        
        order : int 
            1 (harmonic), 2 (cubic), 3, ...
        
        disp_mode : string
            "fd" for finite displacement,
            "random" for simple random displacement,
            "random_normalcoordinate" for random displacement in normal 
            coordinates
        
        nrandom : int
            Number of patterns generated with the random displacement
        
        """
        ### ver.1
        #if order == 1:
        #    structure = self.supercell
        #else:
        #    structure = self.supercell3
        #    #sys.exit()
        #
        ### ver.2
        structure = self.supercell
        
        ### get displacements
        logger.info(f"\n Displacement mode : {disp_mode}\n")
        if disp_mode == "fd":
            file_pattern = self._get_file_pattern(order)
            all_disps = self._get_displacements(
                    "fd", file_pattern=file_pattern, order=order)
        
        elif disp_mode == "random_normalcoordinate":
            logfile = self.out_dirs['harm']['evec'] + '/evec_commensurate.log'
            if wasfinished_alamode(logfile) == False:
                propt = 'evec_commensurate'
                self.write_alamode_input(propt=propt)
                self.run_alamode(propt=propt)
            
            file_evec = self.out_dirs['harm']['evec'] + '/' + self.prefix + '.evec'
            all_disps = self._get_displacements(
                    disp_mode,
                    file_evec=file_evec,
                    temperature=temperature, 
                    number_of_displacements=number_of_displacements,
                    classical=classical)
         
        elif disp_mode == "random":
            all_disps = self._get_displacements(
                    disp_mode,
                    number_of_displacements=number_of_displacements,
                    )
         
        else:
            logger.error(f"\n Error: displacement mode {disp_mode} is not supported")
            sys.exit()
        
        ### error
        if all_disps is None:
            msg = "\n Error: Failed to obtain displacement patterns."
            msg += "\n The structures used in the current calculation may be "
            msg += "incompatible with the ones used in the previous calculation."
            logger.error(msg)
            if can_return_none:
                return None
            else:
                sys.exit()
        
        ## get pristine structure
        structures = {}
        structures['prist'] = structure.copy()
        
        ## set calculator for each pattern
        for i, displacements in enumerate(all_disps):
            scell = structure.copy()
            scell.translate(displacements)
            structures[i+1] = scell
        
        return structures
    
    def _get_displacements(
            self, displacement_mode: None,
            file_pattern=None,
            order=None,
            #
            file_evec=None, number_of_displacements=None, temperature=500., 
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
        try:
            almdisp = AlamodeDisplace(
                    displacement_mode, codeobj,
                    file_evec=file_evec,
                    primitive=self.primitive,
                    verbosity=self.verbosity
                    )
        except:
            ## error which may due to the incompatible structures for 
            ## the structure used and the one used in previous calculations
            almdisp = None
        
        if almdisp is None:
            msg = " Error: Couldn't obtain AlamodeDisplace object properly."
            logger.error(msg)
            return None
        
        msg = "\n"
        msg += " Generate displacement patterns with an Alamode tool\n"
        msg += " Displacement mode : %s" % displacement_mode
        logger.info(msg)
        
        if displacement_mode == 'fd':
            if file_pattern is None:
                logger.error("\n Error: file_pattern must be given.")
                sys.exit()
            
            if order == 1:
                mag = self.magnitude
            elif order == 2:
                mag = self.magnitude2
            
            logger.info(f" Displacement magnitude : {mag:.2f} Ang")
            header_list, disp_list = almdisp.generate(
                    file_pattern=[file_pattern], magnitude=mag,
                    )
            
        elif displacement_mode == 'random_normalcoordinate':
            logger.info(" Temperature : %.1f K" % (temperature))
            header_list, disp_list = almdisp.generate(
                    temperature=temperature,
                    number_of_displacements=number_of_displacements,
                    classical=classical
                    )

        elif displacement_mode == 'random':
            logger.info(" Displacement magnitude : %.2f Ang" % (self.magnitude2))
            header_list, disp_list = almdisp.generate(
                    number_of_displacements=number_of_displacements,
                    magnitude=self.magnitude2
                    )
            
        else:
            logger.warning("\n Error: %s is not supported." % (displacement_mode))
            sys.exit()
        
        all_disps = np.zeros_like(disp_list)
        for i, each in enumerate(disp_list):
            all_disps[i] = np.dot(each, self.supercell.cell)
        
        return all_disps 
    
    def calc_forces(self, order: None, calculator=None, 
            nmax_suggest=100, frac_nrandom=10., 
            temperature=500., classical=False,
            calculate_forces=True, 
            output_dfset=None,
            amin_params={},
            # max_limit_of_estimated_time=24.*30.,
            ):
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
            (``nsuggest``).
            ``nrandom`` = max(``frac_nrandom`` x ``nsuggest``, ``nsuggest``)
        
        output_dfset : int or bool, default=1
            If False or 0, DFSET file is not written.
            If 1, output DFSET file only when it doesn't exist.
            If 2, output DFSET file even when it exists.
        
        """
        # super().__init__()
        
        if output_dfset is not None:
            msg = " Caution: \"output_dfset\" option is no loger used."
        
        self._nmax_suggest = nmax_suggest
        
        self._start_force_calculation(order, calculate_forces)
        
        ### get suggsted structures with ALM
        ### structures : dict of structure objects
        if order <= 2:
            nsuggest = self._get_number_of_suggested_structures(order)
            msg = "\n Number of the suggested structures with ALM : %d" % (nsuggest)
        else:
            ### want to always use LASSO
            nsuggest = self.nmax_suggest + 1
            msg = "\n LASSO is used for high-order FCs."
        logger.info(msg)
        
        ### number of free FCs
        nfcs_list = self._get_number_of_free_fcs(order, verbosity=True)
        nfcs = nfcs_list[min(3, order)]
        natoms = len(self.supercell)
        
        ### output directory
        outdir0 = self._get_base_directory_for_forces(order, nsuggest, nmax_suggest)
        
        if order == 1:
            mag = self.magnitude
        elif order == 2:
            mag = self.magnitude2
        else:
            mag = None
        
        ### Check the directory name for the pristine structure (ver.1.4.0)
        if self.supercell is not None:
            check_directory_name_for_pristine(outdir0, self.supercell)
        
        ### Compressive sensing, LASSO, is used when the number of the 
        ### suggested structure (nsuggest) exceeds nmax_suggest.
        ### When nsuggest > nmax_suggest, max(nsuggest * frac_nrandom,
        ### nmax_suggest) structures will be generated with random displacement
        ### based on the normal coordinates.
        # if order != 1 and nsuggest > nmax_suggest:
        if order == 1:
            fc_type = 'fd'
        elif order == 2:
            fc_type = self.fc3_type
        else:
            fc_type = 'lasso'
        
        if fc_type == 'lasso':
            structures = self._get_suggested_structures_for_lasso(
                order, nfcs, natoms, frac_nrandom, nmax_suggest,
                temperature=temperature, classical=classical)
        else:
            ## if order == 1 or (order == 2 and self.fc3_type == 'fd'):
            structures_tmp = self.get_suggested_structures(order, disp_mode='fd')
            
            ## Adjust keys of structures
            structures = adjust_keys_of_suggested_structures(
                structures_tmp, outdir0, 
                prist_structure=self.supercell,
                mag=mag)
        
        ### If something wrong, return None
        if structures is None:
            return None
        
        ### get AMIN parameters
        pp = default_amin_parameters.copy()
        for key in pp.keys():
            if key in amin_params.keys():
                if amin_params[key] is not None:
                    pp[key] = amin_params[key]
        amin_params_set = pp.copy()
        
        ### start the force calculation
        self._counter_done = 0
        self._counter_calc = 0
        struct_keys = list(structures.keys())
        logger.info("")
        for ii, key in enumerate(struct_keys):
            self._job_for_each_structure(
                ii, structures, outdir0, order, calculator, 
                calculate_forces=calculate_forces, **amin_params_set)
            
        ### output DFSET
        nsuggest =  len(structures)
        # if num_done == len(structures) and output_dfset:
        if self._counter_done == len(structures):
            self._make_dfset_file(order, nsuggest, outdir0)
            
        logger.info("")

    def _write_alamode_input_harmonic(
        self, propt, outdir=None, deltak=None, reciprocal_density=None, **kwargs):
        """ Make an input script for the ALAMODE job for an harmonic property
        """
        from auto_kappa.structure.crystal import get_automatic_kmesh
        
        if propt == "band":
            self.write_alamode_input(
                    propt=propt, deltak=deltak, outdir=outdir, **kwargs)
        elif propt == "dos":
            nkts = get_automatic_kmesh(
                    self.primitive,
                    reciprocal_density=reciprocal_density)
            self.write_alamode_input(
                    propt=propt, nkts=nkts, outdir=outdir, **kwargs)
        elif propt == "fc2":
            self.write_alamode_input(
                    propt=propt, outdir=outdir, **kwargs)
        else:
            msg = "\n Error: %s is not supported." % propt
            logger.error(msg)
    
    def write_alamode_input(
            self, propt=None, order=None, 
            deltak=0.01, kpts=[15,15,15],
            outdir=None, tolerance=1.8897259886e-5,
            **kwargs
            ):
        """ Write Anphon Input
        
        Args
        -----
        propt : string
            "band", "phonons", "kappa", 'kappa_scph', 'kappa_4ph', 'kappa_scph_4ph', 'evec_commensurate'
        
        kpts : list of int, shape=(3)
            k-points
        
        outdir : string
            should be given to use harmonic and cubic FCs obtained using
            supercells of different sizes
        
        deltak : float
            resolution of phonon dispersion [0.01]
        
        kwargs : 
            Any Alamode parameters can be given. The given parameters will be
            updated just before making the input file.
        
        """
        ### get alamode_type and mode
        out = self._get_alamodetype_mode(propt)
        if out is None:
            msg = "\n Error: %s is not supported yet.\n" % propt
        
        alamode_type = out[0]
        mode = out[1]
        
        ### prepare filenames
        files = self._get_filenames(propt, order, **kwargs)
        fc2xml = files['fc2xml']
        fc3xml = files['fc3xml']
        fcsxml = files['fcsxml']
        dfset = files['dfset']
        dir_work = files['dir_work']
        born_xml = files['born_xml']
        
        ## Prepare output directory
        if outdir is not None:
            dir_work = outdir
        os.makedirs(dir_work, exist_ok=True)
        
        #### write output directory in a yaml file
        name = "%s" % (propt)
        if name == "kappa":
            name += "_%dx%dx%d" % (kpts[0], kpts[1], kpts[2])
        
        if self.scell_matrix is not None:
            name += "_%dx%dx%d" % (
                    self.scell_matrix[0][0],
                    self.scell_matrix[1][1],
                    self.scell_matrix[2][2],
                    )
        
        info = {"directory": dir_work.replace(self.base_directory, "."),
                "kind": "ALAMODE",
                f"note": "ALAMODE calculation for {propt}"}
        write_output_yaml(self.yamlfile_for_outdir, name, info)
        
        ## non-analytical term
        if self.nac != 0 and alamode_type.startswith('anphon'):
            
            borninfo = 'BORNINFO'
            
            ## make BORNINFO file
            outfile = dir_work + '/BORNINFO'
            write_born_info(born_xml, outfile=outfile) 
        else:
            borninfo = None
        
        ## set kpmode
        kpmode = None
        if propt in ['band', 'scph']:
            kpmode = 1
        elif propt in ['dos'] or propt.startswith('kappa'):
            kpmode = 2
        elif propt in ['evec_commensurate']:
            kpmode = 0
        
        ### Get ALAMODE input object
        inp = self._get_input_object(
            alamode_type=alamode_type,
            mode=mode, kpmode=kpmode,
            fc2xml=fc2xml, fc3xml=fc3xml, fcsxml=fcsxml, 
            dfset=dfset, borninfo=borninfo)
        
        ####
        self._prefix = inp['prefix']
        
        ################################################
        ### Set parameters for the specific calculation
        self._set_parameters_for_property(
            inp, propt=propt, deltak=deltak, kpts=kpts, order=order)
        
        ### modify path of filenames
        filename_keys = ["fc2xml", "fc3xml", "fcsxml"]
        given_params = kwargs.copy()
        for key in given_params:
            if key in filename_keys:
                if given_params[key][0] == "/":
                    given_params[key] = os.path.relpath(given_params[key], dir_work)
             
        ### make input script for ALAMODE
        inp.update({'tolerance': tolerance})
        inp.update(given_params)
        
        ### Modify parameters for 4 phonon scattering (anphon >= ver.1.9)
        if '4ph' in propt:
            params = self.modify_parameters_for_4ph(inp.as_dict())
            inp = AnphonInput.from_dict(params)
        
        ## Check if the tolerance was changed
        ## Settting of "tolerance" was added from ver.0.4.0.
        filename = dir_work + '/' + propt + '.in'
        if was_primitive_changed(self.unitcell, tol_prev=1e-3*BohrToA, tol_new=1e-5):
            if was_tolerance_changed(filename, inp.as_dict()):
                directory = os.path.dirname(filename)
                backup_previous_results(directory, propt, prefix=self.prefix)
        
        ## Write input file
        if alamode_type.startswith('anphon'):
            ver_alamode = get_version(self.commands['alamode'][alamode_type])
        else:
            ver_alamode = None
        inp.to_file(filename=filename, version=ver_alamode)
        
        # msg = "\n Make an input script for ALAMODE : %s." % (
        #         self.get_relative_path(filename))
        
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
            msg = " Warning: cannot find %s" % fn
            logger.warning(msg)
            return None
    
    def modify_parameters_for_4ph(self, params):
        """ Modify parameters for 4-phonon scattering using anphon >= ver.2.0
        """
        ### fc*xml => fc*file
        params_file = {}
        for key in params.keys():
            if key.startswith('fc') and key.endswith('xml'):
                params_file[key] = params[key]
        for key in params_file:
            params.pop(key, None)
            params[key.replace('xml', 'file')] = params_file[key]
        
        params['nkd'] = None
        params['quartic'] = 1
        params['ismear_4ph'] = 2   # adaptive smearing for 4ph
        params['interpolator'] = 'log-linear' # linear, log-linear, modified-log-linear
        # for adaptive smearing: default 1.0. 
        # Smaller value makes the calculation faster, but less accurate.
        params['adaptive_factor'] = 0.1
        return params
    
    def _get_alamodetype_mode(self, propt):
        
        ### get alamode_type and mode 
        flag = False
        for k1 in self.propts_mode_type:
            for k2 in self.propts_mode_type[k1]:
                if propt in self.propts_mode_type[k1][k2]:
                    return k1, k2
        if flag == False:
            return None
    
    def get_optimal_nac(
            self, tol_neg_frac=0.03,
            deltak=None, max_num_corrections=None, 
            reciprocal_density=None, negative_freq=None):
        """ Analyze harmonic properties with different methods for NAC when
        negative frequencies were obtained. Parameters are required for
        ``analyze_harmonic_property``.
        
        Args
        -----
        tol_neg_frac : float
            Tolerance of the fractional number of negative frequencies w.r.t the
            total k-points (num_kpt). If the number of negative frequencies
            (num_neg) is less than ``tol_neg_frqc * num_kpt``, the optimal NAC
            option will be searched. 
        """
        
        nac_orig = self.nac
        
        ###
        try:
            file_band_orig = (
                    self.out_dirs["harm"]["bandos"] + "/" +
                    self.prefix + ".bands")
            dump = np.genfromtxt(file_band_orig)
            frequencies = dump[:,1]
            num_neg = len(np.where(frequencies < negative_freq)[0])
            num_tot = len(frequencies)
            num_tol = int(num_tot * tol_neg_frac)
            msg = "\n Number of data : %d (negative), %d (total)" % (
                    num_neg, num_tot)
            logger.info(msg)
            if num_neg > num_tol:
                ### If there exists too many negative frequencies, abort the
                ### optimal NAC option.
                return None
        except Exception:
            pass
        
        ### Store old results
        dir1 = self.out_dirs["harm"]["bandos"]
        outdir = dir1 + "/nac_%d" % nac_orig
        if os.path.exists(outdir) == False:
            os.makedirs(outdir, exist_ok=True)
            for ff in glob.glob(dir1+"/*"):
                if os.path.isfile(ff):
                    shutil.move(ff, outdir)
            ###
            msg = "\n >>> Move files in %s to %s." % (
                    self.get_relative_path(dir1),
                    self.get_relative_path(outdir))
            logger.info(msg)
         
        ### Analyze with different NAC options
        for nac in [2, 1, 3, 0]:
            
            if nac == nac_orig:
                continue
            
            self.nac = nac
            outdir = self.out_dirs["harm"]["bandos"] + "/nac_%d" % nac
            for propt in ["band", "dos"]:
                
                fcsxml_abs = (
                        self.out_dirs["harm"]["force"] + "/" + 
                        self.prefix + ".xml")
                fcsxml = os.path.relpath(fcsxml_abs, outdir)
                
                self.analyze_harmonic_property(
                        propt, 
                        outdir=outdir, 
                        max_num_corrections=max_num_corrections,
                        deltak=deltak,
                        reciprocal_density=reciprocal_density,
                        fcsxml=fcsxml,
                        )
            
            ### log file names
            log_band = outdir + "/band.log"
            log_dos = outdir + "/dos.log"
            
            ### get the minimum frequency
            try:
                out_band = get_minimum_frequency_from_logfile(log_band)
                out_dos = get_minimum_frequency_from_logfile(log_dos)
                fmin = min(
                        out_band["minimum_frequency"],
                        out_dos["minimum_frequency"]
                        )
            except Exception:
                fmin = -1e5
            
            if fmin > negative_freq:
                self.nac = nac_orig
                return nac
        
        self.nac = nac_orig
        return None
    
    def analyze_harmonic_property(
            self, propt,
            outdir=None, max_num_corrections=None, neglect_log=None,
            deltak=None, reciprocal_density=None, **kwargs):
        """ Analyze each harmonic property (why not anharmonic?) 
        
        Args
        -----
        propt : string
            "fc2", "band", or "dos"

        kwargs : 
            Any ALAMODE parameters, including "fcsxml".
        """
        ### make an input script for the ALAMODE job
        self._write_alamode_input_harmonic(
            propt, outdir=outdir, deltak=deltak, 
            reciprocal_density=reciprocal_density, **kwargs)
        
        ### get log file name
        if propt == "fc2":
            logfile = self.out_dirs["harm"]["force"] + "/%s.log" % propt
        else:
            logfile = self.out_dirs["harm"]["bandos"] + "/%s.log" % propt
        
        count = 0
        nprocs = self.commands["alamode"]["nprocs"]
        has_error = False
        while True:
            
            if count == 0:
                if neglect_log is None:
                    neg_log = 0
                else:
                    neg_log = neglect_log
            else:
                neg_log = 1
            
            ### calculate phonon property
            self.run_alamode(
                    propt=propt, neglect_log=neg_log, outdir=outdir)
            count += 1
            
            ### >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
            ### check log file
            flag = helper.should_rerun_alamode(logfile)
            
            ### Check phonon band
            ## This part solves a bug in the old calculation.
            ## In old version, eigenvalues has sometimes contained ``nan`` values.
            if flag == False and propt == "band":

                if outdir is None:
                    file_band = (
                            self.out_dirs["harm"]["bandos"] + 
                            "/%s.bands" % self.prefix)
                else:
                    file_band = outdir + "/%s.bands" % self.prefix
                
                flag = helper.should_rerun_band(file_band)
                
                if flag:
                    logger.warning("\n Error in .bands file")
            
            ### the property has been calculated properly
            if flag == False:
                break
            
            ### modify the MPI and OpenMP conditions or finish the job
            ### Note: ALAMODE job sometimes exceeds the memory
            has_error = False
            if count == 1:
                ### Change to OpenMP parallelization
                #ak_log.rerun_with_omp()
                logger.error("\n Rerun the ALAMODE job with OpenMP.")
                self.commands['alamode']['anphon_para'] == "omp"
            
            elif count > 1 and count < max_num_corrections:
                ### modify the number of threads (even number may be better)
                nprocs /= 2
                if nprocs > 1:
                    nprocs += int(nprocs % 2)
                
                self.commands['alamode']['nprocs'] = nprocs
                msg = "\n Rerun the ALAMODE job with "\
                        "OMP_NUM_THREADS/SLURM_CPUS_PER_TASK = %d" % nprocs
                logger.error(msg)

                if nprocs == 1:
                    count = max_num_corrections
                
            else:
                has_error = True
                break
            ### <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
        
        if has_error:
            if max_num_corrections is not None:
                logger.debug(" %d" % max_num_corrections)
            else:
                logger.debug(" None")
            msg = "\n Error: ALAMODE job for %s "\
                    "has not been finished properly." % propt
            msg += "\n A possible solution is to use a large memory node."
            msg += "\n Abort the job."
            logger.error(msg)
            sys.exit()
    
    def run_alamode(self, propt=None, order=None, neglect_log=0, outdir=None, logfile=None):
        """ Run anphon
        
        Args
        ---------
        propt : string
            "band", "dos", "evec_commensurate", "fc*", "kappa", ...

        neglect_log : int
            If False, anphon will not be run if the same calculation had been
            conducted while, if True, anphon will be run forecely.
        
        """
        ### get alamode_type and mode
        out = self._get_alamodetype_mode(propt)
        
        if out is None:
            logger.error(f" Error: {propt} is not supported yet.")
            sys.exit()
        alamode_type = out[0]
        mode = out[1]
        
        ## change directory
        workdir = self._get_work_directory(propt, order=order)
        
        ### modify work directory
        if outdir is not None:
            workdir = outdir
        
        ### print title
        msg = "\n"
        line = "Run %s for %s:" % (alamode_type, propt)
        msg += " " + line + "\n"
        msg += " " + "-" * (len(line)) + "\n\n"
        msg += " Working directory : " + self.get_relative_path(workdir)
        logger.info(msg)
        
        filename = f"{propt}.in"
        if logfile is None:
            logfile = f"{propt}.log"
        
        ## prepare command and environment
        # print(propt, alamode_type, mode)
        if (alamode_type.startswith('anphon') and 
                self.commands['alamode']['anphon_para'] == 'mpi'):
            _nprocs = self.commands['alamode']['nprocs']
            _nthreads = 1
        else:
            _nprocs = 1
            _nthreads = self.commands['alamode']['nprocs']
        
        val = run_alamode(filename, logfile, workdir=workdir, 
                neglect_log=neglect_log,
                mpirun=self.commands['alamode']['mpirun'], 
                nprocs=_nprocs,
                nthreads=_nthreads,
                command=self.commands['alamode'][alamode_type],
                )
        
        if val == -1:
            logger.info(f" {propt} has already been calculated.")
        elif val == 1:
            msg = "\n Error : ALAMODE job was not finished properly."
            msg += "\n Stop the calculation."
            logger.error(msg)
            sys.exit()
        
        if mode == 'optimize' and propt in ['lasso', 'fc2', 'fc3']:
            
            ### copy the generated FCs file to "result" directory
            self._copy_generated_fcsfiles(propt=propt, order=order)
            
            if propt in ['lasso']:
                self._plot_cvsets(order=order)

    def _copy_generated_fcsfiles(self, propt=None, order=None):
        """ Copyr a FCs file into the "result" directory.
        """
        if propt == 'lasso':
            
            #if newversion == False:
            #    fn1 = self.out_dirs['higher']['lasso']+'/'+self.prefix+'.xml'
            #    fn2 = self.out_dirs['result']+'/'+self.outfiles['higher_xml']
            #else:
            ##
            if order == 2:
                fn1 = self.out_dirs['cube']['lasso']+'/'+self.prefix+'.xml'
                fn2 = self.out_dirs['result']+'/'+self.outfiles['cube_lasso_xml']
            else:
                fn1 = self.out_dirs['higher']['lasso']+'/'+self.prefix+'.xml'
                fn2 = self.out_dirs['result']+'/'+self.outfiles['higher_xml']
        
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
            msg = "\n Error: %s is not supported yet." % propt
            logger.error(msg)
            sys.exit()
        
        ##
        if os.path.exists(fn1) == False:
            
            msg = ' %s does not exist.' % self.get_relative_path(fn1)
            logger.info(msg)
        
        else:
            if os.path.exists(fn2):
                msg = "\n"
                msg += " %s was overwritten." % (self.get_relative_path(fn2))
            else:
                msg = "\n"
                msg += " %s was created." % self.get_relative_path(fn2)
            logger.info(msg)
            shutil.copy(fn1, fn2)
    
    #@property
    #def file_result(self):
    #    return self._file_result

    #@property
    #def file_isotope(self):
    #    return self._file_isotope

    @property
    def scat(self):
        return self._scat
    
    def set_scattering_info(
            self, grain_size=None, temperature=300, calc_type="cubic", verbosity=True):
        
        ### directory name for kappa
        dirs_kappa = self.get_kappa_directories(calc_type=calc_type)
        
        ## get directory name for kappa with the finest k grid
        dir_kappa = None
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
        file_result = dir_kappa + '/%s.result' % (self.prefix)
        file_isotope = dir_kappa + '/%s.self_isotope' % (self.prefix)
        
        msg1 = "\n Scattering info is obtained from the following files:"
        if os.path.exists(file_result) == False:
            msg_error = " Error: %s does not exist." % (
                    self.get_relative_path(file_result))
            logger.error(msg_error)
        else:
            msg1 += "\n - %s." % (self.get_relative_path(file_result))
        
        if os.path.exists(file_isotope) == False:
            msg_error = " Note : %s does not exist." % (
                    self.get_relative_path(file_isotope))
            logger.error(msg_error)
            file_isotope = None
        else:
            msg1 += "\n - %s." % (self.get_relative_path(file_isotope))
        
        if verbosity:
            logger.info(msg1)
        
        ###
        from auto_kappa.alamode.analyzer.scattering import Scattering
        self._scat = Scattering(
                file_result, 
                file_isotope=file_isotope,
                grain_size=grain_size, 
                temperature=temperature
                )
    
    def write_lifetime_at_given_temperature(self, temperature=300, outfile=None, grain_size=None):
        """ Write lifetime at a given temperature in a csv file in a format similar to ALAMODE.
        
        Args
        ------
        temperature : float
            Temperature in K.
        
        outfile : string
            Output file name.
        
        grain_size : float
            Grain size in nm.
            
        """
        self.set_scattering_info(temperature=temperature, grain_size=grain_size, verbosity=False)
        
        dump = {"ik": [], "is": [], 
                "frequency[cm^-1]": [], "lifetime[ps]": [], 
                "|velocity|[m/s]": [], "MFP[nm]": [], "multiplicity": [],
                "kxx": [], "kxy": [], "kxz": [],
                "kyx": [], "kyy": [], "kyz": [],
                "kzx": [], "kzy": [], "kzz": []}
        
        ## print(self.scat.result)
        nk = len(self.scat.result['frequencies'])
        nbands = len(self.scat.result['frequencies'][0])
        for ik in range(nk):
            multiplicity = self.scat.result['multiplicity'][ik]
            for ib in range(nbands):
                
                frequency = self.scat.result['frequencies'][ik][ib] ## cm^-1
                lifetime = self.scat.lifetime[ik][ib]  ## ps
                
                ## velocity
                velo_all = self.scat.result['velocities'][ik][ib]  ## shape(multiplicity, 3)
                velo_all = np.asarray(velo_all)
                assert len(velo_all) == multiplicity
                
                vave = 0.
                for im in range(multiplicity):
                    vave += np.linalg.norm(velo_all[im]) / multiplicity
                
                ## kappa at each mode
                kappa_each = np.zeros((3, 3))
                for i1 in range(3):
                    for i2 in range(3):
                        kappa_each[i1, i2] += (
                            self.scat.kmode[ik,ib,i1,i2] / multiplicity
                        )
                
                dump["ik"].append(ik + 1)
                dump["is"].append(ib + 1)
                dump["frequency[cm^-1]"].append(frequency)
                dump["lifetime[ps]"].append(lifetime)
                dump["|velocity|[m/s]"].append(vave)
                dump["MFP[nm]"].append(lifetime * vave)
                dump["multiplicity"].append(multiplicity)
                for i, j in itertools.product(range(3), repeat=2):
                    dump["k" + "xyz"[i] + "xyz"[j]].append(kappa_each[i, j])
            
        ### write a csv file
        if outfile is None:
            outfile = self.out_dirs['result'] + '/tau_%dK.csv' % (int(temperature))
        
        df = pd.DataFrame(dump)
        with open(outfile, 'w') as f:
            f.write("# Created by auto-kappa in a format similar to ALAMODE\n")
            f.write("# Input file   : %s\n" % self.get_relative_path(self.scat.result.filename))
            if self.scat.file_isotope is not None:
               if os.path.exists(self.scat.file_isotope):
                    f.write("# Itotope      : %s\n" % self.get_relative_path(self.scat.file_isotope))
            if grain_size is not None:
                f.write("# Grain size   : %f nm\n" % grain_size)
            f.write("# Temperature  : %d K\n" % temperature)
            f.write("# kpoint range : 1 %d\n" % nk)
            f.write("# mode   range : 1 %d\n" % nbands)
            df.to_csv(f, index=False, float_format='%.7e')
            msg = " Output %s" % self.get_relative_path(outfile)
            logger.info(msg)
        
    
    def plot_lifetime(
            self, temperatures="300:500", calc_type="cubic", outputfile=True):
        
        data = temperatures.split(":")
        ts = [float(t) for t in data]
        
        dfs = {}
        for t in ts:
            
            self.set_scattering_info(temperature=t, grain_size=None, verbosity=False)
            omegas = self.scat.result['frequencies']
            taus = self.scat.lifetime
            
            nk = len(omegas)
            nbands = len(omegas[0])
            
            dfs[int(t)] = pd.DataFrame()
            dfs[int(t)]['frequency'] = omegas.reshape(nk*nbands)
            dfs[int(t)]['lifetime'] = taus.reshape(nk*nbands)
        
        ### figname
        prefix = self.out_dirs['result'] + '/fig_lifetime'
        if calc_type != "cubic":
            prefix += "_%s" % calc_type
        figname = prefix + ".png"
        
        ### plot a figure
        from auto_kappa.plot.pltalm import plot_lifetime
        out = plot_lifetime(dfs, figname=self.get_relative_path(figname))
        
        return out
    
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
            
            if self.scat.scattering_rates[key] is None:
                continue
            
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
        plot_scattering_rates(
                frequencies, scat_rates, labels, 
                figname=self.get_relative_path(figname))
        
    def _plot_cvsets(self, order=None):
        
        from auto_kappa.plot.pltalm import plot_cvsets
        
        msg = "\n ### Plot CV results ###"
        logger.info(msg)
        
        if order == 2:
            figname = self.out_dirs['result'] + '/fig_cvsets_cube.png'
            plot_cvsets(
                    directory=self.out_dirs['cube']['cv'],
                    figname=figname)
        else:
            figname = self.out_dirs['result'] + '/fig_cvsets_high.png'
            plot_cvsets(
                    directory=self.out_dirs['higher']['cv'], 
                    figname=figname)
        
        logger.info("")
    
    def plot_bandos(self, **kwargs):
        
        ### set figure name
        if 'figname' not in kwargs.keys():
            figname = self.out_dirs['result'] + '/fig_bandos.png'
        else:
            figname = kwargs['figname']
        
        from auto_kappa.plot.bandos import plot_bandos

        ### output title
        line = "Plot band and DOS:"
        msg = "\n " + line + "\n"
        msg += " " + "-" * len(line) + "\n"
        logger.info(msg)
        
        fig = plot_bandos(
                directory=self.get_relative_path(
                    self.out_dirs['harm']['bandos']), 
                prefix=self.prefix, 
                figname=self.get_relative_path(figname),
                **kwargs
                )
    
    def get_kappa_directories(self, calc_type="cubic"):
        
        if calc_type == "cubic":
            ll_tmp = self.out_dirs['cube']["kappa_%s" % self.fc3_type] + "_"
        elif calc_type == "scph":
            ll_tmp = self.out_dirs["higher"]["kappa"] + "_"
        
        line = ll_tmp + "*"
        
        dirs = glob.glob(line)
        dirs_done = {}
        for dd in dirs:
            fn_log = dd + "/kappa.log"
            if os.path.exists(fn_log) == False:
                continue
            if wasfinished_alamode(fn_log):
                label = dd.split(ll_tmp)[1]
                dirs_done[label] = dd
        return dirs_done
    
    def plot_kappa(self, figname=None, calc_type="cubic"):
        """
        Args
        ======

        calc_type : string
            "cubic", "scph"

        """
        
        dirs_kappa = self.get_kappa_directories(calc_type=calc_type)
        
        keys = dirs_kappa.keys()
        
        dfs = {}
        logger.info("")
        for i, key in enumerate(keys):
            try:
                dir1 = self.get_relative_path(dirs_kappa[key])
                msg = " Read %s" % (dir1)
                logger.info(msg)
                dfs[key] = helper.read_kappa(dirs_kappa[key], self.prefix)
            except Exception:
                continue
        
        from auto_kappa.plot.pltalm import plot_kappa
        if figname is None:
            figname = self.out_dirs['result'] + '/fig_kappa.png'
        logger.info("")
        plot_kappa(dfs, figname=self.get_relative_path(figname))
        return 0
    
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
            
            ## Calculate cumulative and spectral kappa at a given temperature
            self.scat.change_temperature(t)
            dfs[int(t)] = self.scat.get_cumulative_kappa(
                    temperature=t, wrt=wrt, xscale=xscale, nbins=nbins)
            
            ## output file for the cumulative and spectral kappa
            try:
                df_each = dfs[int(t)].rename(columns={'xdat': wrt})
                outfile = self.out_dirs['result'] + '/kspec_%s_%dK.csv' % (wrt, int(t))
                
                ### save spectral kappa
                with open(outfile, 'w') as f:
                    
                    ## comment
                    f.write("# Temperature : %d\n" % t)
                    if self.scat.size is not None:
                        f.write("# Grain size  : %f\n" % self.scat.size)
                    
                    df_each.to_csv(f, index=False, float_format='%.6e')
                    msg = " Output %s" % self.get_relative_path(outfile)
                    logger.info(msg)
            except:
                pass
        
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
        plot_cumulative_kappa(
                dfs, xlabel=xlabel, 
                figname=self.get_relative_path(figname), 
                xscale=xscale, ylabel=ylabel1)
    
    def calculate_pes(self, negative_freq=-1e-3):
        """ Calculate PES """
        
        ### start PES
        line = "Calculate PES"
        msg = "\n " + line
        msg += "\n " + "=" * len(line)
        logger.info(msg)
        
        ### get a representative kpoint
        from auto_kappa.alamode.pes import get_representative_kpoint
        
        log_band = self.out_dirs["harm"]["bandos"] + "/band.log"
        log_commensurate = self.out_dirs["harm"]["evec"] + "/evec_commensurate.log"
        log_dos = self.out_dirs["harm"]["bandos"] + "/dos.log"
        
        out = get_representative_kpoint(
                log_band=log_band, log_dos=log_dos,
                log_commensurate=log_commensurate,
                epsilon=1e-5, negative_freq=negative_freq)
        
        if out is None:
            return None
        
        kp_type, kpoint = out
        
        ### harmonic FCs
        fc2xml = self.out_dirs['result'] + "/" + self.outfiles['harm_xml']
        if os.path.exists(fc2xml) == False:
            msg = "\n Error: cannot find %s" % self.get_relative_path(fc2xml)
            logger.error(msg)
            return None
        
        ### move directory
        dir_init = os.getcwd()
        dir_pes = self.out_dirs["harm"]["pes"]
        os.makedirs(dir_pes, exist_ok=True)
        os.chdir(dir_pes)
        
        ### make BORNINFO
        if self.nac != 0:
            xml4born = self.out_dirs["nac"] + "/vasprun.xml"
        else:
            xml4born = None
        
        ### calculate PES
        from auto_kappa.alamode.pes import calculate_pes
        calculate_pes(
                self.primitive, kpoint, 
                outdir=self.out_dirs["harm"]["pes"], fcsxml=fc2xml, 
                command=self.commands["alamode"],
                nac=self.nac, vasp_xml=xml4born)
        
        ### move back to the initial directory
        os.chdir(dir_init)

