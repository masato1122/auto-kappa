"""
To Do
=============
1. alm may be faster than ALM.
* LASSO may need to be used with alm command, but not from python command.

2. Use different (larger) supercell size for 3rd-order force constants

"""
# -*- coding: utf-8 -*-
import os
import numpy as np
import warnings
import logging
import datetime

import ase, pymatgen
from ase.build import make_supercell
import subprocess
import pandas as pd

from . import output_directories, output_files
from .units import AToBohr, BohrToA
from .structure.format import change_structure_format
from .io.vasp import wasfinished, get_dfset, print_vasp_params
from .calculator import run_vasp

from .io.vasp import write_born_info
from .io.alm import AlmInput, AnphonInput

class AlmCalc():
    
    def __init__(self, prim_given, mpid='mp', 
            scell_matrix=None, primitive_matrix=None,
            scell_matrix3=None,
            cutoffs=None, nbody=[2,3,3,2,2], mag=0.01,
            cutoff2=-1, cutoff3=4.3, 
            nac=None,
            ncores=1,
            verbosity=0,
            command={
                'mpirun': 'mpirun', "nprocs": 1, 
                "nthreads": 1, "anphon": "anphon",
                }
            ):
        """ Calculations with ALM and Anphon are managed with this class.
         
        Args
        =======
        prim_given : primitive structure
            Different formats such as pymatgen and ASE are accepted while the
            format is changed to ase-Atoms in this module.

        mpid : string
            Material ID, which is used just to determine the directory name.
            It can be anything such as the composition.
        
        scell_matrix : array of float, shape=(3,3)
            supercell matrix wrt unitcell

        primitive_matrix : array of float, shape=(3,3)
            primitive matrix wrt unitcell

        scell_matrix3 : array of float, shape=(3,3)
            supercell matrix wrt unitcell.
            This supercell is used for cubic FCs.

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
            If nac=0, nac is not considered while, if nac=1, nac is considered.
    
        verbosity : int
            If verbosity=1, details are output while, if 0, not.

        command : dict
            Number of processes and threads are given.
            default={
                'mpirun': 'mpirun', "nprocs": 1, 
                "nthreads": 1, "anphon": "anphon",
                }
        
        Example
        ----------
        >>> almcalc = AlmCalc(
        >>>     primitive,
        >>>     mpid='mpid-149',
        >>>     primitive_matrix=pmat,
        >>>     scell_matrix=smat,
        >>>     cutoff2=-1, cutoff3=4.3,
        >>>     nbody=[2,3,3,2], mag=0.01,
        >>>     nac=1,
        >>>     command={'mpirun':'mpirun', 'nprocs':1, 'nthreads':1,'anphon':'anphon'},
        >>>     verbosity=0
        >>>     )
        >>> calc_force = apdb.get_calculator('force', 
        >>>     "./mp-149/harm/force", [4,4,4],
        >>>     )
        >>> almcalc.calc_forces(1, calc_force)
        >>> almcalc.calc_harmonic_force_constants()
        >>> almcalc.write_alamode_input(propt='dos')
        >>> almcalc.run_anphon(propt='dos')
        """
        ### output directories
        dir_init = os.getcwd()
        self.out_dirs = {}
        for k1 in output_directories.keys():
            values1 = output_directories[k1]
            if type(values1) == str:
                self.out_dirs[k1] = dir_init + '/' + mpid + '/' + values1
            else:
                self.out_dirs[k1] = {}
                for k2 in values1.keys():
                    values2 = values1[k2]
                    self.out_dirs[k1][k2] = dir_init + '/' + mpid + '/' + values2
        
        ### structures
        self._primitive = change_structure_format(prim_given, 'ase')
        self._scell_matrix = scell_matrix
        if scell_matrix3 is None:
            self._scell_matrix3 = scell_matrix
        else:
            self._scell_matrix3 = scell_matrix3
        self._primitive_matrix = primitive_matrix
        
        self._unitcell = None
        self._supercell = {'structure': None, 'type': None}
        self._supercell3 = {'structure': None, 'type': None}
         
        self._num_elems = None
        
        ### set cutoff radii: Ang => Bohr
        self._cutoff2 = cutoff2
        self._cutoff3 = cutoff3
        if cutoffs is not None:
            self.cutoffs = cutoffs * AToBohr 
        else:
            n = self.num_elems
            self.cutoffs = np.asarray([
                np.ones((n,n)) * self._cutoff2 * AToBohr, 
                np.ones((n,n)) * self._cutoff3 * AToBohr
                ])
        ##
        self.nbody = nbody
        self.mag = mag            ### Ang
        
        self.force_calculators = {}
        
        self.nac = nac

        self.lasso = False

        self._nmax_suggest = None

        self._command = command

        self._ncores = ncores
        
        ###
        self._prefix = None      ## prefix for input files
        self._frequency_range = None
        
        self.verbosity = verbosity
    
    @property
    def ncores(self):
        return self._ncores
    
    #@ncores.setter
    #def ncores

    @property
    def command(self):
        return self._command
    
    def update_command(self, val):
        """ Command is update """
        self._command.update(val)
    
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
    
    def set_primitive_matrix(self, matrix):
        self._primitive_matrix = matrix
    
    @property
    def scell_matrix(self):
        return self._scell_matrix

    def set_scell_matrix(self, matrix):
        self._scell_matrix = matrix
    
    @property
    def scell_matrix3(self):
        return self._scell_matrix3

    def set_scell_matrix3(self, matrix):
        self._scell_matrix3 = matrix
    
    @property
    def primitive(self):
        return self._primitive

    @property
    def unitcell(self):
        if self._unitcell is None:
            structure = self.primitive.copy()
            
            ### for ase Atoms
            self._unitcell = make_supercell(
                    structure,
                    np.linalg.inv(self.primitive_matrix)
                    )
                
            ### for pymatgen Structure
            #structure.make_supercell(
            #        np.linalg.inv(self.primitive_matrix)
            #        )
            #self._unitcell = structure
        return self._unitcell

    @property
    def supercell(self):
        """
        self._supercell['structure'] : ASE Atoms
        self._supercell['type'] : string
            None, "ase", or "xml"
            "ase" indicates the structure is obtained with "make_supercell"
            module of ASE.
            "xml" indicates the structure is obtained with "vasprun.xml" for
            force calculation in ./mp-***/harm/force/prist.
        """
        if self._supercell['type'] is not None:
            if self._supercell['type'].lower() == 'xml':
                return self._supercell['structure']
        
        ##
        xml = self.out_dirs['harm']['force'] + '/prist/vasprun.xml'
        if os.path.exists(xml):
            self._supercell['structure'] = ase.io.read(xml, format='vasp-xml')
            self._supercell['type'] = 'xml'
        else:
            if self._supercell['structure'] is None:
                structure = self.unitcell.copy()
                
                ### for ase Atoms
                self._supercell['structure'] = make_supercell(
                        structure, self.scell_matrix)
                self._supercell['type'] = 'ase'
        
        return self._supercell['structure']
    
    @property
    def supercell3(self):
        """
        self._supercell3['structure'] : ASE Atoms
        self._supercell3['type'] : string
            None, "ase", or "xml"
        
        Read descriptions of supercell for the detail.
        """
        if self._supercell3['type'] is not None:
            if self._supercell3['type'].lower() == 'xml':
                return self._supercell3['structure']
        
        ##
        xml = self.out_dirs['cube']['force'] + '/prist/vasprun.xml'
        if os.path.exists(xml):
            self._supercell3['structure'] = ase.io.read(xml, format='vasp-xml')
            self._supercell3['type'] = 'xml'
        else:
            if self._supercell3['structure'] is None:
                structure = self.unitcell.copy()
                
                ### for ase Atoms
                self._supercell3['structure'] = make_supercell(
                        structure, self.scell_matrix3)
                self._supercell3['type'] = 'ase'
        
        return self._supercell3['structure']
    
    @property
    def prefix(self):
        return self._prefix
    
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
        
        ### get displacements
        print("")
        print(" Displacement mode :", disp_mode)
        print("")
        if disp_mode == "fd":
            
            all_disps = get_finite_displacements(
                    structure,
                    order=order, 
                    cutoffs=self.cutoffs[:order],  ## Bohr
                    nbody=self.nbody,
                    mag=self.mag                   ## Ang
                    )
        
        elif disp_mode == "random_normalcoordinate":
            
            ### 
            propt = 'evec_commensurate'
            self.write_alamode_input(propt=propt)
            self.run_anphon(propt=propt)

            ###
            fevec = self.out_dirs['lasso']['evec'] + '/' + self.prefix + '.evec'
            all_disps = self._get_random_displacements_normal_coordinate(
                    fevec,
                    temperature=temperature, 
                    number_of_displacements=number_of_displacements,
                    classical=classical)
         
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
    
    def _get_random_displacements_normal_coordinate(self, file_evec: None,
            number_of_displacements=1, temperature=500., classical=False,
            ):
        from .alamode_tools.VASP import VaspParser
        from .alamode_tools.GenDisplacement import AlamodeDisplace
        
        codeobj = VaspParser()
        codeobj.set_initial_structure(self.supercell)
        
        ###########################################################
        ###########################################################
        ###########################################################
        ##
        ## This part need to be fixed.
        ##
        #print(self.scell_matrix)
        #print(self.primitive.get_positions()[:5])
        #print(self.supercell.get_positions()[:5])
        #print(self.supercell.cell)
        #aa = np.dot(self.supercell.cell, np.linalg.inv(self.primitive.cell))
        #print(aa)
        #exit()
        #file_prim = "./mp-160/harm/force/prist/POSCAR"
        almdisp = AlamodeDisplace(
                'random_normalcoordinate', codeobj,
                file_evec=file_evec,
                #file_primitive=file_prim,
                primitive=self.primitive,
                verbosity=self.verbosity
                )
        if almdisp is None:
            return None
        ###########################################################
        ###########################################################
        ###########################################################

        msg = "\n"
        msg += " Generate random displacements with an Alamode tool\n"
        msg += " " + str(datetime.datetime.now()) + "\n"
        msg += "\n"
        msg += " %d patterns will be generated.\n" % (number_of_displacements)
        msg += "\n"
        print(msg)
        
        header_list, disp_list = almdisp.generate(
                temperature=temperature,
                number_of_displacements=number_of_displacements,
                classical=classical
                )
        print("")
        
        all_disps = np.zeros_like(disp_list)
        for i, each in enumerate(disp_list):
            all_disps[i] = np.dot(each, self.supercell.cell)
        return all_disps 
    
    def calc_forces(self, order: None, calculator: None, 
            directory=None, force=False, 
            nmax_suggest=200, frac_nrandom=0.02, 
            temperature=500., classical=False,
            ):
        """ Calculate forces for harmonic or cubic IFCs.
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
        
        directory : string
            If it's given, calculationts will be done in 
            {directory}/{prist, 1, 2, ...}.
        
        force : bool, [False]
            If False, calculation will be skipped if it had been already done.
            If "1", calculation will be done.
            If "2", only VASP files are created.
            
        nmax_suggest : int, default 200
            If the number of suggested structures exceeds 
        
        frac_nrandom : float, default 0.02
            Ratio of the number of patterns generated with a random displacement
            (``nrandom``) to that of the suggested patterns with ALM
            (``nsuggest``). \\
            ``nrandom``:math:`= \\max(` 
                ``frac_nrandom`` :math:`\\times` ``nsuggest``,
                ``nsuggest`` :math:`)`
        
        """
        print("")
        msg = " Force calculation (order: %s)" % order
        border = "=" * (len(msg) + 2)
        print(msg)
        print(border)
        
        self._nmax_suggest = nmax_suggest
        
        ### get suggsted structures with ALM
        ### structures : dict of structure objects
        lasso = 0
        structures = self.get_suggested_structures(order, disp_mode='fd')
        nsuggest =  len(structures)
        
        ###
        msg = "\nNumber of the suggested structures with ALM : %d\n" % (nsuggest)
        print(msg)
        
        ### output directory
        if directory is not None:
            outdir0 = directory
        elif calculator.directory is not None:
            if order == 1:
                outdir0 = self.out_dirs['harm']['force']
            elif order == 2:
                if nsuggest <= nmax_suggest:
                    outdir0 = self.out_dirs['cube']['force']
                else:
                    outdir0 = self.out_dirs['lasso']['force']
        else:
            warnings.warning(" WARNING: Output directory is not defined.")
        
        ### Compressive sensing, LASSO, is used when the number of the 
        ### suggested structure is too large.
        if order != 1 and nsuggest > nmax_suggest:

            self.lasso = True

            nrandom = int(frac_nrandom * nsuggest + 0.5)
            ngenerated = max(nmax_suggest, nrandom)
             
            msg = "\n"
            msg += " Maximum limit of the number of suggested patterns : %d\n" % (nmax_suggest)
            msg += "\n"
            msg += " ### Number of suggested patterns exceeds the maximum limit.\n"
            msg += "\n"
            msg += " ### Number of patterns ###\n"
            msg += " * fractional number of the random patterns : %.3f\n" % (frac_nrandom)
            msg += " * random    : %d\n" % (nrandom)
            msg += " * generated : %d\n" % (ngenerated)
            print(msg)
            
            structures = self.get_suggested_structures(
                    order, 
                    disp_mode='random_normalcoordinate',
                    number_of_displacements=ngenerated,
                    temperature=temperature,
                    classical=classical
                    )
            
        nsuggest =  len(structures)
        for ii, key in enumerate(list(structures.keys())):
            
            outdir = outdir0 + '/' + str(key)
            if wasfinished(outdir) and force == False:
                print(" %s: skipped" % outdir)
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
            if force == 2:
                os.makedirs(calculator.directory, exist_ok=True)
                calculator.write_input(structure)
            else: 
                run_vasp(calculator, structure, method='custodian')
            print(" %s" % (outdir))
        
        print("")
    
    def calc_harmonic_force_constants(self, output=True):
        """ calculate harmonic FCs and return ALM 
        
        Args
        -----
        output : bool
            If True, fcs.xml file is stored, if False, not.
        
        """
        msg = "\n"
        msg += " ### Calculate harmonic force constants\n"
        msg += " " + str(datetime.datetime.now()) + "\n"
        msg += "\n"
        print(msg)
        
        ### prepare directory and file names
        os.makedirs(self.out_dirs['result'], exist_ok=True)
        
        directory = self.out_dirs['harm']['force']
        
        out_dfset = self.out_dirs['result'] + '/' + output_files['harm_dfset']
        if output:
            out_xml = self.out_dirs['result'] + '/' + output_files['harm_xml']
        else:
            out_xml = None
        
        ## get dfset
        offset_xml = self.out_dirs['harm']['force'] + '/prist/vasprun.xml'
        
        if os.path.exists(out_dfset) == False:
            disps, forces = get_dfset(
                    directory, offset_xml=offset_xml, outfile=out_dfset)
        else:
            natoms = len(self.supercell)
            data = np.loadtxt(out_dfset).reshape((-1, natoms, 6))
            
            nstruct = len(data)
            
            disps = np.zeros((nstruct, natoms, 3))
            forces = np.zeros((nstruct, natoms, 3))
            for ii in range(nstruct):
                disps[ii,:,:] = data[ii,:,:3]
                forces[ii,:,:] = data[ii,:,3:]
            
        ## get structure
        structure = ase.io.read(offset_xml)
        
        ## calculate force constants
        ## out = [fc2_values, elem2_indices]
        os.environ['OMP_NUM_THREADS'] = str(self.ncores)
        out = run_alm(
                structure, 1, self.cutoffs[0], [self.nbody[0]], 
                mode='optimize',
                displacements=disps, forces=forces, 
                outfile=out_xml,
                verbosity=self.verbosity
                )
        os.environ.pop('OMP_NUM_THREADS', None)
        return out
     
    def calc_anharm_force_constants(self, maxorder=5):
        """ calculate cubic IFCs and return ALM 
        
        Args
        ----
        maxorder : int, default 5
            order for lasso calculation
        """
        ### prepare directory and file names
        os.makedirs(self.out_dirs['result'], exist_ok=True)
        
        mode = 'optimize'
        
        if self.lasso == False:

            order = 2
            cutoffs = self.cutoffs[:2].copy()

            directory = self.out_dirs['cube']['force']
            out_dfset = self.out_dirs['result'] + '/' + output_files['cube_dfset']
            out_xml = self.out_dirs['result'] + '/' + output_files['cube_xml']
            offset_xml = self.out_dirs['cube']['force'] + '/prist/vasprun.xml'
            
        else:
            
            order = maxorder
            cutoffs = []
            n = self.num_elems
            for i in range(order):
                if i < len(self.cutoffs):
                    cutoffs.append(self.cutoffs[i])
                else:
                    cutoffs.append(np.ones((n,n)) * self._cutoff3 * AToBohr)
            cutoffs = np.asarray(cutoffs)
            
            ### update cutoffs
            self.cutoffs = cutoffs.copy()
            
            directory = self.out_dirs['lasso']['force']
            out_dfset = self.out_dirs['result'] + '/' + output_files['lasso_dfset']
            out_xml = self.out_dirs['result'] + '/' + output_files['lasso_xml']
            offset_xml = self.out_dirs['lasso']['force'] + '/prist/vasprun.xml'
            
        ### print    
        msg = "\n"
        msg += " ### Calculate anharmonic FCs (order=%d)\n" % (order)
        msg += " " + str(datetime.datetime.now()) + "\n"
        msg += "\n"
        print(msg)
        
        ## set nbody
        nbody = []
        for i in range(order):
            if i < len(self.nbody) - 1:
                nbody.append(self.nbody[i])
            else:
                if i <= 1:
                    nbody.append(i+2)
                elif i == 2:
                    nbody.append(3)
                else:
                    nbody.append(2)
        ## update
        self.nbody = nbody.copy()
        
        ## get dfset
        disps, forces = get_dfset(
                directory, 
                offset_xml=offset_xml, 
                outfile=out_dfset)
        
        ### [lasso] extract only required data
        if self.lasso and self.nmax_suggest is not None:
            if len(disps) > self.nmax_suggest:
                disps = disps[:self.nmax_suggest]
                forces = forces[:self.nmax_suggest]
        
        ## get fc2 and elem2
        fc2_values, elem2_indices = \
                self.calc_harmonic_force_constants(output=False)
        
        ## structure
        structure = ase.io.read(offset_xml)
        
        ## 3rd order model
        ## out = [fcs_values, fcs_indices]
        os.environ['OMP_NUM_THREADS'] = str(self.ncores)
        out = run_alm(
                structure, order, cutoffs, nbody, mode='optimize',
                displacements=disps, forces=forces,
                fc2info=[fc2_values, elem2_indices],
                outfile=out_xml, lasso=self.lasso,
                verbosity=self.verbosity)
        os.environ.pop('OMP_NUM_THREADS', None)
        return out
    
    
    #def write_alm_input(self, propt='lasso', maxorder=5, **kwargs):
    #    """ Write an input file for Alm
    #
    #    Args
    #    ------
    #    propt : string
    #        "cv", "lasso", ...
    #    
    #    """
    #    fc2xml = '../../result/' + output_files['harm_xml']
    #
    #    ### prepare parameters
    #    if propt == 'lasso' or propt == 'cv':
    #        dir_work = self.out_dirs['lasso'][propt]
    #        mode = 'optimize'
    #    
    #    ## set input file for anphon
    #    alminp = AlmInput.from_structure(
    #            self.supercell,
    #            norder=maxorder,
    #            fc2xml=fc2xml,
    #            #nonanalytic=self.nac, borninfo=borninfo
    #        )
    #
    #    outfile = propt + '.in'
    #    alminp.to_file(outfile)
    #    print("HEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEERE") 
    #    exit()

    def write_alamode_input(self, propt='band', kpts=[15,15,15], **kwargs):
        """ Write Anphon Input

        Args
        -----
        propt : string
            "band", "phonons", or "RTA"
        
        kpts : list of int, shape=(3)
            k-points
        """ 
        ### prepare filenames
        if propt == 'band' or propt == 'dos':
            
            dir_work = self.out_dirs['harm']['bandos']
            fcsxml = '../../result/' + output_files['harm_xml']
            mode = 'phonons'
            alamode_type = 'anphon'

        elif propt == 'evec_commensurate':
            
            dir_work = self.out_dirs['lasso']['evec']
            fcsxml = '../../result/' + output_files['harm_xml']
            mode = 'phonons'
            alamode_type = 'anphon'
        
        elif propt == 'kappa':
            
            if self.lasso == False:
                dir_work = self.out_dirs['cube']['kappa']
                fcsxml = '../../result/' + output_files['cube_xml']
            else:
                dir_work = self.out_dirs['lasso']['kappa']
                fcsxml = '../../result/' + output_files['lasso_xml']
            
            mode = 'RTA'
            alamode_type = 'anphon'

        elif propt == 'lasso' or propt == 'cv':
            
            dir_work = self.out_dirs['lasso'][propt]
            dfset = '../../result/' + output_files['lasso_dfset']
            fc2xml = '../../result/' + output_files['harm_xml']
            mode = 'optimize'
            alamode_type = 'alm'
        
        ##
        born_xml = self.out_dirs['nac'] + '/vasprun.xml'
        
        ## prepare directory
        os.makedirs(dir_work, exist_ok=True)
        
        ## non-analytical term
        if self.nac == 1:
            
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
                    nonanalytic=self.nac, borninfo=borninfo
                )
        
        ###
        self._prefix = inp['prefix']
         
        ## set kpoints and write a file
        filename = dir_work + '/' + propt + '.in'
        if propt == 'band':
            
            inp.set_kpoint(deltak=0.01)
            inp.to_file(filename=filename)
        
        elif propt == 'dos':
            
            inp.update({'kpts': kpts})
            
            inp.update({'pdos': 1})
            inp.to_file(filename=filename)
        
        elif propt == 'evec_commensurate':
   
            ### supercell matrix wrt primitive cell
            smat = np.dot(self.scell_matrix, np.linalg.inv(self.primitive_matrix))
            
            ### commensurate points
            from .structure.crystal import get_commensurate_points
            comm_pts = get_commensurate_points(smat)
            inp.update({'printevec': 1})
            inp.set_kpoint(kpoints=comm_pts)
            inp.to_file(filename=filename)
        
        elif propt == 'kappa':
            
            inp.update({'kpts': kpts})
            inp.update({'nac':  self.nac})
            inp.update({'isotope': 2})
            inp.update({'kappa_coherent': 1})
            inp.update({'restart': 0})
            inp.update({'tmin': 50})
            inp.update({'tmax': 1000})
            inp.update({'dt': 50})
            inp.to_file(filename=filename)
        
        elif propt == 'cv' or propt == 'lasso':
            
            """
            To Do
            ======

            * Adjust this part
            * set nbody and cutoffs
            
            """
            inp.update({'norder': 5})
            
            inp.update({'lmodel': 'enet'})
            inp.update({'maxiter': 1000})
            inp.update({'conv_tol': 1e-10})
            inp.update({'l1_ratio': 1.0})
            if propt == 'cv':
                inp.update({'cv': 5})
                inp.update({'cv_minalpha': 1e-8})
                inp.update({'cv_maxalpha': 1e-2})
                inp.update({'cv_nalpha': 50})
            elif propt == 'lasso':
                inp.update({'cv': 0})
                inp.update({'l1_alpha': None})      ### read l1_alpha
            
            inp.to_file(filename=filename)
            
            print("")
            print(" Output", filename)
            print("")
            print("HEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEERE")
            exit()

        else:
            print(" Error: %s is not supported." % propt)
            exit()
    
    def run_anphon(self, propt=None, force=False):
        """ Run anphon
        
        Args
        ---------
        propt : string
            "band", "dos", "evec_commensurate", or "kappa"

        force : bool
            If False, anphon will not be run if the same calculation had been
            conducted while, if True, anphon will be run forecely.
        
        """
        ## change directory
        dir_init = os.getcwd()
        if propt == 'band' or propt == 'dos':
            os.chdir(self.out_dirs['harm']['bandos'])
            
            ### store output file names
            if propt == 'band':
                ext = 'bands'
            else:
                ext = propt
         
        elif propt == 'kappa':
            if self.lasso == False:
                os.chdir(self.out_dirs['cube']['kappa'])
            else:
                os.chdir(self.out_dirs['lasso']['kappa'])
            
        elif propt == 'evec_commensurate':
            os.chdir(self.out_dirs['lasso']['evec'])
            
        else:
            print("ERROR")
            warnings.warn(" WARNNING: %s property is not supported." % (propt))
            exit()
        
        ### print title
        print("")
        msg = " Run anphon for %s " % propt
        border = "-" * (len(msg) + 2)
        print(msg)
        print(border)
        
        filename = "%s.in" % propt
        
        ### get number of threads
        #if nthreads is None:
        #    nthreads = os.cpu_count()
        
        ## prepare command and environment
        logfile = propt + '.log'
        cmd = "%s -n %d %s %s" %(
                self.command['mpirun'],
                self.command['nprocs'], 
                self.command['anphon'],
                filename)
        
        ## If the job has been finished, the same calculation is not conducted.
        ## The job status is judged from *.log file.
        if _anphon_finished(logfile) == False or force:
        
            os.environ['OMP_NUM_THREADS'] = str(self.command['nthreads'])
            
            ## run the job!!
            with open(logfile, 'w') as f:
                proc = subprocess.Popen(
                        cmd.split(), env=os.environ, stdout=f,
                        stderr=subprocess.PIPE)
                proc.wait()
            
        else:
            print(" %s has already been calculated." % propt)
        
        ## Return to the original directory
        os.environ.pop('OMP_NUM_THREADS', None)
        os.chdir(dir_init)

    def plot_bandos(self, **args):
        
        fn_band = self.out_dirs['harm']['bandos'] + '/' + self.prefix + '.bands'
        fn_dos = self.out_dirs['harm']['bandos'] + '/' + self.prefix + '.dos'
        
        options = figure_options()
        
        ### set figure name
        if 'figname' not in args.keys():
            options.figname = self.out_dirs['result'] + '/fig_bandos.png'
        else:
            options.figname = args['figname']
        
        from .plot.bandos import plot_bandos
        from .plot.alamode.band import Band
        from .plot.alamode.dos import Dos
        
        band = Band(file=fn_band)
        if os.path.exists(fn_dos):
            dos = Dos(file=fn_dos)
        else:
            dos = None
        
        ##
        options.prefix = self.prefix
        df = self.frequency_range[1] - self.frequency_range[0]
        options.y0 = self.frequency_range[0] - df * 0.05
        options.y1 = self.frequency_range[1] + df * 0.05
        
        ##
        for key in args.keys():
            setattr(options, key, args[key])
        
        ### output title
        print("")
        msg = " Plot band and DOS"
        border = "-" * (len(msg) + 2)
        print(msg)
        print(border)

        fig = plot_bandos(options.figname, band, dos, options)

    def plot_kappa(self):
        
        if self.lasso == False:
            fn_kappa = self.out_dirs['cube']['kappa'] + '/' + self.prefix + '.kl'
        else:
            fn_kappa = self.out_dirs['lasso']['kappa'] + '/' + self.prefix + '.kl'
        
        figname = self.out_dirs['result'] + '/fig_kappa.png'
        
        df = pd.DataFrame()
        data = np.genfromtxt(fn_kappa)
        
        nt = len(data)
        df['temperature'] = data[:,0]
        dirs = ['x', 'y', 'z']
        for i1 in range(3):
            d1 = dirs[i1]
            for i2 in range(3):
                d2 = dirs[i2]
                lab = 'k%s%s' % (d1, d2)
                num = i1*3 + i2 + 1
                df[lab] = data[:,num]
        kave = (df['kxx'].values + df['kyy'].values + df['kzz'].values) / 3.
        df['kave'] = kave
        
        print("")
        print(" ### Plot kappa")
        from .plot.kappa import plot_kappa
        plot_kappa(df, figname=figname)


class figure_options():
    def __init__(self):
        self.prefix = None
        self.figname = None
        self.plot_pdos = 1
        self.prefix2 = None
        self.pr_ratio = 0
        self.y0 = None
        self.y1 = None
        self.maxdos = None
        
        self.lw = 0.5
        self.yticks = None
        self.myticks = None

        self.fig_width = 3.0
        self.fig_aspect = 0.5
        self.dpi = 300
        self.col = 'blue'
        self.col2 = 'grey'
        self.wspace = 0.05
        self.unit = "cm"

        self.colorbar = 1
        self.legend = 1
        self.legend_los = 'best'

def _read_frequency_range(filename, format='anphon'):
    """ read minimum and maximum frequencies from .bands file created by anphon 
    """
    if format == 'anphon':
        data = np.genfromtxt(filename)
        fmax = np.max(data[:,-1])
        fmin = np.min(data[:,1])
        return fmin, fmax
    else:
        print(" Error: %s is not supported yet." % format)
        exit()

def _anphon_finished(logfile):
    try:
        lines = open(logfile, 'r').readlines()
        n = len(lines)
        for i in range(n):
            line = lines[n-1-i]
            data = line.split()
            if len(data) != 0:
                if 'Job finished' in line:
                    return True
                else:
                    return False
    except Exception:
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


def get_finite_displacements(structure, order, cutoffs=None, nbody=None, mag=None):
    """ Generate displacement patterns with ALM for FCs calculation
    Args
    ======
    cutoff : unit=[Bohr]
    mag : unit=[Ang]
    """
    nbody_new = []
    for i in range(order):
        nbody_new.append(nbody[i])
    
    ## get patterns with ALM
    patterns = run_alm(
            structure, order, cutoffs[:order], nbody_new,
            mode='suggest')
    
    if patterns is None:
        return None
    
    natoms = len(structure)
    
    ### convert patterns to vectors
    all_disps = []
    for pattern in patterns:
        all_disps.append(np.zeros((natoms, 3)))
        for each in pattern:
            if each[2].lower() == 'cartesian':
                all_disps[-1][each[0]] = each[1] * mag
            else:
                print(" Error during getting patterns.")
                exit()
    
    return all_disps

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

