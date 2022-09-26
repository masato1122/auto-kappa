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

import ase, pymatgen
from ase.build import make_supercell
import subprocess
import shutil
import pandas as pd

from .. import output_directories, output_files
from ..units import AToBohr, BohrToA
from ..structure.crystal import change_structure_format, get_formula
from ..io.vasp import wasfinished, get_dfset, read_dfset, print_vasp_params
from ..calculator import run_vasp

from ..io.vasp import write_born_info
from ..io.alm import AlmInput, AnphonInput

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
    
    def __init__(self, prim_given, 
            material_name='mp',
            scell_matrix=None, primitive_matrix=None,
            scell_matrix3=None,
            restart=1,
            cutoffs=None, nbody=[2,3,3,2,2], mag=0.01,
            cutoff2=-1, cutoff3=4.3, 
            order_lasso=5,
            nac=None,
            ncores=1,
            verbosity=0,
            command={
                'mpirun': 'mpirun', "nprocs": 1, 
                "nthreads": 1, "anphon": "anphon", "alm": "alm",
                },
            ):
        """ Calculations with ALM and Anphon are managed with this class.
         
        Args
        =======
        prim_given : primitive structure
            Different formats such as pymatgen and ASE are accepted while the
            format is changed to ase-Atoms in this module.

        material_name : string
            Material ID, which is used just to determine the directory name.
            It can be anything such as the composition.
        
        restart : int, default 1
            The calculation will restart (1) or NOT restart (0) when the
            directory exists.

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
                "nthreads": 1, "anphon": "anphon", "alm": "alm",
                }
        
        Example
        ----------
        >>> almcalc = AlamodeCalc(
        >>>     primitive,
        >>>     material_name='mpid-149',
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
        >>> almcalc.run_alamode(propt='dos')
        """
        ### set name of base directory
        self.set_material_name(material_name, restart)

        ### output directories
        dir_init = os.getcwd()
        self.out_dirs = {}
        for k1 in output_directories.keys():
            values1 = output_directories[k1]
            if type(values1) == str:
                self.out_dirs[k1] = dir_init+'/'+material_name+'/'+values1
            else:
                self.out_dirs[k1] = {}
                for k2 in values1.keys():
                    values2 = values1[k2]
                    self.out_dirs[k1][k2] = dir_init+'/'+material_name+'/' + values2
        
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
        self.mag = mag            ### Ang
        
        self.force_calculators = {}
        
        self.nac = nac

        self.lasso = False

        self._nmax_suggest = None

        self._command = command

        self._ncores = ncores

        self._order_lasso = order_lasso
        
        ###
        self._prefix = get_formula(self.primitive)  ## prefix for input files
        self._frequency_range = None
        
        self.verbosity = verbosity
    
    @property
    def material_name(self):
        return self._material_name

    def set_material_name(self, dir0, restart):
        """
        If restart is off, the directory name needs to be determined charefully.
        """
        if restart == 1:
            self._material_name = dir0
        else:
            if os.path.exists(dir0):
                for i in range(2, 100):
                    dir1 = dir0 + '_%d' % i
                    if os.path.exists(dir1) == False:
                        self._material_name = dir1
                        break
    @property
    def ncores(self):
        return self._ncores
    
    @property
    def prefix(self):
        if self._prefix is None:
            self._prefix = get_formula(self.primitive)
        return self._prefix

    @property
    def order_lasso(self):
        return self._order_lasso

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
    def frequency_range(self):
        if self._frequency_range is None:
            fn = self.out_dirs['harm']['bandos'] + '/' + self.prefix + '.bands'
            if os.path.exists(fn) == False:
                print(" Phonon dispersion needs to be calculated.")
            else:
                fmin, fmax = _read_frequency_range(fn, format='anphon')
                self._frequency_range = [fmin, fmax]
        return self._frequency_range
        
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
            propt = 'evec_commensurate'
            self.write_alamode_input(propt=propt)
            self.run_alamode(propt=propt)

            ###
            fevec = self.out_dirs['lasso']['evec'] + '/' + self.prefix + '.evec'
            all_disps = self._get_displacements(
                    "random_normalcoordinate",
                    file_evec=fevec,
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
    
    #def _get_fd_displacements(self, filename: None):
    #    
    #    from .tools.VASP import VaspParser
    #    from .tools.GenDisplacement import AlamodeDisplace
    #    
    #    codeobj = VaspParser()
    #    codeobj.set_initial_structure(self.supercell)
    #    
    #    almdisp = AlamodeDisplace(
    #            'fd', codeobj,
    #            #file_primitive=file_prim,
    #            primitive=self.primitive,
    #            verbosity=self.verbosity
    #            )
    #    if almdisp is None:
    #        return None
    #    
    #    print(filename)
    #    header_list, disp_list = almdisp.generate(file_pattern=filename,
    #            magnitude=0.1)
    #    
    #    print("CCC")
    #    exit()
    #    
    #    all_disps = np.zeros_like(disp_list)
    #    for i, each in enumerate(disp_list):
    #        all_disps[i] = np.dot(each, self.supercell.cell)
    #    return all_disps 
    
    def _get_displacements(
            self, displacement_mode: None,
            file_pattern=None, 
            file_evec=None, number_of_displacements=1, temperature=500., classical=False,
            ):
        from .tools.VASP import VaspParser
        from .tools.GenDisplacement import AlamodeDisplace
        
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
                displacement_mode, codeobj,
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
        msg += "\n"
        msg += " %d patterns will be generated.\n" % (number_of_displacements)
        msg += "\n"
        print(msg)
        
        if displacement_mode == 'fd':
            
            if file_pattern is None:
                warnings.warn(" Error: file_pattern must be given.")
                exit()

            header_list, disp_list = almdisp.generate(
                    file_pattern=[file_pattern], magnitude=self.mag,
                    )
            
        elif displacement_mode == 'random_normalcoordinate':
            
            header_list, disp_list = almdisp.generate(
                    temperature=temperature,
                    number_of_displacements=number_of_displacements,
                    classical=classical
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
            nmax_suggest=200, frac_nrandom=0.02, 
            temperature=500., classical=False,
            calculate_forces=False, output_dfset=1
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
            line = " Force calculation (order: %s)" % order
        else:
            line = " Generate displacement structures (order: %s)" % order
        
        print("")
        msg = " " + line
        border = "=" * (len(msg) + 2)
        print(msg)
        print(border)

        self._nmax_suggest = nmax_suggest
        
        ### get suggsted structures with ALM
        ### structures : dict of structure objects
        #
        ### ver.1 with ALM library
        #structures = self.get_suggested_structures(order, disp_mode='fd')
        #nsuggest =  len(structures)
        #
        ### ver. 2 from alamode output file
        nsuggest = self._get_number_of_suggested_structures(order)
        
        ###
        msg = "\nNumber of the suggested structures with ALM : %d\n" % (nsuggest)
        print(msg)
        
        ### output directory
        if order == 1:
            outdir0 = self.out_dirs['harm']['force']
        elif order == 2:
            if nsuggest <= nmax_suggest:
                outdir0 = self.out_dirs['cube']['force']
            else:
                outdir0 = self.out_dirs['lasso']['force']
        else:
            warnings.warning(" WARNING: given order is not supported.")
        
        ### Compressive sensing, LASSO, is used when the number of the 
        ### suggested structure (nsuggest) excees nmax_suggest.
        ### When nsuggest > nmax_suggest, max(nsuggest * frac_nrandom,
        ### nmax_suggest) structures will be generated with normal coordinates.
        if order != 1 and nsuggest > nmax_suggest:

            self.lasso = True

            nrandom = int(frac_nrandom * nsuggest + 0.5)
            ngenerated = max(nmax_suggest, nrandom)
             
            msg = "\n"
            msg += " Maximum limit of the number of suggested patterns : %d\n" % (nmax_suggest)
            msg += "\n"
            msg += " ### Number of suggested patterns exceeds the maximum limit.\n"
            msg += "\n"
            #msg += " ### Number of patterns ###\n"
            msg += " Fractional number of the random patterns : %.3f\n" % (frac_nrandom)
            #msg += " * random    : %d\n" % (nrandom)
            #msg += " * generated : %d\n" % (ngenerated)
            print(msg)
            
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
                df_prefix = 'harm'
            elif order > 1:
                if self.lasso:
                    df_prefix = 'lasso'
                else:
                    df_prefix = 'cube'
            
            os.makedirs(self.out_dirs['result'], exist_ok=True)
            offset_xml = outdir0 + '/prist/vasprun.xml'
            outfile = self.out_dirs['result'] + "/DFSET." + df_prefix
            if os.path.exists(outfile) == False or output_dfset == 2:
                get_dfset(outdir0, offset_xml=offset_xml, outfile=outfile)
            else:
                msg = "\n"
                msg += " %s already exists\n" % (outfile)
                print(msg)
        print("")
    
    #def calc_harmonic_force_constants_with_ALM(self, output=True):
    #    """ calculate harmonic FCs with python ALM library and return ALM 
    #    ALM library may incompatible with alm command when lmodel = enet for
    #    some cases. Therefore, ALM library is not used in auto-kappa. 
    #    
    #    Args
    #    -----
    #    output : bool
    #        If True, fcs.xml file is stored, if False, not.
    #    
    #    """
    #    msg = "\n"
    #    msg += " ### Calculate harmonic force constants\n"
    #    msg += "\n"
    #
    #    ### prepare directory and file names
    #    os.makedirs(self.out_dirs['result'], exist_ok=True)
    #    
    #    directory = self.out_dirs['harm']['force']
    #    
    #    out_dfset = self.out_dirs['result'] + '/' + output_files['harm_dfset']
    #    if output:
    #        out_xml = self.out_dirs['result'] + '/' + output_files['harm_xml']
    #    else:
    #        out_xml = None
    #    
    #    ## get dfset
    #    offset_xml = self.out_dirs['harm']['force'] + '/prist/vasprun.xml'
    #    
    #    if os.path.exists(out_dfset) == False:
    #        
    #        disps, forces = get_dfset(
    #                directory, offset_xml=offset_xml, outfile=out_dfset)
    #    
    #    else:
    #        
    #        disps, forces = read_dfset(out_dfset, natoms=len(self.supercell))
    #        
    #    ## get structure
    #    structure = ase.io.read(offset_xml)
    #    
    #    ## calculate force constants
    #    ## out = [fc2_values, elem2_indices]
    #    os.environ['OMP_NUM_THREADS'] = str(self.ncores)
    #    out = run_alm(
    #            structure, 1, self.cutoffs[0], [self.nbody[0]], 
    #            mode='optimize',
    #            displacements=disps, forces=forces, 
    #            outfile=out_xml,
    #            verbosity=self.verbosity
    #            )
    #    os.environ.pop('OMP_NUM_THREADS', None)
    #    return out
    
    
    #def calc_anharm_force_constants(self):
    #    """ calculate cubic IFCs and return ALM 
    #    
    #    Args
    #    ----
    #    maxorder : int, default 5
    #        order for lasso calculation
    #    """
    #    ### prepare directory and file names
    #    os.makedirs(self.out_dirs['result'], exist_ok=True)
    #    
    #    mode = 'optimize'
    #    
    #    if self.lasso == False:
    #
    #        order = 2
    #        
    #        directory = self.out_dirs['cube']['force']
    #        out_dfset = self.out_dirs['result'] + '/' + output_files['cube_dfset']
    #        out_xml = self.out_dirs['result'] + '/' + output_files['cube_xml']
    #        offset_xml = self.out_dirs['cube']['force'] + '/prist/vasprun.xml'
    #        
    #    else:
    #        
    #        order = self.order_lasso
    #        
    #        directory = self.out_dirs['lasso']['force']
    #        out_dfset = self.out_dirs['result'] + '/' + output_files['lasso_dfset']
    #        out_xml = self.out_dirs['result'] + '/' + output_files['lasso_xml']
    #        offset_xml = self.out_dirs['lasso']['force'] + '/prist/vasprun.xml'
    #        
    #    ### print    
    #    msg = "\n"
    #    msg += " ### Calculate anharmonic FCs (order=%d)\n" % (order)
    #    msg += "\n"
    #    print(msg)
    #    
    #    ### set nbody
    #    if len(self.nbody) < order:
    #        self.set_nbody_automatically(order=order)
    #    
    #    nbody = []
    #    for i in range(order):
    #        nbody.append(self.nbody[i])
    #
    #    ### get dfset
    #    if os.path.exists(out_dfset):
    #        
    #        disps, forces = read_dfset(
    #                out_dfset,
    #                natoms=len(self.supercell3)
    #                )
    #    
    #    else:
    #        
    #        disps, forces = get_dfset(
    #                directory, 
    #                offset_xml=offset_xml, 
    #                outfile=out_dfset
    #                )
    #    
    #    ## For LASSO, extract only required data
    #    ## This part is required when too many patterns were calculated, but
    #    ## smaller data are needed for LASSO.
    #    if self.lasso and self.nmax_suggest is not None:
    #        if len(disps) > self.nmax_suggest:
    #            disps = disps[:self.nmax_suggest]
    #            forces = forces[:self.nmax_suggest]
    #    
    #    ## get fc2 and elem2
    #    fc2_values, elem2_indices = \
    #            self.calc_harmonic_force_constants(output=False)
    #    
    #    ## structure
    #    structure = ase.io.read(offset_xml)
    #    
    #    ## 3rd order model
    #    ## out = [fcs_values, fcs_indices]
    #    os.environ['OMP_NUM_THREADS'] = str(self.ncores)
    #    out = run_alm(
    #            structure, order, self.cutoffs[:order], nbody, mode='optimize',
    #            displacements=disps, forces=forces,
    #            fc2info=[fc2_values, elem2_indices],
    #            outfile=out_xml, lasso=self.lasso,
    #            verbosity=self.verbosity)
    #    os.environ.pop('OMP_NUM_THREADS', None)
    #    return out
    
    def write_alamode_input(self, propt=None, kpts=[15,15,15], **kwargs):
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

        ### prepare filenames
        if propt == 'band' or propt == 'dos':
            
            dir_work = self.out_dirs['harm']['bandos']
            fcsxml = '../../result/' + output_files['harm_xml']
            
        elif propt == 'evec_commensurate':
            
            dir_work = self.out_dirs['lasso']['evec']
            fcsxml = '../../result/' + output_files['harm_xml']
            
        elif propt == 'kappa':
            
            if self.lasso == False:
                dir_work = self.out_dirs['cube']['kappa']
                fcsxml = '../../result/' + output_files['cube_xml']
            else:
                dir_work = self.out_dirs['lasso']['kappa']
                fcsxml = '../../result/' + output_files['lasso_xml']
            
        elif propt == 'lasso' or propt == 'cv':
            
            dir_work = self.out_dirs['lasso'][propt]
            dfset = '../../result/' + output_files['lasso_dfset']
            fc2xml = '../../result/' + output_files['harm_xml']
        
        elif propt == 'fc2':
            
            dir_work = self.out_dirs['harm']['force']
            dfset = '../../result/' + output_files['harm_dfset']
            fc2xml = None
        
        elif propt == 'fc3':
            
            dir_work = self.out_dirs['cube']['force']
            dfset = '../../result/' + output_files['cube_dfset']
            fc2xml = '../../result/' + output_files['harm_xml']
            
        elif propt == 'suggest':
            
            if 'order' in kwargs:
                order = kwargs['order']
            
            else:
                print("")
                print(" Order is not given. Set order = 1.")
                order = 1
            
            dfset = None
            fc2xml = None
            
            if order == 1:
                dir_work = self.out_dirs['harm']['suggest']
            elif order == 2:
                dir_work = self.out_dirs['cube']['suggest']
            
        else:
            print("")
            warnings.warn(" Error: %s is not supported yet." % (propt))
            exit()

        ##
        born_xml = self.out_dirs['nac'] + '/vasprun.xml'
        
        ## prepare directory
        os.makedirs(dir_work, exist_ok=True)
        
        ## non-analytical term
        if self.nac == 1 and alamode_type == 'anphon':
            
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
        
        elif propt == 'dos':
            
            if self.frequency_range is not None:
                df = self.frequency_range[1] - self.frequency_range[0]
                fmin = self.frequency_range[0] - df * 0.05
                fmax = self.frequency_range[1] + df * 0.05
            
            inp.update({'emin': fmin})
            inp.update({'emax': fmax})
            inp.update({'kpts': kpts})
            inp.update({'pdos': 1})
        
        elif propt == 'evec_commensurate':
   
            ### supercell matrix wrt primitive cell
            smat = np.dot(self.scell_matrix, np.linalg.inv(self.primitive_matrix))
            
            ### commensurate points
            from ..structure.crystal import get_commensurate_points
            comm_pts = get_commensurate_points(smat)
            inp.update({'printevec': 1})
            inp.set_kpoint(kpoints=comm_pts)
        
        elif propt == 'kappa':
            
            inp.update({'kpts': kpts})
            inp.update({'nac':  self.nac})
            inp.update({'isotope': 2})
            inp.update({'kappa_coherent': 1})
            inp.update({'restart': 0})
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
            if propt == 'cv' or propt == 'lasso':
                order = self.order_lasso
            elif propt == 'fc2':
                order = 1
            elif propt == 'fc3':
                order = 2
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
                inp.update({'maxiter': 2000})    ## smaller than the default, 10000
                inp.update({'conv_tol': 1e-10})  ## strincter than the default, 1e-8
                if propt == 'cv':
                    inp.update({'cv': 5})
                    ### The code set automatically.
                    #inp.update({'cv_maxalpha': 1e-2})
                    #inp.update({'cv_minalpha': 1e-8})
                    inp.update({'cv_nalpha': 50})
                elif propt == 'lasso':
                    
                    alpha = self.get_suggested_l1alpha()
                    
                    inp.update({'cv': 0})
                    inp.update({'l1_alpha': alpha})      ### read l1_alpha
        
        else:
            warnings.warn(" Error: %s is not supported." % propt)
            exit()
        
        inp.update(kwargs)
        inp.to_file(filename=filename)
    
    def get_suggested_l1alpha(self):

        fn = self.out_dirs['lasso']['cv']+'/'+self.prefix+'.cvscore'
        try:
            lines = open(fn, 'r').readlines()
            for ll in lines:
                if "Minimum CVSCORE" in ll:
                    alpha = float(ll.split()[-1])
                    return alpha
            return None
        except Exception:
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
     
    def run_alamode(self, propt=None, force=False, **args):
        """ Run anphon
        
        Args
        ---------
        propt : string
            "band", "dos", "evec_commensurate", or "kappa"

        force : bool
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
            if self.lasso == False:
                workdir = self.out_dirs['cube']['kappa']
            else:
                workdir = self.out_dirs['lasso']['kappa']
            
        elif propt == 'evec_commensurate':
            
            workdir = self.out_dirs['lasso']['evec']
        
        elif propt == 'cv' or propt == 'lasso':
            
            workdir = self.out_dirs['lasso'][propt]
        
        elif propt == 'fc2':

            workdir = self.out_dirs['harm']['force']
        
        elif propt == 'fc3':

            workdir = self.out_dirs['cube']['force']
        
        elif propt == 'suggest':

            if 'order' in args:
                order = args['order']
            else:
                order = 1

            if order == 1:
                workdir = self.out_dirs['harm']['suggest']
            elif order == 2:
                workdir = self.out_dirs['cube']['suggest']
            else:
                print("")
                warnings.warn(" Error: order must be gien properly.")
                exit()

        else:
            print("")
            warnings.warn(" WARNNING: %s property is not supported." % (propt))
            exit()
        
        ### print title
        msg = "\n"
        msg += " Run %s for %s\n" % (alamode_type, propt)
        msg += "-" * (len(msg) + 2) + "\n"
        print(msg)
            
        filename = "%s.in" % propt
        logfile = propt + '.log'
        
        ## prepare command and environment
        val = run_alamode(filename, logfile, workdir=workdir, force=force,
                mpirun=self.command['mpirun'], nprocs=self.command['nprocs'],
                nthreads=self.command['nthreads'],
                command=self.command[alamode_type],
                )
        if val == 2:
            print(" %s has been already calculated." % propt)
        
        
        if mode == 'optimize' and propt in ['lasso', 'fc2', 'fc3']:
            
            ### copy the generated FCs file to "result" directory
            self._copy_generated_fcsfiles(propt=propt)
            
            if propt in ['lasso']:
                self._plot_cvsets()
        
    def _copy_generated_fcsfiles(self, propt=None):
        """ Copyr a FCs file into the "result" directory.
        """
        ext = 'xml'
        
        if propt == 'lasso':
            fn1 = self.out_dirs['lasso']['lasso']+'/'+self.prefix+'.'+ext
            fn2 = self.out_dirs['result']+'/FCs_lasso.'+ext
        
        elif propt == 'fc2':
            fn1 = self.out_dirs['harm']['force']+'/'+self.prefix+'.'+ext
            fn2 = self.out_dirs['result']+'/FCs_harm.'+ext
        
        elif propt == 'fc3':
            fn1 = self.out_dirs['cube']['force']+'/'+self.prefix+'.'+ext
            fn2 = self.out_dirs['result']+'/FCs_cube.'+ext
        
        else:
            print("")
            warnings.warn(" Error: %s is not supported yet." % propt)
            exit()
        
        if os.path.exists(fn1) == False:
            warnings.warn(' %s does not exist.' % fn1)
        else:
            if os.path.exists(fn2):
                msg = "\n"
                msg += " %s was overwritten.\n\n" % (fn2)
            else:
                msg = "\n"
                msg += " %s was created.\n\n" % fn2
            print(msg)
            shutil.copy(fn1, fn2)
    
    def _plot_cvsets(self):
        from ..plot.lasso import plot_cvsets
        figname = self.out_dirs['result'] + '/fig_cvsets.png'
        print("")
        print(" ### Plot CV results ###")
        plot_cvsets(
                directory=self.out_dirs['lasso']['cv'],
                figname=figname
                )
        print("")

    def plot_bandos(self, **args):
        
        ### set figure name
        if 'figname' not in args.keys():
            figname = self.out_dirs['result'] + '/fig_bandos.png'
        else:
            figname = args['figname']
        
        from ..plot.bandos import plot_bandos
        
        ### output title
        print("")
        msg = " Plot band and DOS"
        border = "-" * (len(msg) + 2)
        print(msg)
        print(border)
        
        fig = plot_bandos(
                directory=self.out_dirs['harm']['bandos'], 
                prefix=self.prefix, 
                figname=figname,
                **args
                )
    
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
        from ..plot.kappa import plot_kappa
        plot_kappa(df, figname=figname)


def run_alamode(filename, logfile, workdir='.', force=False,
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
    if _anphon_finished(logfile) == False or force:
    
        os.environ['OMP_NUM_THREADS'] = str(nthreads)
        
        ## run the job!!
        with open(logfile, 'w') as f:
            proc = subprocess.Popen(
                    cmd.split(), env=os.environ, stdout=f,
                    stderr=subprocess.PIPE)
            proc.wait()
        val = 1

    else:
        val = 2
    
    ## Return to the original directory
    os.environ.pop('OMP_NUM_THREADS', None)
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


#def get_finite_displacements(structure, order, cutoffs=None, nbody=None, mag=None):
#    """ Generate displacement patterns with ALM for FCs calculation
#    Args
#    ======
#    cutoff : unit=[Bohr]
#    mag : unit=[Ang]
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
#                all_disps[-1][each[0]] = each[1] * mag
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

