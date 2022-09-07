"""
To Do
=============
1. Unknow error in mp-190. 
    The relaxation calculation was stopped.
    phonon dispersion was not calculated properly.
    
    cutoff : shape=(order, num_elems, num_elems)
        how to determine the order of elements???
        
        example for AlN:
        
        cutoff_radii = [np.ones((2, 2)) * -1, np.ones((2, 2)) * 4.8]
    
    nbody : shape=(maxorder)

    https://kitchingroup.cheme.cmu.edu/blog/2013/02/19/Subclassing-an-ase-calculators-vasp-calculator-to-do-post-analysis/

2. If the number of patterns is too large, use lasso regression with random
    displacement.
3. Use different supercell size for 3rd-order force constants

"""
# -*- coding: utf-8 -*-
import os
import numpy as np
import warnings
import logging

import ase, pymatgen
from ase.build import make_supercell
import subprocess
import pandas as pd

from . import output_directories, calculator_parameters, output_files
from .units import AToBohr, BohrToA
from .structure.format import change_structure_format
from .io.vasp import wasfinished, get_dfset, print_vasp_params
from .calculator import run_vasp

class AlmCalc():
    
    def __init__(self, prim_given, mpid='mp', 
            scell_matrix=None, primitive_matrix=None,
            scell_matrix3=None,
            cutoffs=None, nbody=[2,3], mag=0.01,
            cutoff2=-1, cutoff3=4.3,
            nac=None
            ):
        """
        Args
        =======
        prim_given : primitive structure
            different formats such as pymatgen and ASE are accepted.
        cutoffs : shape=(order, num_elems, num_elems), unit=[Ang]
            cutoff radii for force constants. If not given, it will be set
            automatically with cutoff2 and cutoff3. Default cutoff radii are 
            "-1" and "4.3" Angstrom for harmonic and cubid FCs, respectively.
        cutoff2, cutoff3 : float, unit=[Ang]
            cutoff radii.
            If cutoffs is None, cutoff2 and cutoff3 are used.
        """
        ### output directories
        self.out_dirs = {}
        for k1 in output_directories.keys():
            values1 = output_directories[k1]
            if type(values1) == str:
                self.out_dirs[k1] = './' + mpid + '/' + values1
            else:
                self.out_dirs[k1] = {}
                for k2 in values1.keys():
                    values2 = values1[k2]
                    self.out_dirs[k1][k2] = './' + mpid + '/' + values2
        
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
        if cutoffs is not None:
            self.cutoffs = cutoffs * AToBohr    ### 4.3 Ang == 8.0 Bohr
        else:
            n = self.num_elems
            self.cutoffs = np.asarray([
                np.ones((n,n)) * cutoff2 * AToBohr, 
                np.ones((n,n)) * cutoff3 * AToBohr
                ])
        ##
        self.nbody = nbody
        self.mag = mag            ### Ang
        
        self.force_calculators = {}
        
        self.nac = nac
        
        ###
        self._prefix = None      ## prefix for input files
        self._frequency_range = None
    
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
        
    def get_suggested_structures(self, order: None, random=False, npattern_random=0):
        """ Get structures suggested with ALM or a random displacement 
         
        order : int 
            1 (harmonic) or 2 (cubic)

        random : bool, default=False
            use the random displacement or not

        npattern_random : int
            Number of patterns generated with the random displacement
        
        """
        if order == 1:
            structure = self.supercell
        else:
            structure = self.supercell3
        
        ### get displacements
        if random == False:
            all_disps = get_displacements(
                    structure,
                    order=order, 
                    cutoffs=self.cutoffs[:order],  ## Bohr
                    nbody=self.nbody,
                    mag=self.mag                   ## Ang
                    )
        else:
            print("")
            print(" Random displacement method, which is not implemented yet.")
            print("")
            # AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
            #all_disps = get_random_displacements(
            #        structure,
            #        maxmag=self.mag      ## Ang
            #        )
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
    
    def calc_forces(self, order: None, calculator: None, 
            directory=None, force=False, 
            suggest_number_limit=1000, npattern_random=120
            ):
        """ Calculate forces for harmonic or cubic IFCs.
        VASP output will be stored in self.out_dirs['harm/cube']['force'].
        
        Args
        ======
        order : int
            1 (harmonic) or 2 (cubic)
        
        calculator : ase.calculators.vasp.Vasp
        
        directory : string
            If it's given, calculationts will be done in 
            {directory}/{prist, 1, 2, ...}.
        
        force : bool, [False]
            If it's True, calculations will be done forcelly even if they had 
            been already finished.

        suggest_number_limit : int
            If the number of suggested structures exceeds 
        
        npattern_random : int
            Number of suggested structures for the random displacement
        
        Return
        =======
        lasso : int 
            LASSO is not used (0) or is used (1).
        
        """
        if directory is not None:
            outdir0 = directory
        elif calculator.directory is not None:
            outdir0 = calculator.directory
        else:
            warnings.warning(" WARNING: Output directory is not defined.")
        
        print("")
        msg = " Force calculation (order: %s)" % order
        border = "=" * (len(msg) + 2)
        print(msg)
        print(border)

        ### get suggsted structures with ALM
        ### structures : dict of structure objects
        lasso = 0
        structures = self.get_suggested_structures(order)
        nsuggest =  len(structures)
        
        ###
        msg = "Number of the suggested structures with ALM : %d" % (nsuggest)
        print("\n", msg)
        
        ### LASSO is used when the number of the suggested structure is too large.
        if nsuggest > suggest_number_limit:
            lasso = 1
            msg = "Number of suggested structures, %d, "\
                    "exceeds the limit, %d." % (nsuggest, suggest_number_limit)
            print("\n", msg)
            structures = self.get_suggested_structures(
                    order, random=True, npattern_random=npattern_random
                    )
            
            exit()
        
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
            run_vasp(calculator, structure, method='custodian')
            print(" %s" % (outdir))
        
        print("")
        return lasso
    
    def calc_harmonic_force_constants(self, output=True, verbosity=0):
        """ calculate harmonic FCs and return ALM 
        Args
        =========
        output : bool
            If True, fcs.xml file is stored, if False, not.
        verbosity : int
            If 0, log is printed, if 1, not.
        """
        print("")
        print(" ### Calculate harmonic force constants")
        
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
        fc2_values, elem2_indices = run_alm(
                structure, 1, self.cutoffs[0], [self.nbody[0]], 
                mode='optimize',
                displacements=disps, forces=forces, 
                outfile=out_xml,
                verbosity=verbosity
                )
        
        return fc2_values, elem2_indices

    def calc_cubic_force_constants(self):
        """ calculate cubic IFCs and return ALM 
        """
        ### prepare directory and file names
        os.makedirs(self.out_dirs['result'], exist_ok=True)

        directory = self.out_dirs['cube']['force']
        
        out_dfset = self.out_dirs['result'] + '/' + output_files['cube_dfset']
        out_xml = self.out_dirs['result'] + '/' + output_files['cube_xml']
        
        ## parameters
        nbody = []
        for i in range(2):
            nbody.append(self.nbody[i])
        
        ## get dfset
        offset_xml = self.out_dirs['cube']['force'] + '/prist/vasprun.xml'
        disps, forces = get_dfset(
                directory, offset_xml=offset_xml, outfile=out_dfset)
        
        ## get fc2 and elem2
        fc2_values, elem2_indices = \
                self.calc_harmonic_force_constants(output=False)
        
        ## structure
        structure = ase.io.read(offset_xml)
        
        ## 3rd order model
        fcs_values, fcs_indices = run_alm(
                structure, 2, self.cutoffs[:2], nbody, mode='optimize',
                displacements=disps, forces=forces,
                fc2info=[fc2_values, elem2_indices],
                outfile=out_xml)
        
        return fcs_values, fcs_indices

    def write_anphon_input(self, propt='band', kpts=[15,15,15]):
        
        from .io.vasp import write_born_info
        from .io.alm import AnphonInput
        
        ## prepare filenames
        cwd = os.getcwd()
        if propt == 'band' or propt == 'dos':
            dir_work = cwd + '/' + self.out_dirs['harm']['bandos']
            fcsxml = '../../result/' + output_files['harm_xml']
            mode = 'phonons'
        elif propt == 'kappa':
            dir_work = cwd + '/' + self.out_dirs['cube']['kappa']
            fcsxml = '../../result/' + output_files['cube_xml']
            mode = 'RTA'
        
        ##
        born_xml = cwd + '/' + self.out_dirs['nac'] + '/vasprun.xml'
        
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
        
        ## property
        if propt == 'band':
            kpmode = 1
        elif propt == 'dos' or propt == 'kappa':
            kpmode = 2
        else:
            print(" Error: \"%s\" is not supported." % propt)
            exit()
        
        ## set input file for anphon
        anpinp = AnphonInput.from_structure(
                self.primitive,
                mode=mode,
                kpmode=kpmode,
                fcsxml=fcsxml,
                nonanalytic=self.nac, borninfo=borninfo
            )
        
        ### set primitive cell (dummy)
        anpinp.set_primitive(
                change_structure_format(
                    self.primitive, format="pymatgen-structure"
                    )
                )
        
        ###
        self._prefix = anpinp['prefix']

        ## set kpoints and write a file
        filename = dir_work + '/' + propt + '.in'
        if propt == 'band':
            anpinp.set_kpoint(deltak=0.01)
            anpinp.to_file(filename=filename)
        elif propt == 'dos':
            anpinp.update({'kpts': kpts})
            
            anpinp.update({'pdos': 1})
            anpinp.to_file(filename=filename)
        
        elif propt == 'kappa':
            anpinp.update({'kpts': kpts})
            anpinp.update({'nac':  self.nac})
            anpinp.update({'isotope': 2})
            anpinp.update({'kappa_coherent': 1})
            anpinp.update({'restart': 0})
            anpinp.update({'tmin': 50})
            anpinp.update({'tmax': 1000})
            anpinp.update({'dt': 50})
            anpinp.to_file(filename=filename)
        else:
            print(" Error: %s is not supported." % propt)
            exit()
    
    def run_anphon(self, propt='band', nprocs=1, nthreads=None, force=False):
        """ Run anphon
        Args
        ---------
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
            os.chdir(self.out_dirs['cube']['kappa'])
            
        else:
            print("ERROR")
            exit()
        
        ### print title
        print("")
        msg = " Run anphon with mode = %s " % propt
        border = "-" * (len(msg) + 2)
        print(msg)
        print(border)
        
        filename = "%s.in" % propt
        
        ## get number of threads
        if nthreads is None:
            nthreads = os.cpu_count()
        
        ## prepare command and environment
        logfile = propt + '.log'
        cmd = "%s -np %d anphon %s" %(
                calculator_parameters['mpirun'],
                nprocs, filename)
        
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
            
        else:
            print(" %s has already been calculated." % propt)
        
        ## Return to the original directory
        os.environ.pop('OMP_NUM_THREADS', None)
        os.chdir(dir_init)

    def plot_bandos(self, **args):
        
        cwd = os.getcwd()
        fn_band = (cwd + '/' + self.out_dirs['harm']['bandos'] + '/' + 
                self.prefix + '.bands')
        fn_dos = (cwd + '/' + self.out_dirs['harm']['bandos'] + '/' + 
                self.prefix + '.dos')
        
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
        
        fn_kappa = (os.getcwd() + '/' + 
                self.out_dirs['cube']['kappa'] + '/' +
                self.prefix + '.kl')
        
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

#def calc_force_constants(alm,
#        displacements=None, forces=None,
#        order=None, cutoff=None, nbody=None, outfile=None):
#    """ Calculate force constants
#    Args
#    ======
#    order : int
#    cutoff : float, unit=Bohr
#    nbody : shape=(order)
#    outfile : string
#    """
#    #alm.define(order, cutoff_radii=cutoff, nbody=nbody,
#    #        symmetrization_basis="Cartesian")
#    #alm.displacements = displacements
#    #alm.forces = forces
#    #info = alm.optimize()
#    #if info == 1:
#    #    warnings.warn(" Fitting with ALM was not successful.")
#    if outfile is not None:
#        alm.save_fc(outfile, format='alamode')
#        print(" Output", outfile)
#    return alm

def run_alm(structure, order, cutoffs, nbody, mode=None,
        displacements=None, forces=None, outfile=None,
        fc2info=None, verbosity=0
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
    
    Return
    ===========
    If mode == 'suggest':
    patterns : shape=(nstruct, natoms, 3)

    If mode == 'optimize':
    alm.get_fc(order), which are fc*_values, elem*_indices
    
    """
    from alm import ALM

    if type(structure) == "ase.atoms.Atoms":
        atoms = structure.copy()
    else:
        atoms = change_structure_format(structure, format='ase')
    
    lave = atoms.cell * AToBohr
    xcoord = atoms.get_scaled_positions()
    kd = atoms.get_atomic_numbers()
    
    with ALM(lave, xcoord, kd) as alm:
        
        alm.set_verbosity(verbosity)
        alm.define(order, cutoff_radii=cutoffs, nbody=nbody,
                symmetrization_basis="Cartesian")
        
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
            
            ## freeze fc2
            if fc2info is not None:
                alm.freeze_fc(fc2info[0], fc2info[1])
            
            info = alm.optimize()
            if outfile is not None:
                alm.save_fc(outfile, format='alamode')
                print(" Output", outfile)
            return alm.get_fc(order)
        
        else:
            warnings.warn(" WARNING: mode %s is nto supported" % (mode))

def get_displacements(structure, order, cutoffs=None, nbody=None, mag=None):
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

