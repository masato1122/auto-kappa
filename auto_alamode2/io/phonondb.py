# -*- coding: utf-8 -*-
import os.path
import numpy as np
import warnings
import glob

import ase.io
from phonopy import Phonopy
from ase.calculators.vasp import Vasp
from .. import output_directories
from ..structure.format import change_structure_format
from .vasp import read_incar, read_poscar, read_kpoints, wasfinished

class Phonondb():
    """ Read files in phonondb
    Args
    =========
    mode : string
        "relax", "force", "nac", ...
    """
    def __init__(self, directory, mpid='mp', nac=None):
        
        ###
        #self.out_dirs = {}
        #for k1 in output_directories.keys():
        #    values1 = output_directories[k1]
        #    if type(values1) == str:
        #        self.out_dirs[k1] = './' + mpid + '/' + values1
        #    else:
        #        self.out_dirs[k1] = {}
        #        for k2 in values1.keys():
        #            values2 = values1[k2]
        #            self.out_dirs[k1][k2] = './' + mpid + '/' + values2

        ###
        self.directory = directory
        self._filenames = None
        
        self._parameters = None
        
        self._scell_matrix = None
        self._primitive_matrix =None
        
        self._phonon = None
        self._unitcell = None
        
        ## relaxed structure
        self._primitive = None
        self._supercell = None
         
        self._relaxed_primitive = None
        
        ##
        self._nac = None

    #@property
    #def calculator_params(self):
    #    return self._calculator_params
    #    
    #def set_calculator_params(self, keys):
    #    self._calculator_params.update(keys)
    
    @property
    def nac(self):
        if self._nac is None:
            if 'INCAR-nac' in list(self.filenames_as_dict.keys()):
                self._nac = 1
            else:
                self._nac = 0
        return self._nac
    
    @property
    def filenames(self):
        if self._filenames is None:
            line = self.directory + '/*'
            self._filenames = glob.glob(line)
        return self._filenames
    
    @property
    def filenames_as_dict(self):
        fn_dict = {}
        for fn_each in self.filenames:
            key = fn_each.split('/')[-1]
            fn_dict[key] = fn_each
        return fn_dict
        
    @property
    def scell_matrix(self):
        if self._scell_matrix is None:
            sc, pm = read_phonopy_conf(self.directory + '/phonopy.conf')
            self._scell_matrix = sc
            self._primitive_matrix = pm
        return self._scell_matrix
    
    @property
    def primitive_matrix(self):
        if self._primitive_matrix is None:
            sc, pm = read_phonopy_conf(self.directory + '/phonopy.conf')
            self._scell_matrix = sc
            self._primitive_matrix = pm
        return self._primitive_matrix
    
    @property
    def unitcell(self):
        if self._unitcell is None:
            fn = self.directory + '/POSCAR-unitcell'
            self._unitcell = read_poscar(fn)
        return self._unitcell
    
    def get_unitcell(self, format='phonopy'):
        if self._unitcell is None:
            fn = self.directory + '/POSCAR-unitcell'
            self._unitcell = read_poscar(fn)
        return change_structure_format(self._unitcell, format=format)
    
    @property
    def phonon(self):
        if self._phonon is None:
            self._phonon = Phonopy(
                    self.unitcell, self.scell_matrix, 
                    primitive_matrix=self.primitive_matrix
                    )
        return self._phonon
    
    @property
    def primitive(self):
        self._phonon = self.phonon
        return self._phonon
    
    def get_primitive(self, format='phonopy'):
        self._phonon = self.phonon
        return change_structure_format(self._phonon.primitive, format=format)   

    @property
    def supercell(self):
        if self._phonon is None:
            self._phonon = self.phonon
        return self._phonon.supercell
    
    def get_supercell(self, format='phonopy'):
        if self._phonon is None:
            self._phonon = self.phonon
        return change_structure_format(self._phonon.supercell, format=format)
    
    #@property
    #def relaxed_primitive(self):
    #    if self._relaxed_primitive is None:
    #        fn = self.out_dirs['relax'] + '/vasprun.xml'
    #        self._relaxed_primitive = self._read_vasprun_xml(fn)
    #    return self._relaxed_primitive
    
    def get_kpoints(self, mode=None):
        """
        Return
        -----------
        kpoints : Kpoints obj for the given mode
        """
        key = 'KPOINTS-' + mode
        if key in list(self.filenames_as_dict.keys()):
            fn = self.filenames_as_dict[key]
            kpoints = read_kpoints(fn)
            return kpoints
        else:
            return None
    
    def _read_vasprun_xml(self, filename):
        """ 
        filename : vasprun.xml filename
        """
        return ase.io.read(filename, format='vasp-xml', index=-1)
    
    #def run_relaxation(self, nprocs=1,
    #        mpirun=calculator_parameters['mpirun'],
    #        vasp=calculator_parameters['vasp'],
    #        params={'ediffg': -0.0001,}
    #        ):
    #    """ Run structure optimization and return contents of vasprun.xml
    #    """
    #    from ase.calculators.vasp.create_input import GenerateVaspInput
    #    
    #    ### Check if the job has been done or not.
    #    print("")
    #    print(" ### Structure relaxation")
    #    
    #    outdir = self.out_dirs['relax']
    #    if wasfinished(outdir, filename='vasprun.xml') == False:
    #        
    #        print(" calculating...")
    #        calc = self.get_calculator(mode='relax')
    #        calc.command = "%s -n %d %s" % (mpirun, nprocs, vasp)
    #        calc.directory = self.out_dirs['relax']
    #        
    #        ######################################
    #        ### update parameters ################
    #        ######################################
    #        GenerateVaspInput.set(calc, **params)
    #        calc.results.clear()
    #        
    #        ###
    #        prim = change_structure_format(self.primitive, format='ase')
    #        prim.calc = calc
    #        ene = prim.get_potential_energy()
    #    else:
    #        print("")
    #        print(" Relaxation calculation had been already done.")
    #        xml = outdir + '/vasprun.xml'
    #        atoms = ase.io.read(xml, format='vasp-xml')
    #        ene = atoms.get_potential_energy()
    #    
    #    print(" Energy: %.3f eV" % (ene))
    #  
    #def calculate_born_effective_charge(self, nprocs=1,
    #        mpirun=calculator_parameters['mpirun'],
    #        vasp=calculator_parameters['vasp']
    #        ):
    #    """ Calculate Born effective charge
    #    """
    #    print("")
    #    print(" ### Born effective charge")
    #    ### Check if the job has been done or not.
    #    outdir = self.out_dirs['nac']
    #    if wasfinished(outdir, filename='vasprun.xml') == False:
    #        
    #        calc = self.get_calculator(mode='nac')
    #        if calc is None:
    #            print("")
    #            print(" Born effective charge is not calculated.")
    #            print("")
    #            return None
    #         
    #        print(" calculating...")
    #         
    #        calc.command = "%s -n %d %s" % (mpirun, nprocs, vasp)
    #        calc.directory = self.out_dirs['nac']
    #        
    #        prim = change_structure_format(self.relaxed_primitive, format='ase')
    #        prim.calc = calc
    #        prim.get_potential_energy()
    #    
    #    else:
    #        print(" Born effective change had already been calculated.")
    #    print("")
    # 
    #def get_calculator(self, mode=None):
    #    """ Calculator with ASE
    #    mode : string
    #        'relax', 'force', or 'nac'
    #    """
    #    calc = Vasp(
    #            setups=calculator_parameters['setups'],
    #            xc=calculator_parameters['xc'])
    #    
    #    ## get kpoints
    #    key = 'KPOINTS-' + mode
    #    if key in list(self.filenames_as_dict.keys()):
    #        fn_kpts = self.filenames_as_dict[key]
    #        calc.read_kpoints(fn_kpts)
    #    else:
    #        print(" Cannot find %s" % key)
    #        #print(" kpoints needs to be given automatically.")
    #        #from ..calculator import get_recommended_kpoints
    #        #if mode.lower() == 'relax' or mode.lower() == 'nac':
    #        #    structure = self.primitive
    #        #elif mode.lower() == 'force':
    #        #    structure = self.supercell
    #        return None
    #     
    #    ## get incar
    #    key = 'INCAR-' + mode
    #    if key in list(self.filenames_as_dict.keys()):
    #        fn_incar = self.filenames_as_dict[key]
    #        calc.read_incar(fn_incar)
    #    else:
    #        print(" Cannot find %s" % key)
    #        return None
    #    
    #    return calc

    def _check_file(self, filename):
        if os.path.exists(filename) == False:
            print(" %s does not exist." % filename)

def read_phonopy_conf(filename):
    """ Read phonopy.conf in phonondb 
    Args
    ========
    filename : string, phonopy.conf

    Return
    =======
    sc_mat : shape=(3,3), int
    prim_mat : shape=(3,3), flaot
        default: identity matrix
    """
    from phonopy.cui.settings import PhonopyConfParser
    confparser = PhonopyConfParser(filename=filename)
    params = confparser.get_configures()
    
    from fractions import Fraction
    matrix = np.zeros((2,3,3,))
    matrix[0] = np.identity(3)
    matrix[1] = np.identity(3)
    for ii, target in enumerate(['dim', 'primitive_axis']):
        if target in params.keys():
            for i, num in enumerate(params[target].split()):
                i1 = int(i/3)
                i2 = int(i%3)
                matrix[ii][i1,i2] = float(Fraction(num))
    ##
    scell = matrix[0]
    primat = matrix[1]
    return scell, primat

