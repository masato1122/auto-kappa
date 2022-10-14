#
# alm.py
#
# This file helps to generate input scripts for Alamode.
#
# Copyright (c) 2022 Masato Ohnishi
#
# This file is distributed under the terms of the MIT license.
# Please see the file 'LICENCE.txt' in the root directory
# or http://opensource.org/licenses/mit-license.php for information.
#
from typing import Any, Dict, List, Optional, Union
from monty.json import MSONable
import warnings
import re
import numpy as np
import pymatgen
from pymatgen.core.structure import IStructure
from ase import Atoms

from auto_kappa.units import AToBohr
from auto_kappa.structure import (change_structure_format, get_primitive_structure)
from auto_kappa.structure.crystal import get_formula

try:
    import f90nml
except ImportError:
    f90nml = None

class AlmInput(MSONable, dict):
    """ Class for writing and reading alm input file.
    See the following URL for more detailed description.
    https://alamode.readthedocs.io/en/latest/almdir/inputalm.html
    
    Reference
    -------------
    pymatgen.io.shengbte.Control
    """
    required_params = {
            'always':{
                'prefix', 'mode', 'nat', 'nkd', 'kd',     # general
                'norder', 'cell', 'cutoff', 'position'    # others
                },
            'norder>=4':{
                'nbody'
                },
            'optimize':{
                'always':{
                    'dfset'
                    },
                'lmodel!=ols':{
                    'always':{
                        'fc2xml', 'conv_tol'
                        },
                    'cv!=5':{
                        'cv_minalpha', 'cv_maxalpha', 'cv_nalpha'
                        }
                    }
                }
            }
    
    ##'l1_ratio', 'l1_alpha', 
    ##'cv_maxalpha', 'cv_minalpha',
    
    general_keys = [
            'prefix',       # None, string
            'mode',         # None, string (suggest or optimize)
            'nat',          # None, integer
            'nkd',          # None, integer
            'kd',           # None, array of strings
            'tolerance',    # 1e-3, double
            'printsym',     # 0,    integer
            'fcsym_basis',  # Lattice, string
            'magmom',       # 0...0, array of double
            'noncollinear', # 0,     integer
            'periodic',     # 1 1 1, array of integers
            'nmaxsave',     # min(5, norder), integer
            'hessian',      # 0, integer
            ]
    interaction_keys = {
            'norder',       # None, integer
            'nbody',        # [2,3,...,norder+1], array of integers
            }
    optimize_keys = {
            ## always
            'lmodel',         # least-squares, string
            'dfset',          # None,  string
            'ndata',          # None,  integer
            'nstart',         # 1,     integer
            'nend',           # ndata, integer
            'skip',           # None,  int-int
            'iconst',         # 11,    integer
            'rotaxis',        # None,  string
            'fc2xml',         # None,  string
            'fc3xml',         # None,  string
            'sparse',         # 0,     integer
            'sparsesolver',   # SimplicialLDLT, string
            ## for non-ols
            'maxiter',        # 10000,       integer
            'conv_tol',       # 1e-8,        double 
            'l1_ratio',       # 1.0 (LASSO), double
            'l1_alpha',       # 0.0,      double
            'cv',             # 0,        integer
            'dfset_cv',       # dfset,    string
            'ndata_cv',       # None,     integer
            'nstart_cv',      # 1,        integer
            'nend_cv',        # ndata_cv, integer
            'cv_minalpha',    # cv_maxalpha*1e-6,  double
            'cv_maxalpha',    # set automatically, double
            'cv_nalpha',      # 50,  integer
            'enet_dnorm',     # 1.0, double
            'solution_path',  # 0,   integer
            'debias_ols',     # 0,   integer
            ## for lmode == 'enet'
            'standardize'     # 1,   integer
            }
    
    def __init__(self, **kwargs):
        """
        Args
        ------
        kwargs (dict): other alm parameters. Some parameters have to be given to
            run alm.
            - mode
        """
        super().__init__()
        self.update(kwargs)
    
    @property
    def primitive(self):
        return self.get_primitive()

    def get_primitive(self):
        return get_primitive_structure(self.structure)

    @classmethod
    def from_dict(cls, alm_dict: Dict):
        """ Write a ALM input file
        """
        return cls(**alm_dict)
    
    #@classmethod
    #def from_file(cls, filepath: str):
    
    @classmethod
    def from_structure_file(cls, filename, norder=None, **kwargs):
        """ Make ALM input parameters from a structure file
        Args
        -------
        filename (str): filename
            Format must be supported by pymatgen's IStructure module: POSCAR,
            .cif. etc.
        """
        ## set norder
        kwargs['norder'] = norder
        
        ## set a structure
        cls.structure = IStructure.from_file(filename)
        alm_dict = dict(
                cls.from_structure(cls.structure, norder=norder)
                )
        
        alm_dict.update(kwargs)
        return AlmInput(**alm_dict)
    
    @classmethod
    def from_structure(cls, structure, norder=None, **kwargs):
        """ Creat ALM variables from a structure object or a structure file
        
        Args
        --------
        structure (str, pymatgen's (I)Structure, ase's Atoms, or structure file 
        such as POSCAR, .cif, etc.):
            If param is string (structure file), the file will be read with
            pymatgen.loadfn.
        guess_primitive (bool):
            If it's True, the primitive structure will be guessed with spglib.
        """
        alm_dict = {}
        
        ## set norder
        kwargs['norder'] = norder
        
        ## set a structure
        cls.structure = change_structure_format(
                structure, format='pmg-IStructure')
        
        ## parameters determined by the structure
        params_str = get_alamode_variables_from_structure(
                cls.structure, norder=norder)
        alm_dict.update(params_str)

        ## set given parameters
        alm_dict.update(kwargs)
        
        ## set prefix
        if 'prefix' not in alm_dict.keys():
            alm_dict['prefix'] = get_formula(cls.structure)
        else:
            if alm_dict['prefix'] is None:
                alm_dict['prefix'] = get_formula(cls.structure)
        
        return AlmInput(**alm_dict)
    
    def check_parameters(self):
        
        ## check parameters which a re always required
        for param in self.required_params['always']:
            _check_dict_contents(self.as_dict(), param)
        
        ## check 'norder'
        if 'norder' in self.as_dict().keys():
            _check_dict_contents(self.as_dict(), 'norder')
        
        ## Set nbody when higher-order IFCs are calculated
        try:
            if self['norder'] >= 4:
                if _check_dict_contents(
                        self.as_dict(), 'nbody', message=False) == False:
                    self['nbody'] = list(np.arange(2, self['norder']+1))
        except:
            pass
        
        ## When LMODE!=OLS
        try:
            if _lmodel_is_ols(str(self['lmodel'].lower())) == False:
                for param in self.required_params['lmodel!=ols']:
                    _check_dict_contents(self.as_dict(), param)
        except:
            pass
        
        ## check required parameters for 'mode == optimize'
        mode = self.as_dict()['mode'].lower()
        if mode == 'optimize':
            for param in self.required_params[mode]['always']:
                _check_dict_contents(self.as_dict(), param)

    def to_file(self, filename: str = None):
        """ Make ALM input file
        Args
        -------
        filename (str): ALM input file name
        """
        self.check_parameters()
        
        ## output file name
        if filename is None:
            filename = self['mode']+'.in'
        
        ## general parameters
        general_dict = _get_subdict(self, self.general_keys)
        general_nml = f90nml.Namelist({"general": general_dict})
        alm_str = str(general_nml) + "\n"
        
        ## interaction parameters
        interaction_dict = _get_subdict(self, self.interaction_keys)
        interaction_nml = f90nml.Namelist({"interaction": interaction_dict})
        alm_str += str(interaction_nml) + "\n"
        
        ### cell parameters
        lines = _write_cell(self['cell'])
        for line in lines:
            alm_str += line + "\n"

        ### cutoff parameters
        lines = _write_cutoff(self['cutoff'])
        for line in lines:
            alm_str += line + "\n"
         
        ## optimize parameters
        optimize_dict = _get_subdict(self, self.optimize_keys)
        optimize_nml = f90nml.Namelist({"optimize": optimize_dict})
        alm_str += str(optimize_nml) + "\n"
         
        ### position parameters
        lines = _write_positions(self['position'])
        for line in lines:
            alm_str += line + "\n"
        
        ## characters to be replaced
        for char in [',', '\'']:
            alm_str = alm_str.replace(char, "")
        
        ## output an ALM input file
        with open(filename, "w") as file:
            file.write(alm_str)

    def as_dict(self):
        return dict(self)
    
    #def _get_formula(self):
    #    words = self.structure.get_primitive_structure().formula.split()
    #    if len(words) == 1:
    #        return re.sub(r"[123456789]", "", words[0])
    #    else:
    #        prefix = ""
    #        for i in range(len(words)):
    #            prefix += words[i].replace("1", "")
    #        return prefix
    
    #def suggest_cutoff(self):
    #    cutoff = suggest_alm_cutoff(self.structure)
    

## This part should be modified (rewrite).
def _get_subdict(master_dict, subkeys):
    """Get a set of keys """
    """Helper method to get a set of keys from a larger dictionary"""
    return {
            k: master_dict[k] 
            for k in subkeys 
            if k in master_dict and master_dict[k] is not None
            }

def _lmodel_is_ols(lmodel):
    if str(lmodel.lower()) in ['least-squares', 'ls', 'ols', '1']:
        return True
    else:
        return False

##
def _write_cell(cell):
    """
    Args
    ------
    cell (ndarray), shape=(3,3), unit=[Bohr]
    """
    ## set length unit
    lengths = np.zeros(3)
    for j in range(3):
        lengths[j] = np.linalg.norm(cell[j])
    lunit = np.min(lengths)
    
    lines = []
    lines.append("&cell")
    lines.append(" %18.13f"%(lunit))
    for i in range(3):
        lines.append(" %18.13f" * 3 % tuple(cell[i]/lunit))
    lines.append("/")
    return lines

def _write_cutoff(cutoff):
    """
    Args
    -------
    cutoff (dict of array of float):
        keys are two strings connected with hyphene, e.g. 'Si-Si'
    """
    lines = []
    lines.append("&cutoff")
    for k1 in cutoff.keys():
        lines.append("    %7s  "%(k1))
        for val in cutoff[k1]:
            if val is None:
                lines[-1] += "None "
            else:
                lines[-1] += "%18.13f "%(val)
    lines.append("/")
    return lines

def _write_positions(positions):
    """
    Args
    -------
    position (ndarray), shape=(natoms, 3):
        scaled positions
    """
    lines = []
    lines.append("&position")
    for ia, pos in enumerate(positions):
        lines.append(" %3d "%(pos[0]))
        lines[-1] += " %18.13f" * 3 % tuple(pos[1:4])
    lines.append("/")
    return lines

def get_alamode_variables_from_structure(structure, norder=None):
    """ Get alm variables determined by the structure
    """
    params = {}
    
    ## number of atoms
    params['nat'] = len(structure)
    
    ## get list of chemical symbols
    sym_all = []
    for each in structure.species:
        sym_all.append(each.name)
    sym_list = sorted(set(sym_all), key=sym_all.index)
    
    ## get number of atoms for each species
    params['nkd'] = []
    params['kd'] = []
    params['nkd'] = len(sym_list)
    for sym in sym_list:
        params['kd'].append(sym)
    
    ## cell (with the unit of Bohr) 
    params['cell'] = structure.lattice.matrix * AToBohr
    
    ## set positions
    params['position'] = []
    for ia, atom in enumerate(structure):
        sym = atom.specie.name
        idx = sym_list.index(sym) + 1
        p = atom.frac_coords
        params['position'].append([idx, p[0], p[1], p[2]])
    
    ## cutoff (unit: Bohr)
    params['cutoff'] = {}
    for s1 in sym_list:
        for s2 in sym_list:
            comb = s1+'-'+s2
            params['cutoff'][comb] = []
            if norder is not None:
                for ii in range(norder):
                    params['cutoff'][comb].append(None)
    
    return params

#def suggest_alm_cutoff(structure):
#    sym_all = structure.species
#    exit()

class AnphonInput(MSONable, dict):
    """ Class for writing and reading anphon input file.
    See the following URL for more detailed description.
    https://alamode.readthedocs.io/en/latest/anphondir/inputanphon.html
    
    Reference
    -------------
    pymatgen.io.shengbte.Control
    """
    required_params = {
            'always':{
                'prefix', 'mode', 'nkd', 'kd', 'fcsxml',    # general
                'cell', 'kpmode'
                },
            }
    general_keys = [
            'prefix',       # None, string
            'mode',         # None, string
            'nkd',          # None, integer
            'kd',           # None, Array of strings
            'mass',         # weight of elements, Array of double
            'fcsxml',       # None, string
            'fc2xml',       # None, string
            'fc3xml',       # None, string
            'tolerance',    # 1e-6, double
            'printsym',     # 0,    integer
            'nonanalytic',  # 0,    integer
            'na_sigma',     # 0.0,  double
            'borninfo',     # None, string
            'bornsym',      # 0,    integer
            'tmin',         # 0,    double
            'tmax',         # 1000, double
            'dt',           # 10,   double
            'emin',         # 0,    double
            'emax',         # 1000, double
            'delta_e',      # 10,   double
            'ismear',       # -1,   integer
            'epsilon',      # 10.0, double
            'bconnect',     # 0,    integer
            'classical',    # 0,    integer
            'trisym',       # 1,    integer
            'restart',      # 1 or 0, integer
            ]
    scph_keys = [
            'kmesh_interpolate', # None, Array of integers 
            'kmesh_scph',        # None, Array of integers
            'self_offdiag',      # 0,    Integer
            'tol_scph',          # 1e-10, double
            'mixalpha',          # 0.1,   double
            'maxiter',           # 1000,  integer
            'lower_temp',        # 1,     integer
            'warmstart',         # 1,     integer
            'ialgo',             # 0,     integer
            'restart_scph',      # 1 or 0, integer
            ]
    analysis_keys = [
            'gruneisen', 
            'printevec', 
            'printxsf',
            'printvel',
            'printmsd',
            'pdos',
            'tdos',
            'sps', 
            'printpr', 
            'kappa_coherent', 
            'kappa_spec',  
            'isotope',  
            'isofact',        # Automatically calculated from the KD-tag
            'anime',          # None, array of doubles 
            'anime_frames',   # 20, integer
            'anime_cellsize', # None, Array of integers
            'anime_format',   # xyz,  string
            ]

    all_keys = (
            general_keys + scph_keys + 
            analysis_keys + ['cell', 'kpoint', 'kpmode']
            )
    
    ## added keys for k-point 
    kpoint_keys = ['kpmode', 'kpoints', 'kpath', 'kpts', 'deltak']
    
    def __init__(self, **kwargs):
        """
        Args
        ------
        kwargs (dict): other anphon parameters. Some parameters have to be 
            given to run alm.
            - mode
        """
        super().__init__()
        self.update(kwargs)

        self._primitive = None
    
    def as_dict(self):
        return dict(self)
    
    @property
    def primitive(self):
        if self._primitive is None:
            self._primitive = self.get_primitive()
        return self._primitive
    
    def set_primitive(self, structure):
        self._primitive = structure

    def get_primitive(self):
        return get_primitive_structure(self.structure)

    @classmethod
    def from_dict(cls, alm_dict: Dict):
        """ Write a ALM input file
        """
        return cls(**alm_dict)
    
    @classmethod
    def from_structure_file(cls, filename, **kwargs):
        """ Make ANPHON input parameters from a structure file
        Args
        -------
        filename (str): filename
            Object should be supported by pymatgen's IStructure module: POSCAR,
            .cif. etc.
        """
        ## get a structure
        cls.structure = IStructure.from_file(filename)
        
        ## get the primitive structure
        primitive = cls.get_primitive(cls)
        
        ## get structure parameters with the primitive structure
        anp_dict = dict(cls.from_structure(primitive))
        
        ## update parameters
        anp_dict.update(kwargs)
        
        return AnphonInput(**anp_dict)
    
    @classmethod
    def from_structure(cls, structure, **kwargs):
        """ 
        Args
        --------
        structure (structure object)
        """
        cls.structure = structure
        
        ## Get parameters determined by structure
        alm_dict = dict(AlmInput.from_structure(structure))
        
        ## Extract only available keys
        anp_dict = {}
        for k1 in alm_dict.keys():
            if k1 in cls.all_keys:
                anp_dict[k1] = alm_dict[k1]
        
        ## Parameters for k-point
        if 'kpmode' not in kwargs.keys():
            kpmode = None
        
        ## update
        anp_dict.update(kwargs)
        
        return AnphonInput(**anp_dict)
     
    #@classmethod
    #def from_file(cls, filepath: str):
    
    def set_kpoint(self, **kwargs):
        """ Set suggeted k-point parameters. This module may not well written.
        
        ## kpoint_keys = ['kpmode', 'kpoints', 'kpath', 'kpts', 'deltak']
        If kpmode == 0, separated kpoints should be given
        If kpmode == 1, deltak (default: 0.01) and kpath will be set.
        If kpmode == 2, deltak (default: 0.2) and kpts will be set.
        """
        ## set kpmode
        if _check_dict_contents(
                self.as_dict(), 'kpmode', message=False) == False:
            if _check_dict_contents(kwargs, 'kpmode'):
                self['kpmode'] = kwargs['kpmode']
            else:
                ValueError(' kpmode must be given.')
        
        #exit()
        ## Check if kpmode key is set or not
        set_kpmode = False
        if 'kpmode' not in self.as_dict().keys():
            set_kpmode = True
        else:
            if self.as_dict()['kpmode'] is None:
                set_kpmode = True
        
        ## set kpmode
        if set_kpmode:
            _check_dict_contents(self.as_dict(), 'mode')
            if self['mode'].lower() == 'phonons':
                self['kpmode'] = 1
            elif self['mode'].lower() == 'scph':
                self['kpmode'] = 1
            elif self['mode'].lower() == 'rta':
                self['kpmode'] = 2
        
        ## set deltak
        if _check_dict_contents(kwargs, 'deltak', message=False):
            self['deltak'] = kwargs['deltak']
        else:
            if self['kpmode'] == 1:
                self['deltak'] = 0.01
            elif self['kpmode'] == 2:
                self['deltak'] = 0.2
        
        ## set k-points
        if self['kpmode'] == 0:
            if _check_dict_contents(kwargs, 'kpoints'):
                self['kpoints'] = kwargs['kpoints']
            else:
                self['kpoints'] = [[0,0,0]]
        elif self['kpmode'] == 1:
            if 'kpath' not in self.keys():
                print(" kpath is set automatically.")
                self['kpath'] = get_kpoint_path(
                        self.primitive, deltak=self['deltak'])
        elif self['kpmode'] == 2:
            lengths = self.primitive.lattice.reciprocal_lattice.lengths
            kpts = []
            for i, leng in enumerate(lengths):
                n = max(2, int(np.ceil(leng/self['deltak'])))
                kpts.append(n)
            self['kpts'] = kpts
        
        self.update(kwargs)

    #def guess_cutoff(self):
    #    pass

    def check_parameters(self):
        
        ## check required parameters
        for param in self.required_params['always']:
            if param not in self.as_dict():
                warnings.warn(
                    "Required parameter '{}' is not specified!".format(param)
                    )
        ##
        if self['kpmode'] == 0:
            _check_dict_contents(self.as_dict(), 'kpoints')
        
        ## This part helps to set BONRINFO parameter, when nonanalytic term will
        ## be considered.
        try:
            if self['nonanalytic'] != 0:
                if _check_dict_contents(self.as_dict(), 'borninfo') == False:
                    self['borninfo'] = 'BORNINFO'
        except:
            pass

    def to_file(self, filename: str = None):
        """ Make ANPHON input file
        Args
        -------
        filename (str): ANPHON input file name
        """
        self.check_parameters()
        
        ## output file name
        if filename is None:
            if self['mode'].lower() == 'phonons':
                if self['kpmode'] == 0: 
                    filename = 'eigen.in'
                elif self['kpmode'] == 1: 
                    filename = 'band.in'
                elif self['kpmode'] == 2: 
                    filename = 'dos.in'
            else:
                filename = self['mode']+'.in'
        
        ## general parameters
        general_dict = _get_subdict(self, self.general_keys)
        general_nml = f90nml.Namelist({"general": general_dict})
        anp_str = str(general_nml) + "\n"
        
        ## cell parameters
        lines = _write_cell(self['cell'])
        for line in lines:
            anp_str += line + "\n"
        
        ## k-point
        lines = _write_kpoint(self.as_dict())
        for line in lines:
            anp_str += line + "\n"
        
        ## analysis parameters
        analysis_dict = _get_subdict(self, self.analysis_keys)
        analysis_nml = f90nml.Namelist({"analysis": analysis_dict})
        anp_str += str(analysis_nml) + "\n"
        
        ## characters to be replaced
        for char in [',', '\'']:
            anp_str = anp_str.replace(char, "")
        
        ## output an ALM input file
        with open(filename, "w") as file:
            file.write(anp_str)
    
def _write_kpoint(params, kpoint=None):
    lines = []
    lines.append('&kpoint')
    lines.append('    %d'%(params['kpmode']))
    if params['kpmode'] == 0:
        for k in params['kpoints']:
            lines.append("    %18.13f " * 3 % tuple(k))
    elif params['kpmode'] == 1:
        for k in params['kpath']:
            lines.append("")
            for j in range(2):
                key = list(k[j].keys())[0]
                lines[-1] += "    %5s" % key
                lines[-1] += " %10.7f " * 3 % tuple(k[j][key])
            lines[-1] += " %d " % k[2]
    elif params['kpmode'] == 2:
        lines.append("    ")
        lines[-1] += "%d " * 3 % tuple(params['kpts'])

    lines.append('/')
    return lines

def _check_dict_contents(dictionary, key, message=True):
    flag = True
    if key not in dictionary:
        flag = False
    else:
        if dictionary[key] is None:
            flag = False
    ##
    if flag == False and message:
        warnings.warn(
                "Required parameter '{}' is not specified!".format(key)
                )
    return flag

def get_kpoint_path(primitive, deltak=0.01):
    """
    Args
    --------
    primitive (pymatgen's (I)Structure object):
        primitive structure
    
    Return
    --------
    kpoints (array of array of dict)
        key of dictionary is a name of symmetric point and the contained value
        is the corresponding k-point.
    """
    import seekpath
    structure = [primitive.lattice.matrix,
            primitive.frac_coords,
            primitive.atomic_numbers]
    kpath = seekpath.get_path(
            structure, with_time_reversal=True, recipe='hpkot')
    
    ## extract required data
    kpoints = []
    for each in kpath['path']:
        kpoints.append([])
        k1 = np.asarray(kpath['point_coords'][each[1]])
        k0 = np.asarray(kpath['point_coords'][each[0]])
        kvec = k1 - k0 
        kleng = np.linalg.norm(kvec)
        nk = int(np.ceil(kleng / deltak))
        kpoints[-1] = [{each[0]: k0}, {each[1]: k1}, nk]
    return kpoints


#def write_alamode_input(mode: None, structure=None, prefix="prefix", 
#        dfset=None, **args):
#    """ Write an input script of Alamode
#    
#    Args
#    -----
#    mode : string, default=None
#        "suggest", "optimize", "phonons"
#
#    """
#    modes_alm = ['suggest', 'optimize']
#    modes_anphon = ['phonons', 'rta', 'scph']
#    
#    if mode in modes_alm:
#        inp = AlmInput.from_structure(
#                structure,
#                mode=mode,
#                )
#    
#
