from monty.json import MSONable
import warnings
import numpy as np
import pymatgen
from ..units import AToBohr

#__all__ = ['alm_variables', 'str_vars', 'int_vars', 'double_vars',
#        'array_str_vars', 'array_int_vars', 'array_double_vars', 
#        'hyphen_int_vars']

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
            'optimize':{
                'always':{
                    'dfset'
                    },
                'non-ols':{
                    'always':{
                        'fc2xml', 'conv_tol'
                        },
                    'cv!=5':{
                        'cv_minalpha', 'cv_maxalpha', 'cv_nalpha'
                        }
                    }
                }
            }
    #'l1_ratio', 'l1_alpha', 
    #'cv_maxalpha', 'cv_minalpha',
    
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
            'always':{
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
                },
            'non-ols':{
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
                },
            'enet':{
                'standardize'     # 1,   integer
                }
            }
    
    def __init__(
            self, structure=None, filename=None, **kwargs
            ):
        """
        Args
        ------
        structure (pymatgen.core.structure.Structure object): structure
        filename (str): structure file name
            format must be supported by pymatgen.loadfn module
        kwargs (dict): other alm parameters. Some parameters have to be given to
            run alm.
            - mode
        """
        ## read a structure
        self.structure = None
        if structure is not None:
            self.structure = structure
        else:
            if filename is not None:
                self.structure = pymatgen.loadfn(filename)
        
        ## parameters determined by the structure
        if self.structure is not None:
            params = get_alm_variables_from_structure(self.structure)
            self.update(params)
        
        ## update
        self.update(kwargs)
        
    def set_prefix(self):
        self['prefix'] = self.structure.get_primitive_structure().formula
    
    #@classmethod
    #def from_dict(cls):
    
    #@classmethod
    #def from_file(cls, filepath: str):
    
    @classmethod
    def from_structure(cls, param, guess_primitive=False):
        """ Creat ALM variables from a structure
        Args
        --------
        param (str) or (pymatgen Structure):
            If param is string (structure file), the file will be read with
            pymatgen.loadfn.
        guess_primitive (bool):
            If it's True, the primitive structure will be guessed with spglib.
        """
        return 0
        #exit()
    
    def to_file(self, filename=None):
        if filename is None:
            self.outfile = self['mode']+'.in'
        ##
        for k1 in self.required_params.keys():
            for param in self.required_params[k1]:
                if param not in self.as_dict():
                    warnings.warn(
                        "Required parameter '{}' not specified!".format(param)
                    )
        ##
        general_dict = _get_subdict(self, self.general_keys)
        print(general_dict)  

    #def get_structure(self):
    
    def as_dict(self):
        return dict(self)

## This part should be modified (rewrite).
def _get_subdict(master_dict, subkeys):
    """Helper method to get a set of keys from a larger dictionary"""
    return {
            k: master_dict[k]
            for k in subkeys
            if k in master_dict and master_dict[k] is not None
    }

def get_alm_variables_from_structure(structure):
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
    
    ##
    params['nkd'] = []
    params['kd'] = []
    for sym in sym_list:
        params['nkd'].append(sym_all.count(sym))
        params['kd'].append(sym)
    
    ## cell (with the unit of Bohr) and positions
    params['cell'] = structure.lattice.matrix * AToBohr
    params['position'] = structure.frac_coords
    
    return params

class AnphonInput():
    """ Class for writing and reading anphon input file.
    See the following URL for more detailed description.
    https://alamode.readthedocs.io/en/latest/anphondir/inputanphon.html
    
    Reference
    -------------
    pymatgen.io.shengbte.Control
    """
    general_keys = [
            'bconnect',  
            'borninfo',  
            'bornsym',  
            'classical',  
            'emin',  
            'epsilon',  
            'fc2xml',  
            'fcsxml',  
            'ismear',  
            'kd',  
            'mass',  
            'mode',  
            'na_sigma',  
            'nkd',  
            'nonanalytic',  
            'prefix',  
            'printsym',  
            'restart',  
            'tmin',  
            'tolerance',  
            'trisym'
            ]
    scph_keys = [
            'ialgo',  
            'kmesh_interpolate',  
            'kmesh_scph',  
            'lower_temp',  
            'maxiter',  
            'mixalpha',  
            'restart_scph',  
            'self_offdiag',  
            'tol_scph',  
            'warmstart',  
            ]
    analysis_keys = [
            'anime',  
            'anime_frames',  
            'anime_cellsize',  
            'gruneisen',  
            'isofact',  
            'isotope',  
            'kappa_coherent',  
            'kappa_spec',  
            'pdos',  
            'printevec',  
            'printmsd',  
            'printpr',  
            'printvel',  
            'printxsf',
            'sps',
            'tdos'
            ]
    
    #@classmethod
    #def from_dict(cls):
    
    #@classmethod
    #def from_file(cls, filepath: str):
    
    #@classmethod
    #def from_structure(cls, structure):

    #@classmethod
    #def to_file(cls):
    
    #def get_structure(self):
    
    #def as_dict(self)
     
