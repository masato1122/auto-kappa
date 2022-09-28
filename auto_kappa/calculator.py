# -*- coding: utf-8 -*-
import os.path
import numpy as np
from pymatgen.io.vasp.inputs import Kpoints

from ase.calculators.vasp import Vasp
from ase.calculators.vasp.create_input import GenerateVaspInput
import re

def run_vasp_with_custodian(calc, atoms, max_errors=10):
    """ Run a VASP job with Custodian
    Args
    ------
    calc : ASE calculator
    atoms : ASE Atoms obj
    """
    from custodian.custodian import Custodian
    from custodian.vasp.handlers import (
            VaspErrorHandler, UnconvergedErrorHandler)
    from custodian.vasp.jobs import VaspJob
    
    ### prepare a directory and write files
    os.makedirs(calc.directory, exist_ok=True)
    calc.write_input(atoms)
    
    ### change directory
    cwd = os.getcwd()
    outdir = calc.directory
    os.chdir(outdir)
    
    ### set handlers and run VASP
    handlers = [VaspErrorHandler()]
    ##, UnconvergedErrorHandler()]
    
    vjob = VaspJob(calc.command.split(), auto_npar=True)
    
    c = Custodian(handlers, [vjob], max_errors=max_errors)
    c.run()
    
    ### return to the initial directory
    os.chdir(cwd)
    return 1

def run_vasp(calc, atoms, method='custodian', max_errors=10):
    """ Run a VASP job 
    Args
    -----
    calc : ASE VASP calculator
    atoms : ASE Atoms obj
    method : string
        "ase" or "custodian"
    max_errors : int
        parameter for "custodian"
    """
    if calc.directory is None:
        warnings.warn(" WARNING: "\
                "output directory is not set in the calculator.")
        exit(0)
    
    if 'ase' in method.lower():
        
        atoms.calc = calc
        atoms.get_potential_energy()
        return 1
    
    elif 'custodian' in method.lower():
        
        value = run_vasp_with_custodian(calc, atoms, max_errors=max_errors)
        return value
    
    else:
        print(" Error: method %s is not supported." % (method))
        exit(0)

def get_vasp_calculator(mode, atoms=None, directory=None, kpts=None,
        encut_scale_factor=1.3,
        auto_lreal_scell_size=65,
        setups='recommended', xc='pbesol',
        ):
    """ Get VASP parameters for the given mode. Parameters are similar to those
    used for phonondb.
    
    Args
    -------
    mode : string
        "relax", "force", "nac", or "md"
    
    atoms : ASE Atoms object
    
    directory : string
        output directory
    
    kpts : list of float, shape=(3)
    
    How to Use
    -----------
    >>> from auto_kappa.calculator import get_vasp_calculator
    >>> 
    >>> atoms = ase.io.read('POSCAR.primitive', format='vasp')
    >>> mode = 'relax'
    >>>
    >>> #atoms = ase.io.read('POSCAR.supercell', format='vasp')
    >>> #mode = 'force'
    >>>
    >>> calc = get_vasp_calculator(mode,
    >>>     directory='./out',
    >>>     kpts=[10,10,10])
    >>> calc.command = "mpirun -n 2 vasp"
    >>> calc.write_input(structure)
    
    """
    from auto_kappa import default_vasp_parameters
    
    calc = Vasp(setups=setups, xc=xc)
    
    ### initialization
    calc.initialize(atoms)
    
    ### set defualt parameters
    params = default_vasp_parameters[mode.lower()].copy()
    
    ### shared parameters
    params.update(default_vasp_parameters['shared'])
    
    ### encut
    enmax = get_enmax(calc.ppp_list)
    params['encut'] = enmax * encut_scale_factor
    
    ### set LREAL
    if len(atoms) >= auto_lreal_scell_size:
        params['lreal'] = 'Auto'
    else:
        params['lreal'] = False
    
    ### kpoints
    if kpts is not None:
        calc.set(kpts=kpts)
        params['gamma'] = True
    
    ### output directory
    if directory is None:
        outdir = 'out_%s' % mode
    else:
        outdir = directory
    calc.directory = outdir
    
    ### Generate
    vasp = GenerateVaspInput.set(calc, **params)
    calc.results.clear()
    return calc

def get_enmax(ppp_list):
    """
    Args
    ----------
    ppp_list : list of string
        ppp_list[*] is a POTCAR file

    Return
    ---------
    Maximum value of ENMAX
    """
    enmaxes = []
    for filename in ppp_list:
        lines = open(filename, 'r').readlines()
        line = [l for l in lines if 'ENMAX' in l]
        if len(line) == 0:
            print(' Error')
            return 300.
        data = re.split(r'\s+|;|=', line[0])
        data2 = [d for d in data if d != '']
        enmaxes.append(float(data2[1]))
    ##
    enmaxes = np.asarray(enmaxes)
    return np.max(enmaxes)


#def get_recommended_kpoints(struct_given):
#    from .structure.format import change_structure_format
#    structure = change_structure_format(
#            struct_given, format='pymatgen-structure')
#    return Kpoints.automatic_density(structure).num_kpts

