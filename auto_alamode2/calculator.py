# -*- coding: utf-8 -*-
import os.path
import numpy as np
from pymatgen.io.vasp.inputs import Kpoints

from . import calculator_parameters

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
    outdir = cwd + '/' + calc.directory
    os.chdir(outdir)
    
    ### set handlers and run VASP
    handlers = [VaspErrorHandler(), UnconvergedErrorHandler()]

    vjob = VaspJob(
            calc.command.split(),
            auto_npar=True,
            )
    
    c = Custodian(handlers, [vjob], max_errors=max_errors)
    c.run()
    
    ### return to the initial directory
    os.chdir(cwd)

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
        warnings.warning(" WARNING: "\
                "output directory is not set in the calculator.")
        exit(0)
    
    #print("")
    if 'ase' in method.lower():
        
        #print(" ### Run VASP with ASE")
        atoms.calc = calc
        atoms.get_potential_energy()
    
    elif 'custodian' in method.lower():
        
        #print(" ### Run VASP with Custodian")
        run_vasp_with_custodian(calc, atoms, max_errors=max_errors)
    
    else:
        print(" Error: method %s is not supported." % (method))
        exit(0)
    #print("")

def get_vasp_parameters(mode, directory=None, atoms=None, kpts=None,
        encut_scale_factor=1.3,
        auto_lreal_scell_size=65,
        setups='recommended', xc='pbesol',
        ):
    """ Get VASP parameters for the given mode. Parameters are similar to those
    used for phonondb.
    
    Args
    -------
    #calc : ASE calculator
    mode : string
        "relax", "force", "nac", or "md"
    directory : string
        output directory
    atoms : ASE Atoms object
    kpts : list of float, shape=(3)
    """
    from . import default_vasp_parameters
    
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
    GenerateVaspInput.set(calc, **params)
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

#def get_calculator(
#        xc=calculator_parameters['xc'], 
#        outdir='./out',
#        print_comment=True
#        ):
#    ##skipifexists=True,
#    """ Calculator with ASE 
#    Example
#    ==============
#    >>> calc = get_calculator()
#    
#    >>> calc.read_kpoints(fn_kp)
#    or
#    >>> Kpoints.automatic_density(structure)
#    
#    >>> calc.read_incar(fn_incar)
#
#    >>> atoms.calc = calc
#    
#    structure : Pymatgen Structure
#    atoms : ASE Atoms
#    """
#    from ase.calculators.vasp import Vasp
#    cmd = "%s -np %d %s" % (
#            calculator_parameters['mpirun'], nprocs, 
#            calculator_parameters['vasp'])
#    
#    ##
#    calc = Vasp(
#            setups=calculator_parameters['setups'], 
#            command=cmd, directory=outdir, xc=xc)
#    ##
#    #xml = outdir + '/vasprun.xml'
#    #if os.path.exists(xml) and skipifexists:
#    #    lines = open(xml, 'r').readlines()
#    #    for il in range(len(lines)):
#    #        num = len(lines) - 1 - il
#    #        data = lines[num].split()
#    #        if len(data) != 0:
#    #            if data[0] == "</modeling>":
#    #                if print_comment:
#    #                    print(" Job in %s was finished. Skip..." % (outdir))
#    #                return None
#    #            else:
#    #                return calc
#    return calc
#
#def get_recommended_kpoints(struct_given):
#    from .structure.format import change_structure_format
#    structure = change_structure_format(
#            struct_given, format='pymatgen-structure')
#    return Kpoints.automatic_density(structure).num_kpts

