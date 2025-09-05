#
# vasp.py
#
# This file helps to perform VASP jobs.
#
# Copyright (c) 2022 Masato Ohnishi
#
# This file is distributed under the terms of the MIT license.
# Please see the file 'LICENCE.txt' in the root directory
# or http://opensource.org/licenses/mit-license.php for information.
#
import os
import glob
import numpy as np
import xmltodict
import warnings
from pymatgen.io.vasp.outputs import Vasprun
from phonopy.interface.vasp import read_vasp
from pymatgen.io.vasp.inputs import Incar, Kpoints
import ase.io

from auto_kappa.units import BohrToA, RyToEv
from auto_kappa.structure.two import get_normal_index

import logging
logger = logging.getLogger(__name__)

def get_dfset(directory, offset_xml=None, outfile=None, nset=None, fd2d=False, use_mlips=False):
    """ Get dataset of displacements and forces from many vasprun.xml files.
    
    Args
    -----
    directory : string
        vasprun.xml can be found under \${directory}/\*/ 
        while \${directory}/prist/vasprun.xml is ignored.
    
    offset_xml : string
        vasprun.xml name for offset
    
    outfile : string
        if given, dfset is saved.
    
    nset : int
        Number of data set. If not given, all data will be read.
    
    fd2d : bool
        If True, the displacement method was finite difference for 2D systems.
        The translational displacement is allowed along the out-of-plane direction.
    """
    from ase.geometry import get_distances
    
    if use_mlips:
        line = directory + '/*/forces.xyz'
    else:
        line = directory + '/*/vasprun.xml'
    all_files = glob.glob(line)
    
    ## get keys
    all_keys = []
    for fn in all_files:
        dirname = os.path.dirname(fn)
        key = fn.split('/')[-2]
        if key == 'prist':
            continue
        # For MLIPS, check forces.xyz instead of vasprun.xml
        if use_mlips:
            finished = wasfinished_mlips(dirname)
        else:
            finished = wasfinished(dirname)
        
        if not finished:
            continue
        all_keys.append(key)
    all_keys = sorted(all_keys, key=int)
    
    ## get offset data
    try:
        if use_mlips:
            atoms0 = ase.io.read(offset_xml)
        else:
            atoms0 = ase.io.read(offset_xml, format='vasp-xml')
    except Exception as e:
        warnings.warn(f" Cannot read {offset_xml}: {e}")
        return None

    forces0 = atoms0.get_forces()
    positions0 = atoms0.get_positions()
    
    if fd2d:
        norm_idx = get_normal_index(atoms0)
    else:
        norm_idx = None
    
    ## 
    def _get_diagonal_elements(flarge):
        n = len(flarge)
        fdiag = np.zeros((n,3))
        for i in range(n):
            fdiag[i,:] = flarge[i,i,:]
        return fdiag
    
    ## read all dfset
    all_disps = []
    all_forces = []
    all_energies = []
    all_lines = []
    # count = 0
    for ii, key in enumerate(all_keys):
        if use_mlips:
            fn = "%s/%s/forces.xyz" % (directory, key)
            atoms1 = ase.io.read(fn)
            print("Reading forces.xyz file: ", fn)
        else:
            fn = "%s/%s/vasprun.xml" % (directory, key)
            atoms1 = ase.io.read(fn, format='vasp-xml')
            print("Reading vasprun.xml file: ", fn)        
        
        disp_large, _ = get_distances(
                positions0, p2=atoms1.get_positions(),
                cell=atoms0.cell, pbc=True)
        displacements = _get_diagonal_elements(disp_large) / BohrToA  ## Bohr
        forces = (atoms1.get_forces() - forces0) * BohrToA / RyToEv   ## Bohr/Ry
        
        if fd2d:
            displacements = _adjust_displacements_with_trans_disp(
                displacements, normal_index=norm_idx)
        
        ene = atoms1.get_potential_energy()
        all_energies.append(ene)     ## eV
        
        all_disps.append(displacements)
        all_forces.append(forces)
        
        ##
        rel_path = (("./" + os.path.relpath(fn, os.getcwd())
                      if os.path.isabs(fn) else fn))
        ## get lines
        all_lines.append("# Filename: %s, Snapshot: %d, E_pot (eV): %.7f" % (
            rel_path, ii+1, ene))
        for (d, f) in zip(displacements, forces):
            line = "%16.13f " * 3 % tuple(d)
            line += "  "
            line += "%20.13e " * 3 % tuple(f)
            all_lines.append(line)
        ##
        # count += 1
        # if nset is not None:
        #     if count == nset:
        #         break
    ##
    all_disps = np.asarray(all_disps)
    all_forces = np.asarray(all_forces)
    
    ##
    if outfile is not None:
        ofs = open(outfile, 'w')
        ofs.write("\n".join(all_lines))
        ofs.close()
        
        msg = "\n Output " + outfile
        logger.info(msg)

    return [all_disps, all_forces]

def _adjust_displacements_with_trans_disp(displacements, normal_index=None, tolerance=1e-8):
    """ Adjust displacements with translational displacement.
    The translational displacement is allowed along the normal direction.
    
    Args
    -----
    displacements : ndarray
        Displacement matrix of shape (nset, natoms, 3)
    
    normal_index : int or None
        If not None, the index of the normal direction.
        The translational displacement is allowed along this direction.
    """
    if normal_index is None:
        return displacements
    
    disps_norm = displacements[:, normal_index]
    
    vmin = np.min(disps_norm)
    vmax = np.max(disps_norm)
    if vmax - vmin < 1e-7:
        return displacements
    
    ## Gest the best shift with the tolerance
    from collections import Counter
    shift_candidates = - disps_norm
    rounded = np.round(shift_candidates / tolerance) * tolerance
    count = Counter(rounded)
    best_shift, max_zeros = count.most_common(1)[0]
    
    ## Gest the exact best shift
    imin = np.argmin(abs(disps_norm + best_shift))
    best_shift = - disps_norm[imin]
    
    ## Make new displacements
    disp_new = displacements.copy()
    disp_new[:, normal_index] += best_shift
    return disp_new
    

def read_dfset(filename, natoms=None, nstructures=None):
    """ Read DFSET with Alamode format and return displacements and forces.
    If natoms and nstrctures are set, size of matrices are determined with
    them.
    """
    data = np.loadtxt(filename)
    ntot = len(data)
    if nstructures is not None:
        natoms = int(ntot / nstructures)
    
    ### Obtain natoms automatically
    if natoms is None:
        lines = open(filename, 'r').readlines()
        nstructures = 0
        for line in lines:
            if "#" == line[0]:
                nstructures += 1
        if ntot%nstructures != 0:
            warnings.warn(" Numbers of atoms and structures are "\
                    "inconsistent in %s" % filename)
        natoms = int(ntot / nstructures)
    
    nstructures = int(ntot / natoms)
    
    ###
    data = data.reshape((-1, natoms, 6))
    
    ### extract displacements and forces
    disps = np.zeros((nstructures, natoms, 3))
    forces = np.zeros((nstructures, natoms, 3))
    for ii in range(nstructures):
        disps[ii,:,:] = data[ii,:,:3]
        forces[ii,:,:] = data[ii,:,3:]
    
    return disps, forces

def wasfinished_mlips(directory):
    """Check if MLIPS calculation is finished by checking if forces.xyz exists and has content.
    
    Args:
        directory: Path to the directory to check
        
    Returns:
        bool: True if forces.xyz exists and has content, False otherwise
    """
    fn_target = directory + '/forces.xyz'
    try:
        with open(fn_target, 'r') as f:
            lines = f.readlines()
        # Check if file has at least 3 lines (natoms, comment, first atom)
        if len(lines) >= 3:
            return True
        else:
            return False
    except Exception:
        return False

def wasfinished(directory, filename='vasprun.xml', tar=None):
    """ 
    Args
    ------
    directory : string
        Path for the directory to be checked. Note that the path should be that
        inside the tar.gz file when tar is not None.
    """
    fn_target = directory + '/' + filename
    if tar is None:
        try:
            lines = open(fn_target, 'r').readlines()
        except Exception:
            return False
    else:
        lines_tmp = tar.extractfile(fn_target).readlines()
        lines = [ll.decode('utf-8') for ll in lines_tmp]
     
    n = len(lines)
    for i in range(n):
        num = n - 1 - i
        data = lines[num].split()
        if len(data) != 0:
            if data[0] == '</modeling>':
                return True
            else:
                return False
    return False

def read_poscar(filename):
    """ Read filename
    filename : string, POSCAR filename
    """
    structure = None
    try:
        structure = read_vasp(filename)
    except Exception:
        warnings.warn(" Warning: cannot find %s" % filename)
    return structure

def read_incar(filename):
    incar = None
    try:
        incar = Incar.from_file(filename)
    except Exception:
        warnings.warn(" Warning: cannot find %s" % filename)
    return incar

def read_kpoints(filename):
    kpoints = None
    try:
        kpoints = Kpoints.from_file(filename)
    except Exception:
        warnings.warn(" Warning: cannot find %s" % filename)
    return kpoints

def write_born_info(filename, outfile='BORNINFO'):
    """
    Args
    -------
    filename (str) : vasprun.xml
    """
    lines = []
    
    vasprun = Vasprun(filename, parse_potcar_file=False)
    dielectric_tensor = vasprun.epsilon_static
    born_charges = get_born_charges(filename)
    
    ## dielectric tensor
    for j in range(3):
        lines.append("%18.13f " * 3 % tuple(dielectric_tensor[j]))
    
    ## Born effective charge
    for ia in range(len(born_charges)):
        for j in range(3):
            lines.append("%18.13f " * 3 % tuple(born_charges[ia][j]))
    
    lines.append("")
    f = open(outfile, 'w')
    f.write('\n'.join(lines))

def get_born_charges(filename):
    """
    Args
    ---------
    filename (str) : vasprun.xml
    """
    out = None
    with open(filename, 'r') as f:
        out = xmltodict.parse(f.read())
    array = out['modeling']['calculation']['array']
    
    ## Check
    if array['@name'] != 'born_charges':
        warnings.warn(' Cannot find born_charges in %s' % filename)
        return None
    
    ## number of atoms
    natoms = len(array['set'])
    
    ## Read contents in "born_charges"
    borns = []
    for i in range(natoms):
        
        ### modified on April 17, 2023
        if natoms == 1:
            lines = array['set']['v']
        else:
            lines = array['set'][i]['v']
        
        born = np.zeros((3,3))
        for i1, line in enumerate(lines):
            data = line.split()
            for i2 in range(3):
                born[i1,i2] = float(data[i2])
        borns.append(born)
    return borns

def print_vasp_params(params):
    """ print the contents of a dictionary """
    msg = ""
    for name in params.keys():
        msg += ("\n " +  str(name.upper().ljust(13)) + " = " + 
                str(params[name]))
    logger.info(msg)

def read_outcar(filename):
    """ 
    Args
    -----

    filename : string
        path for OUTCAR
    """
    dump = {}
    lines = open(filename, 'r').readlines()
    for line in lines:
        data = line.split()
        if "Total CPU time used (sec)" in line:
            idx = -1
            key = "cpu_time(min)"
            scale = 1./60.
        elif "Maximum memory used (kb)" in line:
            idx = -1
            key = "max_memory(GB)"
            scale = 1e-3
        else:
            continue

        try:
            dump[key] = float(data[idx]) * scale
        except Exception:
            pass
    
    return dump
