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
import glob
import numpy as np
import xmltodict
import warnings
from pymatgen.io.vasp.outputs import Vasprun
from phonopy.interface.vasp import read_vasp
from pymatgen.io.vasp.inputs import Incar, Kpoints
import ase.io
from auto_kappa.units import BohrToA, RyToEv

import logging
logger = logging.getLogger(__name__)

def get_dfset(directory, offset_xml=None, outfile=None, nset=None):
    """ Get dataset of displacements and forces from many vasprun.xml files.
    
    Args
    -----
    directory : string
        vasprun.xml can be found under \${directory}/\*/ 
        while \${directory}/prist/vasprun.xml is neglected.
    
    offset_xml : string
        vasprun.xml name for offset
    
    outfile : string
        if given, dfset is saved.
    
    nset : int
        Number of data set. If not given, all data will be read.
    
    """
    from ase.geometry import get_distances
    
    line = directory + '/*/vasprun.xml'
    all_dirs = glob.glob(line)
    
    ## get keys
    all_keys = []
    for dd in all_dirs:
        key = dd.split('/')[-2]
        if key == 'prist':
            continue
        all_keys.append(key)
    all_keys = sorted(all_keys, key=int)
    
    ## get offset data
    try:
        atoms0 = ase.io.read(offset_xml, format='vasp-xml')
    except Exception:
        warnings.warn(" Cannot read %s" % (offset_xml))
        return None

    forces0 = atoms0.get_forces()
    positions0 = atoms0.get_positions()
    
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
    count = 0
    for ii, key in enumerate(all_keys):
        fn = "%s/%s/vasprun.xml" % (directory, key)
        atoms1 = ase.io.read(fn, format='vasp-xml')
        
        disp_large, _ = get_distances(
                positions0, p2=atoms1.get_positions(),
                cell=atoms0.cell, pbc=True)
        displacements = _get_diagonal_elements(disp_large) / BohrToA  ## Bohr
        forces = (atoms1.get_forces() - forces0) \
                * BohrToA / RyToEv                         ## Bohr/Ry
        
        ene = atoms1.get_potential_energy()
        all_energies.append(ene)     ## eV
        
        all_disps.append(displacements)
        all_forces.append(forces)

        ## get lines
        all_lines.append("# Filename: %s, Snapshot: %d, E_pot (eV): %.7f" % (
            fn, ii+1, ene))
        for (d, f) in zip(displacements, forces):
            line = "%16.13f " * 3 % tuple(d)
            line += "  "
            line += "%20.13e " * 3 % tuple(f)
            all_lines.append(line)
        ##
        count += 1
        if nset is not None:
            if count == nset:
                break
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

def wasfinished(directory, filename='vasprun.xml'):
    try:
        fn = directory + '/' + filename
        lines = open(fn, 'r').readlines()
    except Exception:
        return False

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

#
#class Vasprun:
#    """ This class reads vasprun.xml.
#
#    modeling
#        - generator
#        - incar
#        - structure
#        - varray
#        - kpoints
#            - generation
#            - varray
#            - i
#            - 0
#        - parameters
#        - atominfo
#        - calculation
#            - scstep
#            - structure
#            - varray
#            - energy
#            - time
#            - eigenvalues
#            - separator
#            - dos
#    """
#    def __init__(self, filename):
#        
#        self.filename = filename
#        self._number_of_atoms = None
#        self._born_charges = None
#        self._dielectric_tensor = None   # dielectric tensor
#        self._nelect = None
#        self._totalsc = None
#        self._nkpoints = None
#        
#        ## dict
#        self._dict_kpoints = None
#        self._dict_electronic = None
#         
#        with open(filename, 'r') as f:
#            self._doc = xmltodict.parse(f.read())
#    
#    @property
#    def nkpoints(self):
#        return self.get_nkpoints()
#    
#    def get_nkpoints(self):
#        if self._nkpoints is None:
#            self.set_kpoints()
#            self._nkpoints = len(self._dict_kpoints['kpointlist'])
#        return self._nkpoints
#
#    def set_kpoints(self):
#        if self._dict_kpoints is None:
#            kpoints = self._doc['modeling']['kpoints']
#            params = {}
#            
#            ## kpoints/generation
#             
#            ## kpoints/varray
#            varray = kpoints['varray']
#            for each in varray:
#                name = each['@name']
#                values = each['v']
#                n1 = len(values)
#                
#                if name == 'kpointlist':
#                    params[name] = np.zeros((n1,3))
#                    for i1 in range(n1):
#                        data = values[i1].split()
#                        for j, kk in enumerate(data):
#                            params[name][i1,j] = float(kk)
#                
#                elif name == 'weights':
#                    params[name] = np.zeros(n1)
#                    for i1, dd in enumerate(values):
#                        params[name][i1] = float(dd)
#                
#                elif name == 'tetrahedronlist':
#                    params[name] = []
#                    for i1 in range(n1):
#                        params[name].append([])
#                        data = values[i1]['#text'].split()
#                        for j, ii in enumerate(data):
#                            params[name][-1].append(int(ii))
#                
#            ## kpoints/i
#            out = kpoints['i']
#            params[out['@name']] = float(out['#text'])
#        
#            self._dict_kpoints = params
#        
#        return self._dict_kpoints
#
#
#    @property
#    def totalsc(self):
#        return self.get_totalsc()
#
#    def get_totalsc(self):
#        if self._totalsc is None:
#            out = self._doc['modeling']['calculation']['time']
#            if out['@name'] == 'totalsc':
#                times = []
#                for a in out['#text'].split():
#                    times.append(float(a))
#                self._totalsc = times
#        return self._totalsc 
#
#    @property
#    def number_of_atoms(self):
#        return self.get_number_of_atoms()
#    
#    def get_number_of_atoms(self):
#        if self._number_of_atoms is None:
#            natoms = int(self._doc['modeling']['atominfo']['atoms'])
#            self._number_of_atoms = natoms
#        return self._number_of_atoms
#    
#    @property
#    def nelect(self):
#        return self.get_nelect()
#
#    def get_nelect(self):
#        if self._dict_electronic is None:
#            self._set_electronic()
#        self._nelect = self._dict_electronic['NELECT']
#        return self._nelect
#    
#    def _set_electronic(self):
#        if self._nelect is None:
#            out = self._doc['modeling']['parameters']['separator']
#            for each1 in out:
#                if each1['@name'] == 'electronic':
#                    params = {}
#                    for each2 in each1['i']:
#                        key = each2['@name']
#                        if '@type' in list(each2.keys()):
#                            if each2['@type'] == 'string':
#                                params[key] = each2['#text']
#                            elif each2['@type'] == 'int':
#                                params[key] = int(each2['#text'])
#                        else:
#                            params[key] = float(each2['#text'])
#                    break
#        self._dict_electronic = params
#    
#    @property
#    def dielectric_tensor(self):
#        return self.get_dielectric_tensor()
#
#    def get_dielectric_tensor(self):
#        
#        if self._dielectric_tensor is not None:
#            return self._dielectric_tensor
#            
#        array = self._doc['modeling']['calculation']['varray']
#        
#        ## read 'epsilon' in vasprun.xml
#        for each in array:
#            if each['@name'] == 'epsilon':
#                eps = []
#                for i, line in enumerate(each['v']):
#                    data = line.split()
#                    eps.append(np.asarray([float(data[j]) for j in range(3)]))
#        self._dielectric_tensor = np.asarray(eps)
#        return self._dielectric_tensor
#        
#    @property
#    def born_charges(self):
#        return self.get_born_charges()
#
#    def get_born_charges(self):
#        
#        if self._born_charges is None:
#            array = self._doc['modeling']['calculation']['array']
#            
#            ## Check
#            if array['@name'] != 'born_charges':
#                warnings.warn(' Cannot find born_charges in %s' % self.filename)
#                return None
#            
#            ## Read contents in "born_charges"
#            borns = []
#            for i in range(self.number_of_atoms):
#                lines = array['set'][i]['v']
#                born = np.zeros((3,3))
#                for line in lines:
#                    data = line.split()
#                    for j in range(3):
#                        born[j] = float(data[j])
#                borns.append(born)
#            
#            self._born_charges = np.asarray(borns)
#        
#        return self._born_charges
#
#    def write_born_info(self, filename : str = 'BORNINFO'):
#
#        lines = []
#        
#        ## dielectric tensor
#        for j in range(3):
#            lines.append("%18.13f " * 3 % tuple(self.dielectric_tensor[j]))
#        
#        ## Born effective charge
#        for ia in range(self.number_of_atoms):
#            for j in range(3):
#                lines.append("%18.13f " * 3 % tuple(self.born_charges[ia,j]))
#        
#        lines.append("")
#        f = open(filename, 'w')
#        f.write('\n'.join(lines))
         

