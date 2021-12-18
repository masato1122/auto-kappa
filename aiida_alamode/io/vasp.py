# Copyright (c) aiida-alamode Team.
# Distributed under the terms of the MIT License.
# Written by M. Ohnishi

import numpy as np
import xmltodict
import warnings
from pymatgen.io.vasp.outputs import Vasprun

def write_born_info(filename, outfile='BORNINFO'):
    """
    Args
    -------
    filename (str) : vasprun.xml
    """
    lines = []

    vasprun = Vasprun(filename)
    get_born_charges(filename)
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
        lines = array['set'][i]['v']
        born = np.zeros((3,3))
        for line in lines:
            data = line.split()
            for j in range(3):
                born[j] = float(data[j])
        borns.append(born)
    return borns



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
         

