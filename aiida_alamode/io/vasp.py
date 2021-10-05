# Copyright (c) aiida-alamode Team.
# Distributed under the terms of the MIT License.
# Written by M. Ohnishi

import numpy as np
import xmltodict
import warnings

class Vasprun:
    """ This class reads vasprun.xml.
    """
    def __init__(self, filename):
        
        self.filename = filename
        self._number_of_atoms = None
        self._born_charges = None
        self._dielectric_tensor = None        # dielectric tensor
         
        with open(filename, 'r') as f:
            self._doc = xmltodict.parse(f.read())
    
    @property
    def number_of_atoms(self):
        return self.get_number_of_atoms()
    
    def get_number_of_atoms(self):
        if self._number_of_atoms is None:
            natoms = int(self._doc['modeling']['atominfo']['atoms'])
            self._number_of_atoms = natoms
        return self._number_of_atoms
    
    @property
    def dielectric_tensor(self):
        return self.get_dielectric_tensor()

    def get_dielectric_tensor(self):
        
        if self._dielectric_tensor is not None:
            return self._dielectric_tensor
            
        array = self._doc['modeling']['calculation']['varray']
        
        ## read 'epsilon' in vasprun.xml
        for each in array:
            if each['@name'] == 'epsilon':
                eps = []
                for i, line in enumerate(each['v']):
                    data = line.split()
                    eps.append(np.asarray([float(data[j]) for j in range(3)]))
        self._dielectric_tensor = np.asarray(eps)
        return self._dielectric_tensor
        
    @property
    def born_charges(self):
        return self.get_born_charges()

    def get_born_charges(self):
        
        if self._born_charges is None:
            array = self._doc['modeling']['calculation']['array']
            
            ## Check
            if array['@name'] != 'born_charges':
                warnings.warn(' Cannot find born_charges in %s' % self.filename)
                return None
            
            ## Read contents in "born_charges"
            borns = []
            for i in range(self.number_of_atoms):
                lines = array['set'][i]['v']
                born = np.zeros((3,3))
                for line in lines:
                    data = line.split()
                    for j in range(3):
                        born[j] = float(data[j])
                borns.append(born)
            
            self._born_charges = np.asarray(borns)
        
        return self._born_charges

    def write_born_info(self, filename : str = 'BORNINFO'):

        lines = []
        
        ## dielectric tensor
        for j in range(3):
            lines.append("%18.13f " * 3 % tuple(self.dielectric_tensor[j]))
        
        ## Born effective charge
        for ia in range(self.number_of_atoms):
            for j in range(3):
                lines.append("%18.13f " * 3 % tuple(self.born_charges[ia,j]))
        
        lines.append("")
        f = open(filename, 'w')
        f.write('\n'.join(lines))
         

