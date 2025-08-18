# 
# suggest.py
#
# This script parses files for "suggest" mode
#
# Created on August 15, 2025
# Copyright (c) 2025 Masato Ohnishi
#
# This file is distributed under the terms of the MIT license.
# Please see the file 'LICENCE.txt' in the root directory
# or http://opensource.org/licenses/mit-license.php for information.
# 
import os
# import sys
import numpy as np
import re

import ase
# from ase.geometry import get_distances
# from auto_kappa.io.alm import read_structure_from_file
from auto_kappa.structure import change_structure_format
from auto_kappa.units import BohrToA

import logging
logger = logging.getLogger(__name__)

class Pattern:
    """ Class to handle ".pattern_{HARMONIC/ANHARM*}" files
    
    How to use
    -----------
    >>> from auto_kappa.io.alm import read_structure_from_file
    >>> structure = read_structure_from_file("suggest.in")
    >>> from auto_kappa.io.suggest import Pattern
    >>> pattern = Pattern('Si.pattern_ANHARM3')
    >>> pattern.set_structure(structure)
    >>> mag = 0.03
    >>> pattern.get_displacements(mag=mag)
    >>> pattern.get_suggested_structures(mag=mag)
    """
    def __init__(self, filename, structure=None):
        self._filename = filename
        self.basis = None
        self._patterns = None
        self._structure = structure
        if os.path.exists(filename):
            self.set_pattern_info()
    
    @property
    def filename(self):
        return self._filename
    @property
    def patterns(self):
        if self._patterns is None:
            self.set_pattern_info()
        return self._patterns
    @property
    def structure(self):
        """ Pristine structure (supercell) """
        if self._structure is None:
            raise ValueError("\n Structure is not set. Please provide a structure.")
        return self._structure
    
    def set_structure(self, structure):
        """ Set the pristine structure (supercell) """
        self._structure = structure
    
    def set_pattern_info(self):
        self.basis, self._patterns = _parse_pattern_file(self.filename)
    
    def print_patterns(self):
        logger.info(f"\n Basis: {self.basis}")
        for pattern in self.patterns:
            msg = f" {pattern['index']:2d} | order {pattern['order']} | "
            for i, (atom_index, direction) in enumerate(pattern['data']):
                msg += f"{atom_index:2d}: (" + "%2d " * 3 % tuple(direction) + ")"
                if i < len(pattern['data']) - 1:
                    msg += ",  "
            logger.info(msg)

    def get_distances(self):
        """ Get distances between displaced atoms. When multiple atoms are displaced,
        the maximum atomic distance is used.
        """
        if self.structure is None:
            logger.error("\n Error: Structure is not given.")
            return None
        
        if not isinstance(self.structure, ase.Atoms):
            structure = change_structure_format(self.structure, format='ase')
        else:
            structure = self.structure.copy()
        
        Ds = structure.get_all_distances(mic=True) # distances
        
        ds = []
        for pattern in self.patterns:
            # index = pattern['index']
            order = pattern['order']
            data = pattern['data']
            if order > 1:
                indices = []
                for each in data:
                    indices.append(each[0] - 1) # atomic index
                dmax = _get_maximum_distance(Ds, indices)
                ds.append(dmax)
        
        return np.asarray(ds, dtype=float)

    def get_displacements(self, mag=None):
        """ Get displacement patterns from the structure. 
        Args
        ----
        mag : float
            Magnitude of the displacement with the unit of Angstrom.
        """
        if self.structure is None:
            msg = "\n Error: Structure is not given."
            logger.error(msg)
            return None
        
        if mag is None:
            msg = "\n Error: Displacement magnitude is not given."
            logger.error(msg)
            return None
        
        natoms = len(self.structure)
        disp_all = []
        for pattern in self.patterns:
            ds = np.zeros((natoms, 3), dtype=float)
            for atom_index, direction in pattern['data']:
                ## displace the "atom_index"-th atom
                ds[atom_index-1] += direction * mag
            disp_all.append(ds)
        return disp_all

    def get_suggested_structures(self, mag=None):
        """ Get suggested structures from the displacement patterns.
        Args
        ----
        mag : float
            Magnitude of the displacement with the unit of Angstrom.
        """
        if mag is None:
            msg = "\n Error: Displacement magnitude is not given."
            logger.error(msg)
            return None
        
        ## Pristine structure
        prist = change_structure_format(self.structure, format='ase')

        ## Make structures with displacement patterns
        structures = []
        disp_all = self.get_displacements(mag=mag)
        for ds in disp_all:
            disp_atoms = prist.copy()
            mod_positions = prist.get_positions() + ds
            disp_atoms.set_positions(mod_positions)
            structures.append(disp_atoms)
        return structures

def _get_maximum_distance(distances, atomic_indices):
    """ Get the maximum distance between atoms in the given indices. 
    
    Args
    -----
    distances : np.ndarray, shape=(natoms, natoms)
        Array of distances between atoms.
    atomic_indices : list, shape=(n)
        List of atomic indices to consider.
    """
    d_pairs = []
    n = len(atomic_indices)
    for i1 in range(n):
        iatom1 = atomic_indices[i1]
        for i2 in range(i1 + 1, n):
            iatom2 = atomic_indices[i2]
            d_pairs.append(distances[iatom1, iatom2])
    return np.max(d_pairs)

def _parse_pattern_file(filename):
    
    basis = None
    patterns = []
    
    with open(filename, "r", encoding="utf-8") as f:
        for line in f:
            
            line = line.strip()
            if not line:
                continue  # Skip empty lines
            # print(line)
            
            # line "Basis : C"
            if line.startswith("Basis"):
                parts = line.split(":")
                basis = parts[-1]
            
            # line "n:    m"
            elif ":" in line:
                parts = line.split(":")
                index = int(parts[0].strip())
                order = int(parts[1].strip())
                block_data = []
                n_lines_read = 0
                # if order not in patterns:
                #     patterns[order] = []
            
            # Data line
            else:
                values = list(map(int, line.split()))
                atom_index = int(values[0])
                direction = np.array(values[1:], dtype=float)
                block_data.append((atom_index, direction))
                n_lines_read += 1
                
                # Save the block data
                if n_lines_read == order:
                    patterns.append({'index': index, 'order': order, 'data': block_data})
                    block_data = []
    
    return basis, patterns

class SuggestLogParser:
    """ Parse the log file for "suggest" mode (suggest.log)
    
    How to use
    -----------
    >>> filename = 'suggest.log'
    >>> suggest_log = SuggestLogParser(filename)
    >>> suggest_log.print_neighboring_info()
    >>> print(suggest_log.get_neighboring_distances())
    >>> suggest_log.structure
    """
    def __init__(self, filename):
        self._filename = filename
        
        if os.path.exists(filename):
            self.set_neighboring_info()
            self._contents = parse_suggest_log(filename)
        else:
            self._neighboring_info = None
            self._contents = None
    
    @property
    def filename(self):
        return self._filename
    @property
    def contents(self):
        return self._contents
    @property
    def structure(self):
        return self.contents.get('structure', None)
    @property
    def space_group(self):
        return self.contents.get('space_group', None)
    
    def set_neighboring_info(self):
        self._neighboring_info = parse_neighbor_data(self.filename)
    
    def get_neighboring_info(self):
        if self._neighboring_info is None:
            self.set_neighboring_info()
        return self._neighboring_info

    def print_neighboring_info(self, unit='Ang'):
        for entry in self._neighboring_info:
            msg = "%3d (%2s): " % (entry['center_atom']['index'], entry['center_atom']['name'])
            d = entry['distance']
            if unit.upper().startswith('A'):
                d *= BohrToA
            msg += f"{entry['shell_index']}-th, {d:7.3f} {unit} | "
            for i in range(entry['count']):
                msg += f"{entry['neighbors'][i]['index']:2d} ({entry['neighbors'][i]['name']:2s})"
                if i < entry['count'] - 1:
                    msg += ", "
            logger.info(msg)
        
    def get_neighboring_distances(self, unit='Ang'):
        """ Get distances between neighboring atoms with the unit of Angstrom
        
        Args
        ----
        unit : string
            "Angstrom" or "Bohr"
        
        Return
        -------
        dict
            A dictionary mapping atom names to their neighboring distances.
        """
        neighboring_distances = {}
        for entry in self.get_neighboring_info():
            name = entry['center_atom']['name']
            if name not in neighboring_distances:
                neighboring_distances[name] = []
            d = entry['distance']  # Bohr
            if unit.upper().startswith('A'):
                d *= BohrToA
            neighboring_distances[name].append(d)
        return neighboring_distances


def _get_title_line_index(lines, title):
    for il, line in enumerate(lines):
        if line.strip().lower() == title.lower():
            return il
    return None

def _get_lattice_vector(lines):
    """ Get the lattice vector from the log file. 
    
    Return
    ------
    np.ndarray
        The lattice vector in Angstrom.
    """
    il0 = _get_title_line_index(lines, "Lattice Vector")
    if il0 is None:
        raise ValueError("\n Lattice vector section not found in the log file.")
    
    vectors = []
    for i in range(3):
        parts = lines[il0 + 1 + i].strip().split()
        for j in range(3):
            vectors.append(float(parts[j]))
    return np.array(vectors).reshape(3, 3) * BohrToA

def _get_atomic_species(lines):
    """ Get atomic species """
    il0 = _get_title_line_index(lines, "Atomic Species:")
    if il0 is None:
        raise ValueError("\n Atomic species section not found in the log file.")

    species = {}
    flag = True
    count = 0
    while flag:
        parts = lines[il0 + 1 + count].strip().split()
        if len(parts) != 2:
            break
        species[int(parts[0])] = str(parts[1])
        count += 1
    return species

def _get_scaled_positions(lines):
    """ Get atomic positions """
    il0 = _get_title_line_index(lines, "Atomic positions in fractional basis and atomic species")
    if il0 is None:
        raise ValueError("\n Atomic positions section not found in the log file.")

    positions = []
    species_index = []
    flag = True
    count = 0
    while flag:
        parts = lines[il0 + 1 + count].strip().split()
        if len(parts) != 5:
            break
        positions.append([])
        for j in range(3):
            positions[-1].append(float(parts[j + 1]))
        species_index.append(int(parts[-1]))
        count += 1
    return positions, species_index

def _get_space_group_info(lines):
    """ Get the space group from the log file. 
    
    Return
    ------
    str
        The space group.
    """
    for il in range(len(lines)):
        line = lines[il].strip()
        if line.startswith("Space group:"):
            line = line.replace("(", " ").replace(")", " ")
            parts = line.split()
            mat_sg = (parts[-2], int(parts[-1]))
            return mat_sg
    return None

def parse_suggest_log(filename):
    
    with open(filename, 'r') as f:
        lines = f.readlines()
    
    contents = {}
    
    ### SYSTEM
    contents['lattice_vector'] = _get_lattice_vector(lines) # Angstrom
    contents['atomic_species'] = _get_atomic_species(lines)
    contents['scaled_positions'], contents['species_index'] = _get_scaled_positions(lines)
    
    ## Make an ase.Atoms object
    flag = True
    for key in contents:
        if contents[key] is None:
            flag = False
    if flag:
        atomic_species = contents['atomic_species']
        symbols = [atomic_species[i] for i in contents['species_index']]
        atoms = ase.Atoms(
            cell=contents['lattice_vector'],
            scaled_positions=contents['scaled_positions'],
            symbols=symbols
        )
        contents['structure'] = atoms
    
    ### SYMMETRY
    contents['space_group'] = _get_space_group_info(lines)
    
    return contents


def _parse_current_atoms(line):
    """ Parse the current atom line to extract the atom type and its index.
    Example: "1 ( Te)"
    """
    parts = re.sub(r'[()]', ' ', line).split()
    atoms = {"index": int(parts[0]), "name": str(parts[1])}
    return atoms

def _parse_distance(line):
    """ Parse the distance line to extract the distance value.
    Example: "1  1.234 ( 2 )"
    """
    parts = re.sub(r'[()]', ' ', line).split()
    distance = {"index": int(parts[0]), "distance": float(parts[1]), "count": int(parts[2])}
    return distance

def _parse_neighbors(line):
    """ Parse the neighbors line to extract the neighbor atoms.
    Example: "1 ( Te) 2 ( Pb) ..."
    """
    parts = re.sub(r'[()]', ' ', line).split()
    neighbors = []
    for i in range(0, len(parts), 2):
        neighbor = {"index": int(parts[i]), "name": str(parts[i + 1])}
        neighbors.append(neighbor)
    return neighbors

def parse_neighbor_data(filename):
    
    with open(filename, 'r') as f:
        lines = f.readlines()
    
    line_start = "Format [N th-nearest shell, distance (Number of atoms on the shell)]"
    line_end = "Time Elapsed:"
    
    pattern_elem_in_braket = r'\(\s*\w+\s*\)' # "( Pb )", "(Te)", ...
    pattern_title = r'\d+\s*' + pattern_elem_in_braket
    pattern_subtitle = r'\d+\s*\d+\.\d+\s*' + pattern_elem_in_braket + r'\s*-\s*'
    pattern_repeat = r'\d+\s*(\(\s*\w+\s*\)\s*)+'
    pattern1 = pattern_title + r'\s*:\s*' + pattern_subtitle + pattern_repeat
    pattern2 = pattern_subtitle + pattern_repeat
    pattern3 = pattern_repeat
    
    results = []
    current_atom = None
    
    flag = False
    for line in lines:
        
        if line_start in line:
            flag = True
            continue
        if not flag:
            continue
        if line_end in line:
            break
        
        line = line.strip()
        
        # Start of a new central atom
        if re.match(pattern1, line):
            line_current_atom = line.split(':')[0].strip()
            line_distance = line.split(':')[1].split('-')[0].strip()
            line_neighbors = line.split(':')[1].split('-')[1].strip()
            current_atom = _parse_current_atoms(line_current_atom)
            dist_info = _parse_distance(line_distance)
            neighbors = _parse_neighbors(line_neighbors)

            results.append({
                "center_atom": current_atom,
                "shell_index": dist_info["index"],
                "distance": dist_info["distance"], # Bohr
                "count": dist_info['count'],
                "neighbors": neighbors
            })
        
        # Line that starts a new shell entry
        elif re.match(pattern2, line):
            
            line_distance = line.split('-')[0].strip()
            line_neighbors = line.split('-')[1].strip()
            
            dist_info = _parse_distance(line_distance)
            neighbors = _parse_neighbors(line_neighbors)

            results.append({
                "center_atom": current_atom,
                "shell_index": dist_info["index"],
                "distance": dist_info["distance"],
                "count": dist_info['count'],
                "neighbors": neighbors
            })
            
        # Continuation of previous neighbor line
        elif re.match(pattern3, line):
            neighbors = _parse_neighbors(line)
            results[-1]['neighbors'].extend(neighbors)
        # else:
        #     if len(line.split()) > 0:
        #         raise ValueError(" Unexpected line format: {}".format(line))
    
    return results


# def main(options):
    
#     file_log = 'suggest.log'
#     suggest_log = SuggestLogParser(file_log)
#     # print(suggest_log.get_neighboring_distances())
    
#     input_file = 'suggest.in'
#     structure = read_structure_from_file(input_file)
    
#     patterns = Pattern("TePb.pattern_ANHARM3", structure=structure)
#     patterns.print_patterns()
    

# if __name__ == '__main__':
    
#     parser = argparse.ArgumentParser(description='Input parameters')
    
#     parser.add_argument('-f', '--filename', dest='filename', type=str,
#                         default="./suggest.log", help="input file name")
    
#     args = parser.parse_args()

#     main(args)
