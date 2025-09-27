# 
# system.py
# 
# This script parses "system" section of ALAMODE log file.
# 
# Created on September 14, 2025
# Copyright (c) 2025 Masato Ohnishi
#
# This file is distributed under the terms of the MIT license.
# Please see the file 'LICENCE.txt' in the root directory
# or http://opensource.org/licenses/mit-license.php for information.
# 
import numpy as np
import ase

from auto_kappa.units import BohrToA
from auto_kappa.almlog.utils import (
    parse_data_line, 
    replace_symbols_to_blank,
    extract_afterward_lines, 
    get_position_list
)

def _get_vectors(lines, title, idx_init=0, idx_last=2, dtype=float):
    extracted_lines = extract_afterward_lines(lines, title=title, num_lines=3)
    vec_list = []
    for line in extracted_lines:
        parts = line.strip().split()
        try:
            vec_list.append([dtype(v) for v in parts[idx_init:idx_last+1]])
        except ValueError:
            print(f" Error reading {title}.")
            pass
    return np.array(vec_list)

def _get_species_dict(lines, num_species):
    title = "Atomic species"
    extracted_lines = extract_afterward_lines(lines, title=title, num_lines=num_species)
    species_dict = {}
    for line in extracted_lines:
        if not line.strip():
            continue
        line = replace_symbols_to_blank(line)
        parts = line.strip().split()
        try:
            index = parts[0].strip()
            element = parts[1].strip()
            species_dict[index] = element
        except (IndexError, ValueError):
            print(f" Error reading atomic species: {line}")
            pass
    return species_dict

def _get_system_structure(lines):
    
    ## Lattice vectors
    try:
        lattice_vector = _get_vectors(lines, title="Lattice Vector") * BohrToA
    except ValueError:
        lattice_vector = None
    
    ## Atomic positions
    try:
        title = "Atomic positions in fractional basis and atomic species"
        indices, frac_coords, species_numbers = get_position_list(lines, title)
    except TypeError:
        indices, frac_coords, species_numbers = None, None, None
    
    ## Atomic species
    if species_numbers is not None:
        num_species = len(set(species_numbers))
        species_dict = _get_species_dict(lines, num_species)
        chemical_symbols = [species_dict[num] for num in species_numbers]
    try:
        structure = ase.Atoms(
            symbols=chemical_symbols,
            scaled_positions=frac_coords,
            cell=lattice_vector,
            pbc=True
        )
    except (ValueError, TypeError):
        structure = None
        print(" Error creating structure.")
    
    return structure

def read_system(lines):
    """ Read system information from log file lines (System section)
    """
    info = {}
    
    for il, line in enumerate(lines):
        
        line = line.strip().lower()
        if not line:
            continue
        
        if line.startswith("cell volume"):
            info['cell_volume'] = parse_data_line(line, index=-1, dtype=float)
    
    ## Reciprocal lattice vectors
    try:
        info['reciprocal_lattice_vector'] = _get_vectors(lines, title="Reciprocal Lattice Vector")
    except ValueError:
        info['reciprocal_lattice_vector'] = None
    
    ## Get structure
    info['structure'] = _get_system_structure(lines)    
    
    return info
