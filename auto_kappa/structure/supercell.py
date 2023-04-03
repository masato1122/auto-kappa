#
# supercell.py
#
# This script suggests a supercell based on the given maximum number of atoms.
# Please note that this script is a modified script of cells.py in Phonopy.
# Please cite Phonopy when you use this script.
#
# Copyright (c) 2023 Masato Ohnishi
#
# This file is distributed under the terms of the MIT license.
# Please see the file 'LICENCE.txt' in the root directory
# or http://opensource.org/licenses/mit-license.php for information.
#
import numpy as np
from auto_kappa.structure.crystal import change_structure_format
import spglib

def estimate_supercell_matrix(structure, max_num_atoms=120, max_iter=100):
    
    cell = (
            structure.cell,
            structure.get_scaled_positions(),
            structure.numbers)
    
    spglib_dataset = spglib.get_symmetry_dataset(cell)
    
    spg_num = spglib_dataset["number"]
    num_atoms = len(spglib_dataset["std_types"])
    
    ### in Phonopy, which might be wrong
    #lengths = _get_lattice_parameters(spglib_dataset["std_lattice"])
    ### corrected
    lengths = _get_lattice_parameters(spglib_dataset["std_lattice"].T)
    
    if spg_num <= 74:     # Triclinic, monoclinic, and orthorhombic
        multi = _get_multiplicity_abc(
            num_atoms, lengths, max_num_atoms, max_iter=max_iter
        )
    elif spg_num <= 194:  # Tetragonal and hexagonal
        multi = _get_multiplicity_ac(
            num_atoms, lengths, max_num_atoms, max_iter=max_iter
        )
    else:  # Cubic
        multi = _get_multiplicity_a(
            num_atoms, lengths, max_num_atoms, max_iter=max_iter
        )

    return np.diag(multi)

def _get_lattice_parameters(lattice):
    """Return basis vector lengths.

    Parameters
    ----------
    lattice : array_like
        Basis vectors given as column vectors
        shape=(3, 3), dtype='double'

    Returns
    -------
    ndarray, shape=(3,), dtype='double'

    """
    return np.array(np.sqrt(np.dot(lattice.T, lattice).diagonal()), dtype="double")

#def get_sphericity(lattice):
#    """
#    Args
#    -----
#    lattice : float, shape=(3,3)
#        row vectors
#    """
#    ### total surface area
#    Atot = 0.
#    for i in range(3):
#        i1 = i
#        i2 = (i1+1) % 3
#        v1 = lattice[i1]
#        v2 = lattice[i2]
#        area = np.linalg.norm(np.cross(v1, v2))
#        Atot += area * 2.0
#
#    ### calculate the volume
#    volume = abs(np.dot(
#        np.cross(lattice[0], lattice[1]), 
#        lattice[2]
#        ))
#
#    ### calculate the sphericity
#    chi = np.power(np.pi, 1./3.) * np.power(6.*volume, 2./3.) / Atot
#    return chi

#def _get_next_modification_index(lengths, multi):
#    """ Return the modified index
#
#    Args
#    -------
#    lengths : float, shape=(n)
#    multi : int, shape=(n)
#
#    Return
#    -------
#    mod_index : int, (< n)
#    """
#    stds = []
#    for imod in range(len(lengths)):
#        multi_next = multi.copy()
#        multi_next[imod] += 1
#        lengths_sc = np.asarray(lengths) * np.asarray(multi_next)
#        stds.append(np.std(lengths_sc))
#    return np.argmin(stds)    

def _get_multiplicity_abc(num_atoms, lengths, max_num_atoms, max_iter=20):
    multi = [1, 1, 1]

    for i in range(max_iter):
        l_super = np.multiply(multi, lengths)
        
        ### ver.1 original (Phonopy)
        min_index = np.argmin(l_super)
        multi[min_index] += 1
        if num_atoms * np.prod(multi) > max_num_atoms:
            multi[min_index] -= 1
        
        #### ver.2 modified
        #mod_index = _get_next_modification_index(lengths, multi)
        #multi[mod_index] += 1
        #if num_atoms * np.prod(multi) > max_num_atoms:
        #    multi[mod_index] -= 1

    return multi


def _get_multiplicity_ac(num_atoms, lengths, max_num_atoms, max_iter=20):
    multi = [1, 1]
    a = lengths[0]
    c = lengths[2]

    for i in range(max_iter):
        l_super = np.multiply(multi, [a, c])
        
        #### ver.1 Phonopy
        min_index = np.argmin(l_super)
        multi[min_index] += 1
        if num_atoms * multi[0] ** 2 * multi[1] > max_num_atoms:
            multi[min_index] -= 1
        
        ### ver.2 modified
        #mod_index = _get_next_modification_index([a,c], multi)
        #multi[mod_index] += 1
        #if num_atoms * multi[0] ** 2 * multi[1] > max_num_atoms:
        #    multi[mod_index] -= 1

    return [multi[0], multi[0], multi[1]]


def _get_multiplicity_a(num_atoms, lengths, max_num_atoms, max_iter=20):
    multi = 1
    for i in range(max_iter):
        multi += 1
        if num_atoms * multi**3 > max_num_atoms:
            multi -= 1

    return [multi, multi, multi]

