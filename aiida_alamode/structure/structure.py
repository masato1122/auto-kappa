# -*- coding: utf-8 -*-
import numpy as np
from alm import ALM

from phonopy.structure.cells import get_smallest_vectors


def analyze_cutoff(supercell, primitive, symprec=1e-3):
    


    print("AAAAAAAA")
    exit()



    #print(supercell) 
    #for d in dir(supercell):
    #    print(d)
    
    species = supercell.species
    #print(type(species[0]))
    #for d in dir(species[0]):
    #    print(d)
    
    #from pymatgen.core.bonds import get_bond_length
    #get_bond_length(species[0], species[1])
    
    #for d in dir(supercell):
    #    print(d)

    from pymatgen.core.structure import Structure
    from pymatgen.symmetry.analyzer import SpacegroupAnalyzer 
    from pymatgen.symmetry.site_symmetries import get_site_symmetries
    #import (
    #        get_symmetrized_structure)
    
    structure = Structure(
            supercell.lattice.matrix,
            supercell.species,
            supercell.frac_coords
            )
    
    #print(structure.sites[0])
    #exit()
    
    #SpacegroupAnalyzer(supercell)
    pointops = get_site_symmetries(structure)
    print(pointops[0][0])
    print(pointops[0][1])
    print("")
    print(pointops[1][0])
    print(pointops[1][1])
    #print(len(supercell))
    #print(len(pointops[0]))
    
    #for ia in range(len(supercell)):
    #    print(len(pointops[ia]))

    #for d in dir(pointops):
    #    print(d)

    exit()

    pcell = primitive.lattice.matrix
    
    scaled_pos = np.dot(np.linalg.inv(pcell), supercell.cart_coords.T).T
    idx_prim = []
    for ia, spos in enumerate(scaled_pos):
        if -1. <= np.min(spos) and np.max(spos) <= 1.:
            idx_prim.append(ia)
    
    print(idx_prim)
    print(scaled_pos[idx_prim])
     
    get_smallest_vector_of_atom_pair(
            np.arange(len(supercell)), idx_prim,
            supercell, symprec=symprec)
    
    print(len(supercell), len(idx_prim))


#def get_smallest_vector_of_atom_pair(
#        idx_scell, idx_prim, supercell, symprec=1e-3):
#    """Return smallest vectors of an atom pair in supercell."""
#    
#    s_pos = supercell.frac_coords    
#    
#    svecs, multi = get_smallest_vectors(
#            supercell.lattice.matrix,
#            [s_pos[idx_scell]], 
#            [s_pos[idx_prim]], 
#            store_dense_svecs=True, 
#            symprec=symprec
#            )
#    
#    print(svecs, multi)
#    return svecs[0]

