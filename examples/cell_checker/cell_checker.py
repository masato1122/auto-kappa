# -*- coding: utf-8 -*-
import os
import os.path
import numpy as np
from optparse import OptionParser

from auto_kappa.structure.crystal import change_structure_format
import yaml

import ase
import ase.io
from phonopy.structure.cells import get_supercell

def read_phonondb(dir_phdb):
    """ Read data in Phonondb and return unitcell and primitive and supercell
    matrices.
    """
    from auto_kappa.io.phonondb import Phonondb
    
    ### unitcell
    filename = dir_phdb + '/POSCAR-unitcell'
    atoms = ase.io.read(filename, format='vasp')
    
    ### phonondb
    phdb = Phonondb(options.dir_phdb)
    
    unit_ase = phdb.get_unitcell(format='ase')
    mat_u2p = phdb.primitive_matrix
    mat_u2s = phdb.scell_matrix
    
    return unit_ase, mat_u2p, mat_u2s

def make_structures(unit, prim_mat, sc_mat, format='ase'):
    """ Make and return primitive and supercells
    
    unit :
        unitcell
    
    prim_mat : float, shape=(3,3)
        primitive matrix

    sc_mat : float, shape=(3,3)
        supercell matrix
    """
    from phonopy import Phonopy
    
    unit_pp = change_structure_format(unit, format='phonopy')
    phonon = Phonopy(unit_pp, sc_mat, primitive_matrix=prim_mat)
    
    prim = change_structure_format(phonon.primitive, format=format)
    sc = change_structure_format(phonon.supercell, format=format)
    structures = {'unitcell': unit, "primitive": prim, "supercell": sc}
    
    return structures

def _get_relaxed_structure(filename, prim_mat):
    """ Make and return the unitcell

    Args
    ----
    filename : string
        ./relax/relax.yaml
    
    prim_mat : float, shape=(3,3)
        primitive matrix
    """

    with open(filename) as file:
        obj = yaml.safe_load(file)
    
    lattice = np.asarray(obj['lattice'])
    scaled_positions = np.asarray(obj['positions'])
    symbols = obj['species']
    cell_type = obj['cell_type_for_relaxation']
    
    positions = scaled_positions @ lattice
    
    ### Atoms obj
    atoms = ase.Atoms(cell=lattice, pbc=True)
    for (pos, el) in zip(positions, symbols):
        atoms.append(ase.Atom(position=pos, symbol=el))
    
    ###
    if cell_type.lower()[0] == "p":
        mat_p2u = np.linalg.inv(prim_mat)
        prim_pp = change_structure_format(atoms, format='phonopy')
        unit = change_structure_format(
                get_supercell(prim_pp, mat_p2u),
                format='ase')
        
    else:
        unit = atoms.copy()
    
    from ase.build.tools import sort
    unit_sort = sort(unit, tags=unit.get_chemical_symbols())
    return unit_sort

def get_structures_apdb(dir_apdb, prim_mat):
    """ Read and return structures from calculated results
    """
    ### relaxed structure
    filename = dir_apdb + "/relax/relax.yaml"
    try:
        unitcell = _get_relaxed_structure(filename, prim_mat)
    except Exception:
        try:
            filename = dir_apdb + "/relax/CONTCAR"
            prim = ase.io.read(filename, format='vasp')
            prim_pp = change_structure_format(prim, format='phonopy')
            unitcell = change_structure_format(
                    get_supercell(prim_pp, np.linalg.inv(prim_mat)),
                    format="ase"
                    )
        except Exception:
            try:
                filename = dir_apdb + "/relax/full-2/CONTCAR"
                prim = ase.io.read(filename, format='vasp')
                prim_pp = change_structure_format(prim, format='phonopy')
                unitcell = change_structure_format(
                        get_supercell(prim_pp, np.linalg.inv(prim_mat)),
                        format="ase"
                        )
            except Exception:
                unitcell = None
    
    ### supercell
    filename = dir_apdb + "/harm/force/prist/POSCAR"
    try:
        supercell = ase.io.read(filename, format='vasp')
    except Exception:
        supercell = None
    
    structures = {'unitcell': unitcell, "supercell": supercell}

    return structures

def _get_lengths_of_translational_vectors(cell):
    """ Calculate and return the lengths and angles of the translational
    vectors.

    Args
    -----
    cell : float, shpae=(3,3)
    """
    lengths = np.zeros(3)
    for i in range(3):
        lengths[i] = np.linalg.norm(cell[i])
    
    ###
    def get_angle(a, b):
        angle = np.arccos(np.dot(a, b) / (np.linalg.norm(a) * np.linalg.norm(b)))
        return angle * 180/np.pi

    angles = np.zeros(3)
    angles[0] = get_angle(cell[0], cell[1])
    angles[1] = get_angle(cell[1], cell[2])
    angles[2] = get_angle(cell[2], cell[0])
    return lengths, angles

def _print_cell(cell):
    for j in range(3):
        print(" %7.2f" * 3 % tuple(cell[j]), end=" ")
    print("")

def _sort_cell(cell):
    lengths, _ = _get_lengths_of_translational_vectors(cell)
    i0 = np.argmin(lengths)
    i2 = np.argmax(lengths)
    if i0 == i2:
        i0 = 0
        i1 = 1
        i2 = 2
    else:
        for j in range(3):
            if j != i0 and j != i2:
                i1 = j
    cell_new = np.zeros((3,3))
    for inew, ii in enumerate([i0, i1, i2]):
        cell_new[inew] = cell[ii]
    return cell_new

def _compare_cells(cell1, cell2, tol_len=0.1, tol_ang=1.0):
    
    cell1 = _sort_cell(cell1)
    cell2 = _sort_cell(cell2)
    
    lengths1, angles1 = _get_lengths_of_translational_vectors(cell1)
    lengths2, angles2 = _get_lengths_of_translational_vectors(cell2)
    
    ### check lengths
    slen1 = lengths1 / lengths1[0]
    slen2 = lengths2 / lengths2[0]
    
    diff_max = np.max(abs(slen1 - slen2))
    
    flag = 0
    if diff_max > tol_len:
        print(" Error: lengths are different.", diff_max)
        print("  Correct : ", end="")
        #_print_cell(cell1)
        print(" %7.2f " * 3 % tuple(slen1))
        print("  Used    : ", end="")
        print(" %7.2f " * 3 % tuple(slen2))
        #_print_cell(cell2)
        flag += 1
     
    ### check angles
    diff_max = np.max(abs(angles1 - angles2))
    if diff_max > tol_ang:
        print(" Error: angles are different.", diff_max)
        print("  Correct : ", end="")
        print(" %7.2f " * 3 % tuple(angles1))
        print("  Used    : ", end="")
        print(" %7.2f " * 3 % tuple(angles2))
        flag += 10
    
    return flag

def compare_structures(structures_phdb, structures_apdb):
    
    flags = {}
    for type in ['unitcell', 'supercell']:
        
        unit1 = structures_phdb[type]
        unit2 = structures_apdb[type]

        flags[type] = 0
        
        if unit2 is None:
            if type == 'unitcell':
                print(" No relaxed structure")
            else:
                print(" No supercell")
            continue
        
        else:
            flags[type] += _compare_cells(
                    unit1.cell.array,
                    unit2.cell.array
                    )

    return flags

def main(options):
    
    unit, prim_mat, sc_mat = read_phonondb(options.dir_phdb)    
    
    structures_phdb = make_structures(unit, prim_mat, sc_mat)
    
    structures_apdb = get_structures_apdb(options.dir_apdb, prim_mat)  
    
    flags = compare_structures(structures_phdb, structures_apdb)
    
    import shutil
    for type in ['supercell', 'unitcell']:
        
        ### If thermal conductivity has been calculated.
        dir_cube = options.dir_apdb + "/cube"
        if os.path.exists(dir_cube) == False:
            outdir = "WORNGCELL"
        else:
            outdir = "WORNGCELL2"
        
        ###
        flag = flags[type]
        if flag != 0:
            
            print(" Error:", flag)    
            #if flag == 1:
            #    outdir += "/%s/length" % type
            #elif flag == 10 or flag == 20:
            #    outdir += "/%s/angle" % type
            #elif flag == 11 or flag == 21:
            #    outdir += "/%s/length_angle" % type
            #else:
            #    outdir += "/others"
            
            ### move directory
            os.makedirs(outdir, exist_ok=True)
            shutil.move(options.dir_apdb, outdir)
            break

if __name__ == '__main__':
    parser = OptionParser()
    
    parser.add_option("--dir_phdb", dest="dir_phdb", type="string",
            default=None, help="directory for Phonondb")
    
    parser.add_option("--dir_apdb", dest="dir_apdb", type="string",
            default="./mp-2552", 
            help="directory for anharmonic phonon database")
    
    (options, args) = parser.parse_args()
    
    main(options)

