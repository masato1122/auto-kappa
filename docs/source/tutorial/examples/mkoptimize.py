import numpy as np
import ase.io
from auto_kappa.io import AlmInput

filename = './FILES/POSCAR-supercell'
supercell = ase.io.read(filename, format='vasp')
elems = list(set(supercell.get_chemical_symbols()))

dfset = "./FILES/DFSET.cube"
fc2xml = "./FILES/FC2.xml"
borninfo = "./FILES/BORNINFO"

def get_cutoff_for_FC3(elems, cutoff3=4.3):
    cutoff3_bohr = cutoff3 * 1.88973
    cutoffs = {}
    for i1 in range(len(elems)):
        ## It may be safer to start from 0, not from "i1".
        for i2 in range(len(elems)):
            lab = "%s-%s" % (elems[i1], elems[i2])
            cutoffs[lab] = [None, cutoff3_bohr]
    return cutoffs

cutoff_radii = get_cutoff_for_FC3(elems, cutoff3=4.3)

alminp = AlmInput.from_structure(
        supercell,
        mode='optimize',
        norder=2, cutoff=cutoff_radii,
        dfset=dfset,
        fc2xml=fc2xml,
        nonanalytic=1, borninfo=borninfo,
        )

outfile = 'fc3.in'
alminp.to_file(filename=outfile)
print(" Output", outfile)
