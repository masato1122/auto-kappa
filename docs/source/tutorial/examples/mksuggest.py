import numpy as np
import ase.io
from auto_kappa.io import AlmInput

filename = './FILES/POSCAR-supercell'
supercell = ase.io.read(filename, format='vasp')
elems = list(set(supercell.get_chemical_symbols()))

### 1. suggest for FC2 (norder=1)
alminp = AlmInput.from_structure(
        supercell,
        mode='suggest',
        norder=1,
        )

outfile = "suggest_fc2.in"
alminp.to_file(filename=outfile)
print(" Output", outfile)

## 2. suggest for FC3 (norder=2)
## Note that AlmInput required the length unit of Bohr while auto-kappa required
## the length unit of Angstrom.
norder = 2
alminp.update({'norder': norder})
alminp.update({'nbody': [2, 3]})

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
alminp.update({'cutoff': cutoff_radii})

outfile = "suggest_fc3.in"
alminp.to_file(filename=outfile)
print(" Output", outfile)

