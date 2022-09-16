# -*- coding: utf-8 -*-
import numpy as np
from optparse import OptionParser

import ase.io
from ase.build import make_supercell
from auto_alamode2.apdb import ApdbVasp
from auto_alamode2.structure.crystal import get_standerdized_cell
    
parser = OptionParser()
parser.add_option("-f", "--filename", dest="filename", type="string",
        default="POSCAR.orig", help="input file name")
(options, args) = parser.parse_args()

### prepare structures
prim = ase.io.read(options.filename, format='vasp')
unit = get_standerdized_cell(prim)

### get the primitive matrix
pmat = np.dot(np.linalg.inv(unit.cell), prim.cell)

### set ApdbVasp
apdb = ApdbVasp(unit, primitive_matrix=pmat)


### calculator
kpts = [4, 4, 4]
calc = apdb.get_calculator('relax', './out', kpts)

###
calc.write_input(prim)

