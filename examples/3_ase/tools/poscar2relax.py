# -*- coding: utf-8 -*-
import os
import numpy as np
from optparse import OptionParser

import ase.io
from ase.build import make_supercell
from auto_kappa.apdb import ApdbVasp
from auto_kappa.structure.crystal import (
        get_primitive_structure_spglib,
        get_standardized_structure,
        )

parser = OptionParser()
parser.add_option("-f", "--filename", dest="filename", type="string",
        default="POSCAR.unitcell", help="input file name")
parser.add_option("-o", "--outdir", dest="outdir", type="string",
        default="./out", help="output directory name")
(options, args) = parser.parse_args()

### prepare structures
unit = ase.io.read(options.filename, format='vasp')

prim = get_primitive_structure_spglib(unit)

### get the primitive matrix
pmat = np.dot(np.linalg.inv(unit.cell), prim.cell)

### set ApdbVasp
apdb = ApdbVasp(unit, primitive_matrix=pmat, scell_matrix=np.identity(3))

### calculator
os.makedirs(options.outdir, exist_ok=True)
kpts = [4, 4, 4]
calc = apdb.get_calculator('relax', options.outdir, kpts)

###
calc.write_input(prim)

