# -*- coding: utf-8 -*-
from pymatgen.io.vasp.sets import MPRelaxSet
from pymatgen.core import Structure

s = Structure.from_file("POSCAR.GaSmO")
#print(s)

mpr = MPRelaxSet(s)

mpr.write_input('./out')

