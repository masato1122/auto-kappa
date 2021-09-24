# -*- coding: utf-8 -*-
import numpy as np
from optparse import OptionParser

from aiida_alamode.io import AlmInput, AnphonInput
#get_alm_variables

filename = 'POSCAR-unitcell'

#params = get_alm_variables(filename=filename, guess_primitive=True)
#print(params)

alminp = AlmInput(filename=filename)

alminp.from_structure(filename)

print(alminp.set_prefix())

alminp.to_file(filename='alm.in')


#for d in dir(alminp):
#    print(d)


