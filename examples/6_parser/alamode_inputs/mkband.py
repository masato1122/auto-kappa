#!/usr/bin/env python
# File:    mkband.py
import ase.io
from auto_kappa.structure import get_primitive_structure_spglib
from auto_kappa.io import AnphonInput

import warnings
warnings.filterwarnings("ignore", category=DeprecationWarning)

def main():
    
    filename = './FILES/POSCAR-supercell'
    fc2xml = "./FILES/FC2.xml"
    borninfo = "./FILES/BORNINFO"
    
    supercell = ase.io.read(filename, format='vasp')
    primitive = get_primitive_structure_spglib(supercell)
    
    anpinp = AnphonInput.from_structure(
        primitive,
        mode='phonons',
        kpmode=1,
        fcsxml=fc2xml,
        nonanalytic=2,
        borninfo=borninfo,
    )
    anpinp.set_kpoint(deltak=0.01)
    
    outfile = 'band.in'
    anpinp.to_file(outfile)
    print(" Output", outfile)

if __name__ == '__main__':
    main()
