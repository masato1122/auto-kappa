#!/usr/bin/env python
# File:    mkkappa.py
import ase.io
from auto_kappa.io import AnphonInput
from auto_kappa.structure import get_primitive_structure_spglib
from auto_kappa.cui.suggest import klength2mesh

def main():
    filename = './FILES/POSCAR-supercell'
    fcsxml = "./FILES/FC3_lasso.xml"
    borninfo = "./FILES/BORNINFO"
    k_length = 20
    
    structure = ase.io.read(filename, format='vasp')
    primitive = get_primitive_structure_spglib(structure)
    
    kpts = klength2mesh(k_length, primitive.cell.array)
    
    anpinp = AnphonInput.from_structure(
        primitive,
        mode='rta',
        kpmode=2,
        fcsxml=fcsxml,
        nonanalytic=2,
        borninfo=borninfo,
        kpts=kpts,
        dt=50,
        tmin=50,
        tmax=1000,
        kappa_coherent=1,
    )
    
    outfile = 'kappa.in'
    anpinp.to_file(outfile)
    print(" Output", outfile)

if __name__ == '__main__':
    main()
