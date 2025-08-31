from pymatgen.io.vasp import Poscar
from auto_kappa.io import AnphonInput

def main():
    filename = './FILES/POSCAR-supercell'
    fc2xml = "./FILES/FC2.xml"
    borninfo = "./FILES/BORNINFO"
    
    supercell = Poscar.from_file(filename).structure
    primitive = supercell.get_primitive_structure()
    
    anpinp = AnphonInput.from_structure(
        primitive,
        mode='phonons',
        kpmode=1,
        fcsxml=fc2xml,
        nonanalytic=1,
        borninfo=borninfo,
    )
    anpinp.set_kpoint(deltak=0.01)
    
    outfile = 'band.in'
    anpinp.to_file(outfile)
    print(" Output", outfile)

if __name__ == '__main__':
    main()
