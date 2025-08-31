from pymatgen.io.vasp import Poscar, Kpoints
from auto_kappa.io import AnphonInput

def main():
    filename = './FILES/POSCAR-supercell'
    fcsxml = "./FILES/FC3_lasso.xml"
    borninfo = "./FILES/BORNINFO"
    
    structure = Poscar.from_file(filename).structure
    primitive = structure.get_primitive_structure()
    
    kpts = Kpoints.automatic_density(primitive, 1000).kpts[0]
    
    anpinp = AnphonInput.from_structure(
        primitive,
        mode='rta',
        kpmode=2,
        fcsxml=fcsxml,
        nonanalytic=1,
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
