from aiida_alamode.io import AnphonInput

filename = 'POSCAR-supercell'
fcsxml = 'IFCs.xml'

## to calculate band
anpinp = AnphonInput.from_structure_file(
        filename,
        mode='phonons',
        kpmode=1,
        fcsxml=fcsxml
        )
anpinp.set_kpoint(deltak=0.01)
anpinp.to_file()

## to calculate DOS
anpinp = AnphonInput.from_structure_file(
        filename,
        mode='phonons',
        kpmode=2,
        fcsxml=fcsxml
        )
anpinp.set_kpoint(deltak=0.1)
#anpinp.update({'kpts':[10, 10, 10]})
anpinp.to_file()

## to calculate thermal conductivity
anpinp = AnphonInput.from_structure_file(
        filename,
        mode='RTA',
        kpmode=2,
        fcsxml=fcsxml
        )
anpinp.set_kpoint()
anpinp.to_file()
