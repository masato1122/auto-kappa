from aiida_alamode.io import AlmInput

filename = 'POSCAR-supercell'
dfset = 'dfset.txt'

## mode "suggest"
alminp = AlmInput.from_structure_file(
        filename,
        mode='suggest',
        norder=2
        )
#alminp.to_file()
#alminp.to_file(filename='filename.in')

print(alminp.as_dict())
exit()

from aiida_alamode.structure.structure import analyze_cutoff
analyze_cutoff(alminp.structure, alminp.primitive)


