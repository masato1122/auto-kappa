from aiida_alamode.io import AlmInput

filename = 'POSCAR-supercell'
dfset = 'dfset.txt'

## mode "suggest"
alminp = AlmInput.from_structure_file(
        filename,
        mode='suggest',
        norder=2
        )
alminp.to_file()
#alminp.to_file(filename='filename.in')

## mode "optimize"
## Note: Cutoff should be given for high-order IFCs.
alminp = AlmInput.from_structure_file(
        filename,
        mode='optimize',
        norder=2,
        cutoff={'Si-Si': [None, None]},
        dfset=dfset
        )
alminp.to_file()

