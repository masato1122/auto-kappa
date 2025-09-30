
### Description
# volume_relaxation : relaxation calc. with the equation of state
# max_natoms        : maximum number of atoms in a supercell
# nmax_suggest      : maximum number of displacement patterns for the FD method

### Abbreviation
# FC : force constant
# FD : finite displacement


nprocs=2

file_structure=./POSCAR.Si
mpid=mp-149   ## material ID for Silicon

akrun \
    --file_structure $file_structure \
    --outdir $mpid \
    --nprocs $nprocs \
    --mpirun mpirun \
    --command_vasp vasp \
    --command_alm alm \
    --command_anphon anphon \
    --calculate_forces 1 \
    --volume_relaxation 1 \
    --max_natoms 150 \
    --nmax_suggest 1

