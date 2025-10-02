
### Description
# volume_relaxation : relaxation calc. with the equation of state
# max_natoms        : maximum number of atoms in a supercell
# nmax_suggest      : maximum number of displacement patterns for the FD method

### Abbreviation
# FC : force constant
# FD : finite displacement

nprocs=2

### Input for each material
file_structure=./POSCAR.Si
mpid=mp-149   ## material ID for Silicon

### Commands
c_mpirun=mpirun
c_vasp=vasp
c_alm=alm
c_anphon=anphon

### Variables
vol_relax=1     ## Structural relaxation with the equation of state
nmax_suggest=1  ## Maximum number of displacement patterns for the FD method. If 1, always use LASSO approach.
max_natoms=150  ## Maximum number of atoms in a supercell

### Check command existence
for command in $c_mpirun $c_vasp $c_alm $c_anphon; do
    if ! command -v $command &> /dev/null; then
        echo "$command could not be found"
        exit
    fi
done

### Run auto-kappa
akrun \
    --file_structure $file_structure \
    --outdir $mpid \
    --nprocs $nprocs \
    --mpirun         $c_mpirun \
    --command_vasp   $c_vasp \
    --command_alm    $c_alm \
    --command_anphon $c_anphon \
    \
    --volume_relaxation $vol_relax \
    --max_natoms   $max_natoms \
    --nmax_suggest $nmax_suggest \
    --calculate_forces 1
