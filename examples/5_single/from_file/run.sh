

nprocs=2

file_structure=./POSCAR.Si
mpid=mp-149   ## material ID for Silicon

akrun_develop \
    --file_structure $file_structure \
    --outdir $mpid \
    --nprocs $nprocs \
    --mpirun mpirun \
    --command_vasp vasp \
    --command_alm alm \
    --command_anphon anphon \
    --calculate_forces 0


