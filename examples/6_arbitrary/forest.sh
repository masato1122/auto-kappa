#!/bin/sh
#PBS -q workq
#PBS -l nodes=1:ppn=32
#PBS -j oe
#PBS -N test-si

export LANG=C
export OMP_NUM_THREADS=1
cd $PBS_O_WORKDIR
rm test-si.o*

export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:${CONDA_PREFIX}/lib

ncores=32

poscar=POSCAR.init
material_name=Si

akrun \
    --file_structure $poscar \
    --material_name $material_name \
    --ncores $ncores \
    --mpirun mpirun \
    --command_vasp vasp \
    --command_alm alm \
    --command_anphon anphon \
    --verbosity 1 \
    > log/log-${material_name}.txt


