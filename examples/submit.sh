#!/bin/sh
#PBS -q default
#PBS -l nodes=1:ppn=24
#PBS -j oe
#PBS -N auto_kappa

export LANG=C
export OMP_NUM_THREADS=1
nprocs=24
cd $PBS_O_WORKDIR
#rm auto_kappa.o*

export PYTHONPATH=$PYTHONPATH:..

mpid=mp-149    ## Si

dir_db=./phonondb-20180417/${mpid}
python ../scripts/auto_kappa.py \
    --directory $dir_db \
    --mpid $mpid \
    --nprocs $nprocs

