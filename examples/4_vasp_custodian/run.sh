#!/bin/sh
WORKDIR=$(cd $(dirname $0); pwd)
# 
# One can confirm if VASP and Custodian works property with this example. 
# If VASP starts and different files such as OUTCAR, OSZICAR, etc. are created, 
# they are working properly, and one can stop the job.
# The following line in vasp_custodian.py may need to be modified depending on the user's
# environment.
#
# cmd = "mpirun -n 2 vasp"
#
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:${CONDA_PREFIX}/lib

# Prepare output directory and VASP input files.
outdir=./out
if [ ! -e $outdir ]; then
    mkdir $outdir
fi
cp $WORKDIR/FILES/* $outdir

cd $outdir
#
# Run VASP with Custodian.
# Please modify the input parameters if necessary.
#
python $WORKDIR/tools/vasp_custodian.py \
    --mpirun mpirun \
    --command_vasp vasp
