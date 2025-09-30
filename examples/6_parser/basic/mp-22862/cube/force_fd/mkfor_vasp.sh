#!/bin/sh
WORKDIR=$(cd $(dirname $0); pwd)
OFILE=a.sh

for i in `seq 1 50`; do
    
    cd $WORKDIR
    if [ ! -e $i ]; then
        continue
    fi

cd $i
cat >$OFILE<<EOF
#!/bin/sh
#PBS -q workq
#PBS -l nodes=1:ppn=32
#PBS -j oe
#PBS -N job${i}

export LANG=C
export OMP_NUM_THREADS=1
nprocs=1
cd \$PBS_O_WORKDIR
rm job${i}.o*

mpirun -n 32 vasp 
EOF
qsub $OFILE; sleep 0.1
done

