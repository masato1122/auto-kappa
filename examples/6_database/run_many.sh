#!/bin/sh
WORKDIR=$(cd $(dirname $0); pwd)

logdir=./logs
if [ ! -e $logdir ]; then
    mkdir $logdir
fi

#############################################
## parameters depending on the environment
dir_phonondb=/home/apdb/phonondb-20180417
nmax=10          ## maximum number of the loop
nprocs=32        ## number of cores to be used for VASP and Alamode
mpirun="mpirun"
#
DIR_AUTOKAPPA=$HOME/src/auto-kappa
##############################################

###########
imin=140       ## Initial material ID 
imax=160       ## Last material ID
neach=5        ## number of materials in a loop (average : 2-3 node*days/mater)
###########

for i in `seq 0 $nmax`; do

    anphon_para=mpi
    
    count=0
    i0=$imin
    for ii in `seq $imin $imax`; do
        mpid="mp-${ii}"
        dir_db=$dir_phonondb/${mpid}
        i1=$ii
        if [ ! -e $dir_db ]; then
            continue
        fi

        flag=`python $DIR_AUTOKAPPA/tools/check_apdb_log.py -d $mpid`
        if [ $flag != "NotYet" ]; then
            continue
        fi
        echo $mpid
        
        count=`expr $count + 1`
        if [ $count == "$neach" -o $ii == "$imax" ]; then
            break
        fi
    done
    i1=$ii
    imin=`expr $i1 + 1`
    echo $i : $count, ${i0} - ${i1}
    label=${i}_${i0}-${i1}
    
ofile=a.sh
##################################################################
##
## Please modify depending on the environment, job scheduler, etc.
##
##################################################################
cat >$ofile <<EOF
#!/bin/sh
#PBS -q default
#PBS -l nodes=1:ppn=${nprocs}
#PBS -j oe
#PBS -N $label

export LANG=C
export OMP_NUM_THREADS=1
cd \$PBS_O_WORKDIR

rm ${label}.o*

for id in \`seq $i0 $i1\`; do

    mpid="mp-\${id}"
    
    dir_db=$dir_phonondb/\${mpid}
    if [ ! -e \$dir_db ]; then
        continue
    fi

    ### start a calculation
    akrun \\
        --directory \$dir_db \\
        --outdir \$mpid \\
        --nprocs $nprocs \\
        --mpirun ${mpirun} \\
        --anphon_para $anphon_para \\
        > $logdir/log_\${mpid}.txt
    
done
EOF

################################################################
#### Remove the comment-out before "qsub" to submit the job 
#### and, if neccesarry, modify "qsub" command.
qsub $ofile
################################################################

    if [ $i1 == $imax ]; then
        break
    fi
done

