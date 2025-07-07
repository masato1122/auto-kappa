#!/bin/sh
WORKDIR=$(cd $(dirname $0); pwd)

logdir=./logs
if [ ! -e $logdir ]; then
    mkdir $logdir
fi

#############################################
## parameters depending on the environment
dir_phonondb=$HOME/phonondb-20180417
nmax=10          ## maximum number of the loop
ncores=32        ## number of cores to be used for VASP and Alamode
mpirun="mpirun"
#
#DIR_AUTOKAPPA=$HOME/src/auto-kappa
DIR_AUTOKAPPA=/home/ohnishi/work/PhDB/auto-kappa
##############################################

###########
imin=145       ## Initial material ID 
imax=155       ## Last material ID
neach=5        ## number of materials in a loop (average : 2-3 node*days/mater)
###########

##############################
dir_pre_base=../6_database
##############################

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
        
        ### check the previous calculation
        dir_pre=$dir_pre_base/$mpid
        flag=`python $DIR_AUTOKAPPA/tools/check_apdb_log.py -d $dir_pre`
        if [ $flag != "Finished_harmonic" ]; then
            continue
        fi
        echo $mpid $flag
        
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
#PBS -l nodes=1:ppn=${ncores}
#PBS -j oe
#PBS -N $label

export LANG=C
export OMP_NUM_THREADS=1
cd \$PBS_O_WORKDIR

rm ${label}.o*

dir_phonondb=$dir_phonondb
DIR_AUTOKAPPA=$DIR_AUTOKAPPA
dir_pre_base=$dir_pre_base

for id in \`seq $i0 $i1\`; do

    mpid="mp-\${id}"
    
    ### check Phonondb
    dir_db=\$dir_phonondb/\${mpid}
    if [ ! -e \$dir_db ]; then
        continue
    fi
    
    ### check the previous calculation
    dir_pre=\$dir_pre_base/\${mpid}
    flag=\`python \$DIR_AUTOKAPPA/tools/check_apdb_log.py -d \$dir_pre\`
    if [ \$flag != "Finished_harmonic" ]; then
        continue
    fi
    
    ### start a calculation
    akrun \\
        --directory \$dir_db \\
        --material_name \$mpid \\
        --ncores $ncores \\
        --mpirun ${mpirun} \\
        --anphon_para $anphon_para \\
        --volume_relaxation 1 \\
        > $logdir/log_\${mpid}.txt
    
done
EOF

################################################################
#### Remove the comment-out before "qsub" to submit the job 
#### and, if neccesarry, modify "qsub" command.
#qsub $ofile
################################################################

    if [ $i1 == $imax ]; then
        break
    fi
done

