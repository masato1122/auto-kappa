#!/bin/sh
WORKDIR=$(cd $(dirname $0); pwd)

logdir=./logs
if [ ! -e $logdir ]; then
    mkdir $logdir
fi

##########################################
## parameters depending on the environment
dir_phonondb=/home/apdb/phonondb-20180417
imax=50000
nmax=10
ncores=32
mpirun="mpirun"
###########################################

###########
imin=40000   ## need to be modified each series of job submissions
neach=5      ## number of materials in a loop (average : 2-3 node*days/mater)
###########

for i in `seq 0 $nmax`; do
    
    count=0
    i0=$imin
    for ii in `seq $imin $imax`; do
        mpid="mp-${ii}"
        dir_db=$dir_phonondb/${mpid}
        i1=$ii
        if [ ! -e $dir_db ]; then
            continue
        fi
        
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
##################################################
##                                              ##
## Please modify contents in the created file.  ##
##                                              ##
##################################################
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

### environment setting
#LD_LIBRARY_PATH=\$LD_LIBRARY_PATH:/home/ohnishi/.anaconda3/lib
#export LD_LIBRARY_PATH

for id in \`seq $i0 $i1\`; do

    mpid="mp-\${id}"
    
    dir_db=$dir_phonondb/\${mpid}
    if [ ! -e \$dir_db ]; then
        continue
    fi

    ### start a calculation
    akrun \\
        --directory \$dir_db \\
        --material_name \$mpid \\
        --ncores $ncores \\
        --mpirun ${mpirun} \\
        --verbosity 1 \\
        > $logdir/log_\${mpid}.txt
    
    echo \$mpid >> finished.txt

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

