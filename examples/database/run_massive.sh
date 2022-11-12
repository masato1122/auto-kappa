#!/bin/sh
#
# Generate and submit many jobs. Each job calculates "neach" materials.
#
WORKDIR=$(cd $(dirname $0); pwd)

##########################################
dir_phonondb=$HOME/phonondb-20180417   ## directory of Phonondb downloaded
imax=50000                             ## maximum ID to be calculated
##########################################

###########
imin=40000
neach=7
###########

nmax=10
for i in `seq 0 $nmax`; do
    
    count=0
    for ii in `seq $imin $imax`; do
        mpid="mp-${ii}"
        dir_db=$dir_phonondb/${mpid}
        i1=$ii
        if [ ! -e $dir_db ]; then
            continue
        fi
        
        flag=`python check_finished.py --mpid $mpid`
        if [ $flag == "1" ]; then
            continue
        fi
        count=`expr $count + 1`
        if [ $count == "$neach" -o $ii == "$imax" ]; then
            i0=$imin
            break
        fi
    done
    imin=`expr $i1 + 1`
    echo $i : $count, ${i0} - ${i1}
    
    
ofile=a.sh
cat >>$ofile <<EOF
#
# Modify this part to make a job script
#
EOF
#qsub $ofile

    if [ $i1 == $imax ]; then
        break
    fi
done

