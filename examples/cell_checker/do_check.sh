
####################################################
##
## Modify this part depending on your environment
##
dir_autokappa=/home/ohnishi/work/PhDB/auto-kappa
dir0=~/phonondb-20180417
####################################################

tdir=$dir_autokappa/examples/cell_checker

##
## modify the range of ID
##
for id in `seq 0 10000`; do
    
    ##
    ## You may also need to modify this part.
    ##
    dir_each=$dir0/mp-${id}
    dir_apdb=./mp-${id}

    if [ ! -e $dir_each ]; then
        continue
    fi

    if [ ! -e $dir_apdb ]; then
        continue
    fi
    
    echo " === mp-${id} ==="
    python $tdir/cell_checker.py \
        --dir_phdb $dir_each \
        --dir_apdb $dir_apdb
done

