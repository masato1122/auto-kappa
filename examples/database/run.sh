
########################################
##
## Parameters you may need to modify
##
export PYTHONPATH=$PYTHONPATH:../..

mpid=mp-149    ## Si

ncores=2
dir_db=./phonondb-20180417/${mpid}
#########################################

if [ ! -e $dir_db ]; then
    continue
fi

../../scripts/auto_kappa.py \
    --directory $dir_db \
    --mpid $mpid \
    --ncores $ncores \
    --nmax_suggest 200 \
    --frac_nrandom 0.02 \
    --verbosity 1

