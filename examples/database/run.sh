
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:${CONDA_PREFIX}/lib

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

akrun \
    --directory $dir_db \
    --material_name $mpid \
    --ncores $ncores \
    --nmax_suggest 200 \
    --frac_nrandom 0.02 \
    --verbosity 1

