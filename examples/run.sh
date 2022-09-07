export PYTHONPATH=$PYTHONPATH:..

mpid=mp-149    ## Si

#dir_db=~/phonondb-20180417/${mpid}
dir_db=./phonondb-20180417/${mpid}
python auto_kappa.py \
    --directory $dir_db \
    --mpid $mpid \
    --nprocs 2

