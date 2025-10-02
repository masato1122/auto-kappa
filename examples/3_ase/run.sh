#
# This example generates VASP scripts with Auto_kappa. 
# INCAR, POSCAR, POTCAR, and KPOINTS are generated in ./out directory, 
# if the script works properly.
# Auto_kappa must be added to PYTHONPATH before running this example.
#
tdir=./tools

if [ ! -e $VASP_PP_PATH/potpaw_PBE ]; then
    echo " VASP_PP_PATH variable may not set properly."
    echo " Please visit the GitHub page: https://github.com/masato1122/auto-kappa"
    exit
fi

python $tdir/poscar2relax.py \
    -f ./POSCAR-unitcell


