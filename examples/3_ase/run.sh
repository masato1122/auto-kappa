#
# This example generates VASP scripts with Auto_kappa. 
# INCAR, POSCAR, POTCAR, and KPOINTS are generated in ./out directory, 
# if the script works properly.
# Auto_kappa must be added to PYTHONPATH before running this example.
#
tdir=./tools

python $tdir/poscar2relax.py \
    -f ./POSCAR-unitcell


