#
# This example generates VASP scripts with auto_kappa. If INCAR, POSCAR, POTCAR, and KPOINTS
# are generated in ./out directory, they are working properly. 
# auto_kappa must be added to PYTHONPATH before running this example.
#
#
tdir=./tools

python $tdir/poscar2relax.py \
    -f ./POSCAR-unitcell


