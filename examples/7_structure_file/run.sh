#
# This script helps to run akrun command of auto_kappa.
# If this script works property, anharmonic phonon properties of Silicon 
# can be calculated automatically. Then, other materials can be easily 
# calculated if a few variables are modified in this script.
#
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:${CONDA_PREFIX}/lib

nprocs=2

poscar=POSCAR.init   ## file name of the structure
material_name=Si     ## name of the output directory

akrun \
    --file_structure $poscar \
    --material_name $material_name \
    --nprocs $nprocs \
    --mpirun mpirun \
    --command_vasp vasp \
    --command_alm alm \
    --command_anphon anphon 

