#
# Test job with "--analyze_with_larger_supercell" option.
#
# --analyze_with_larger_supercell : int
# If this option has nonzero value, phonon properties are analyzed with larger supercells 
# when negative frequencies exist.
#
# --max_loop_for_largesc : int
#
# --delta_max_natoms : int, default=50
#
#
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:${CONDA_PREFIX}/lib

ncores=2

poscar=POSCAR.init   ## file name of the structure
material_name=Si     ## this option is used only for the name of output directory

akrun \
    --file_structure $poscar \
    --material_name $material_name \
    --ncores $ncores \
    --mpirun mpirun \
    --command_vasp vasp \
    --command_alm alm \
    --command_anphon anphon \
    --volume_relaxation 1 \
    --analyze_with_larger_supercell 1 \
    --max_loop_for_largersc 2

