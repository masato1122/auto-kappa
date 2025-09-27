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

nprocs=2

poscar=POSCAR.Si   ## file name of the structure
mpid=mp-149

akrun \
    --file_structure $poscar \
    --outdir $mpid \
    --nprocs $nprocs \
    --analyze_with_larger_supercell 1 \
    --max_natoms 150 \
    --delta_max_natoms 50 \
    --max_loop_for_largersc 1

