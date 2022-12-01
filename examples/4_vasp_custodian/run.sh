# 
# One can confirm if VASP and Custodian works property with this example. 
# If VASP starts and different files such as OUTCAR, OSZICAR, etc. are created, 
# they are working properly, and one can stop the job.
# The following line in vasp_custodian.py may need to be modified depending on the user's
# environment.
#
# cmd = "mpirun -n 2 vasp"
#
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:${CONDA_PREFIX}/lib
python vasp_custodian.py

