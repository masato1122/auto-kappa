#!/usr/bin/env python
import argparse
from custodian.custodian import Custodian
from custodian.vasp.handlers import (VaspErrorHandler, 
                                     UnconvergedErrorHandler)
from custodian.vasp.jobs import VaspJob

def main(options):
    
    handlers = [
            VaspErrorHandler(), 
            #UnconvergedErrorHandler(),    
            ]
    
    cmd = f"{options.mpirun} -n 2 {options.command_vasp}"
    print(f"Running command: {cmd}")
    
    jobs = VaspJob.double_relaxation_run(cmd.split())
    c = Custodian(handlers, jobs, max_errors=10)
    c.run()

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Input parameters')
    parser.add_argument('--mpirun', dest='mpirun', type=str,
                        default="mpirun", help="MPIRUN command (mpirun, mpiexec.hydra, etc.) [mpirun]") 
    parser.add_argument('--command_vasp', dest='command_vasp', type=str,
                        default="vasp", help="VASP command [vasp]")
    args = parser.parse_args()
    main(args)
