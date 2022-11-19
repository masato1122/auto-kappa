from custodian.custodian import Custodian
from custodian.vasp.handlers import (
        VaspErrorHandler, 
        UnconvergedErrorHandler, 
        )
from custodian.vasp.jobs import VaspJob

handlers = [
        VaspErrorHandler(), 
        #UnconvergedErrorHandler(),    
        ]

#cmd = "mpiexec.hydra -n 2 vasp"
cmd = "mpirun -n 2 vasp"
jobs = VaspJob.double_relaxation_run(cmd.split())

c = Custodian(handlers, jobs, max_errors=10)
c.run()

