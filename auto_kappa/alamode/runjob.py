# -*- coding: utf-8 -*-
#
# runjob.py
#
# This script helps to run ALAMODE job
#
# Copyright (c) 2024 Masato Ohnishi
#
# This file is distributed under the terms of the MIT license.
# Please see the file 'LICENCE.txt' in the root directory
# or http://opensource.org/licenses/mit-license.php for information.
#
import os
import os.path
import logging
import subprocess
import time

try:
    import psutil
except ImportError:
    pass

from auto_kappa.alamode.io import wasfinished_alamode
from auto_kappa.alamode.errors import check_unexpected_errors

logger = logging.getLogger(__name__)

def run_alamode(
        filename, logfile, workdir='.', neglect_log=0, file_err="std_err.txt",
        mpirun='mpirun', nprocs=1, nthreads=1, command='anphon',
        max_num_corrections=None):
    """ Run alamode with a command (alm or anphon)

    Args
    ======
    filename : string
        input script of Alamode in workdir

    logfile : string
        log file name in workdir

    workdir : string
        work directory

    Return
    =======
    int :
        ``-1`` when the job had been finished.
        ``0`` when the job was conducted.
        ``1`` when the job was not finished.

    """
    omp_keys = ["OMP_NUM_THREADS", "SLURM_CPUS_PER_TASK"]
    
    ### change directory
    dir_init = os.getcwd()
    os.chdir(workdir)

    ## If the job has been finished, the same calculation is not conducted.
    ## The job status is determined from *.log file.
    status = None
    count = 0
    while True:
        
        if wasfinished_alamode(logfile) and neglect_log == False:
            status = -1
            break
        
        if count == 0:
            ppn_i = nprocs
            nth_i = nthreads
        else:
            ppn_prev = ppn_i
            nth_prev = nth_i
            
            ppn_i = max(1, int(ppn_prev / 2))
            if ppn_i > 1:
                ppn_i += int(ppn_i % 2)
            
            if ppn_prev == ppn_i:
                status = 1
                break
            
            ###
            msg = "\n Processes per node : %d => %d" % (ppn_prev, ppn_i)
            logger.info(msg)
        
        ### set number of parallelization
        cmd = "%s -n %d %s %s" %(mpirun, ppn_i, command, filename)
        
        ### set OpenMP
        for key in omp_keys:
            os.environ[key] = str(nth_i)
        
        ### run the job
        _run_job(cmd, logfile=logfile, file_err=file_err)
        
        count += 1
        
        if wasfinished_alamode(logfile):
            status = 0
            break
        
        if max_num_corrections is not None:
            if count >= max_num_corrections:
                status = 1
                break
    
    ### set back OpenMP
    for key in omp_keys:
        os.environ[key] = "1"
    
    ### check logfile
    dir_base = dir_init + "/" + workdir.replace(dir_init, ".").split("/")[1]
    check_unexpected_errors(logfile, dir_base=dir_base)
    
    #### Return to the original directory
    os.chdir(dir_init)
    
    return status

def _run_job(cmd, logfile="log.txt", file_err="std_err.txt"):
    """ Run a job with subprocss
    
    Args
    -------

    cmd : string
        command to run a job
    
    """    
    ## run the job!!
    status = None
    with open(logfile, 'w') as f, open(file_err, "w", buffering=1) as f_err:
        
        ### Error file, termination check
        proc = subprocess.Popen(
                cmd, shell=True, env=os.environ,
                stdout=f, stderr=f_err)
        
        count = 0
        mem_max = -1
        while True:
            
            if proc.poll() is not None:
                break
            
            ### get memory info if available
            mem_percentage = 0.
            try:
                mem_info = psutil.virtual_memory()
                mem_max = max(mem_max, mem_info.used)
                mem_tot = mem_info.total
                mem_percentage = mem_info.percentage

                if mem_percentage > 95.:
                    logger.info("\n Warning: memory usage is %.2f%%" % (
                        mem_info.percentage))
                    break

            except Exception:
                pass

            waiting_time = min(10, count)
            time.sleep(waiting_time)
            count += 1
        
        if mem_max > 0.:
            msg = "\n Maximum memory usage : %.3f GB" % (mem_max / 1e9)
            logger.info(msg)

        status = proc.poll()
    
    return 0

