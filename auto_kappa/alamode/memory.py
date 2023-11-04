# -*- coding: utf-8 -*-
#
# memory.py
#
# Copyright (c) 2023 Masato Ohnishi
#
# This file is distributed under the terms of the MIT license.
# Please see the file 'LICENCE.txt' in the root directory
# or http://opensource.org/licenses/mit-license.php for information.
#

import logging

logger = logging.getLogger(__name__)

def get_used_memory():
    """ Get and return the used memory """
    try:
        
        mem = _get_used_memory_psutil()
        return mem
    
    except ModuleNotFoundError:
        
        try:
            mem = _get_used_memory_resource()
            return mem

        except ModuleNotFoundError:
            
            return None

def _get_used_memory_psutil(busy_thread=80):

    import psutil
    mem_percent = psutil.virtual_memory().percent
    if mem_percent > busy_thred:
        msg = "\n Warning: memory usage percentage : %.1f%%" % mem_percent
        logger.warning(msg)
        
    mem = psutil.virtual_memory().used     ## byte
    return mem

def is_node_busy(busy_thred=80):
    """ Check if the node is busy or not. """
    try:
        import psutil
        cpu_percent = psutil.cpu_percent(interval=1)
    except Exception:
        cpu_percent = 0.
    return cpu_percent > busy_thred
    
def _get_used_memory_resource():

    import platform
    import resource
    system = platform.system()
    #r = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    #if system == "Darwin":
    #    mem = r / 1024 / 1024
    #else:
    #    mem = r / 1024
    #return mem
    return None

#print("psutil  :", _get_used_memory_psutil())
#print("resource:", _get_used_memory_resource())

