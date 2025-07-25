# -*- coding: utf-8 -*-
import numpy as np
from optparse import OptionParser

import os
import os.path
import glob

def get_all_directories(root_directory):
    directory_names = []
    for root, dirs, files in os.walk(root_directory):
        for dir_name in dirs:
            directory_names.append(os.path.join(root, dir_name))
    return directory_names

def _check_file_for_string(filename, target_string):
    """ Check whther ``filename`` contains a target string ``target_string`` or
    not.
    """
    with open(filename, 'r') as f:
        content = f.read()
        if target_string in content:
            return True
        else:
            return False

def _contains_vasp_result(dir_name):
    """ Check if the target directory contains VASP result or not.
    """
    fns = glob.glob(dir_name + "/vasprun.xml")
    if len(fns) == 0:
        return False
    else:
        return True

def _is_text_file(filename):
    """ Check if the target file is a text file or not.
    """
    try:
        with open(filename, 'r', encoding='utf-8') as file:
            file.read()
        return True
    except UnicodeDecodeError:
        return False

def _contains_alamode_result(dir_name):
    """ Check if the target directory contains ALAMODE result or not.
    """
    fns = glob.glob(dir_name + "/*")
    
    for fn in fns:
        ##
        if os.path.isfile(fn) == False:
            continue
        if _is_text_file(fn) == False:
            continue
        ##
        if _check_file_for_string(fn, "Program ALM"):
            return True
        elif _check_file_for_string(fn, "Program ANPHON"):
            return True
    
    return False

def _get_data_type(dir_name):
    """ Get data type in a target directory
    """
    if _contains_vasp_result(dir_name):
        data_type = "vasp"
    elif _contains_alamode_result(dir_name):
        data_type = "alamode"
    else:
        data_type = "others"
    return data_type

def _get_kind_of_time(relative_path):
    
    def _get_kind_of_time_eachstructure(data):
        kind = None
        d1 = data[1]
        if d1 == "relax":
            kind = "relax"
            
        elif d1 in ["nac", "harm"]:
            kind = "harmonic"
        
        elif d1 == "cube":
            if len(data) >= 3:
                d2 = data[2]
                if "force" in d2 or "suggest" in d2:
                    kind = "force(cube)"
                elif "kappa" in d2:
                    tmp = d2.split("_")[-1]
                    try:
                        kpts = [int(dd) for dd in tmp.split("x")]
                        kind = "kappa(%s)" % tmp
                    except Exception:
                        kind = "kappa"
        return kind
    
    data = relative_path.split("/")
    if "sc-" not in data[1]:
        kind = _get_kind_of_time_eachstructure(data)
    else:
        data2 = []
        for i in range(1, len(data)):
            data2.append(data[i])
        ##
        label_sc_mat = data[1].split("-")[1]
        try:
            kind = _get_kind_of_time_eachstructure(data2)
            kind += "(SC:%s)" % label_sc_mat
        except Exception:
            #print(" Error with", data2)
            kind = "others"
        
    if kind is None:
        kind = "others"

    return kind

#def _get_time_outcar_in_tar(tar_file_path, file_to_read):
#    """ Read and return the simulation time for VASP calculation compressed as a
#    tar.gz file.
#
#    Args
#    =====
#    tar_file_path : string
#        .tar.gz file name
#    
#    file_to_read : string
#        relative path in .tar.gz directory, which may be "OUTCAR".
#    """
#    import tarfile
#    
#    line_searched = "Total CPU time used"
#    
#    with tarfile.open(tar_file_path, 'r:gz') as tar:
#        try:
#            file_info = tar.getmember(file_to_read)
#            with tar.extractfile(file_info) as f:
#                content = f.read().decode('utf-8')
#                lines = content.split("\n")
#                print(lines)
#                for il, line in enumerate(lines, start=1):
#                    if line_searched in line:
#                        print(line)
#        
#        except KeyError:
#            print(f"File '{file_to_read}' not found in the archive.")
#    
#    print(tar_file_path)
#    exit()

def _get_time_each(dir_name, dtype):
    """ Get simulation time for each directory

    Args
    =====
    dir_name : string
        directory name
    dtype : string
        "vasp" or "alamode"
    """
    from auto_kappa.alamode.log_parser import (
            _get_alamode_runtime, _extract_data)
    
    if dtype == "alamode":
        duration = 0.
        fns = glob.glob(dir_name + "/*")
        for fn in fns:
            if os.path.isfile(fn):
                if _is_text_file(fn):
                    out = _get_alamode_runtime(fn)
                    if out is not None:
                        duration += out['value']    ## sec
    
    else:
        duration = 0.

        ### check OUTCAR
        outcar = dir_name + "/OUTCAR"
        time_single = 0.
        if os.path.exists(outcar):
            value = _extract_data(outcar, "Total CPU time used")
            if value is not None:
                time_single = float(value[0])
        
        ### check error.*.tar.gz
        dir_errors = glob.glob(dir_name + "/error.*.tar.gz")
        
        ### estimated time
        duration = time_single * (len(dir_errors) + 1)
    
    return duration

def get_times(base_dir):

    times = {}

    all_dirs = get_all_directories(base_dir)
    
    for dir_name in all_dirs:

        dtype = _get_data_type(dir_name)
        
        if dtype == "others":
            continue
        
        ##
        kind = _get_kind_of_time(dir_name.replace(base_dir, "."))
        
        duration = _get_time_each(dir_name, dtype)
        
        if duration is None:
            continue

        if kind not in times.keys():
            times[kind] = duration
        else:
            times[kind] += duration
    
    ###
    #total_time = 0.
    labels = list(times.keys())
    durations = []
    for lab in labels:
        durations.append(times[lab])
        #total_time += times[lab]

    #durations.append(total_time)
    #labels.append("total")
    
    return durations, labels

#def plot_times(directory, figname="fig_times.png"):
#    
#    times, labels = get_times(directory)
#
#    from auto_kappa.plot.pltalm import plot_times_with_pie
#    plot_times_with_pie(times, labels, figname=figname)
#
#def main(options):
#    
#    plot_times(options.directory, figname=options.figname)
#
#if __name__ == '__main__':
#    parser = OptionParser()
#    
#    parser.add_option("-d", "--directory", dest="directory", type="string",
#            default=None, help="directory name")
#    
#    parser.add_option("-f", "--figname", dest="figname", type="string",
#            default=None, help="figname name")
#    
#    (options, args) = parser.parse_args()
#    
#    main(options)

