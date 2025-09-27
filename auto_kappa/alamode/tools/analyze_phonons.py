#!/usr/bin/env python
#
# analyze_phonons.py
#
# Simple interface to the command line script "analyze_phonons.cpp".
#
# Copyright (c) 2014, 2015, 2016 Terumasa Tadano
# Modified by Masato Ohnishi to adjust auto-kappa
#
# This file is distributed under the terms of the MIT license.
# Please see the file 'LICENCE.txt' in the root directory
# or http://opensource.org/licenses/mit-license.php for information.
#
# import sys
import os
import os.path
import warnings
# import subprocess

def analyze_obj():
    
    if os.environ['ALAMODE_PATH'] is None:
        print("")
        warnings.warn(" Error: ALAMODE_PATH must be set.")
        print(" ${ALAMODE_PATH}/tools/analyze_phonons must exist.")
        print("")
        exit()
    
    filename = os.environ['ALAMODE_PATH'] + "/tools/analyze_phonons"
    if os.path.exists(filename) == False:
        print("")
        warnings.warn(" Error: %s must exist." % filename)
        print("")
        exit()
        
    obj = filename + " "
    return obj

# def write_lifetime_at_given_temperature(file_result: None, 
#         file_isotope=None, outfile=None,
#         temperature=300, kpoint=None, mode=None, average_gamma=True):
    
#     ### analyze tau
#     calc = "tau"

#     if outfile is None:
#         outfile = "tau_%dK.dat" % int(temperature)
    
#     ### isotope
#     isotope = "0"
#     if file_isotope is not None:
#         if os.path.exists(file_isotope):
#             isotope = "1"
#         else:
#             print("")
#             warnings.warn(" WARNNING: %s does not exist." % file_isotope)
#             print("")
#     if isotope == "0":
#         file_isotope = "none"
    
#     if average_gamma:
#         avg = "1"
#     else:
#         avg = "0"
    
#     if kpoint is None:
#         beg_k = 1
#         end_k = 0
#     else:
#         if len(kpoint.split(':')) == 1:
#             beg_k = int(kpoint)
#             end_k = beg_k
#         elif len(kpoint.split(':')) == 2:
#             arr = kpoint.split(':')
#             beg_k, end_k = int(arr[0]), int(arr[1])
#         else:
#             sys.exit("Invalid usage of kpoint")
    
#     if mode is None:
#         beg_s = 1
#         end_s = 0
#     else:
#         if len(mode.split(':')) == 1:
#             beg_s = int(options.mode)
#             end_s = beg_s
#         elif len(mode.split(':')) == 2:
#             arr = mode.split(':')
#             beg_s, end_s = int(arr[0]), int(arr[1])
#         else:
#             sys.exit("Invalid usage of --mode for --calc=tau")
    
#     command = analyze_obj() + file_result + " "\
#             + str(calc) + " " + str(avg) + " "\
#             + str(beg_k) + " " + str(end_k) + " "\
#             + str(beg_s) + " " + str(end_s) + " " + str(temperature)\
#             + " " + isotope + " " + file_isotope + " > " + outfile
    
#     subprocess.call(command, shell=True)

