# -*- coding: utf-8 -*-
import os.path
import numpy as np
from optparse import OptionParser

def main(options):
    
    filename = 'finished.txt'
    if os.path.exists(filename) == False:
        print("0")
        exit()
    
    lines = open(filename, 'r').readlines()
    for line in lines:
        if options.mpid in line:
            print("1")
            exit()
    
    print("0")

    
if __name__ == '__main__':
    parser = OptionParser()
    parser.add_option("-f", "--mpid", dest="mpid", type="string",
            help="input file name")
    (options, args) = parser.parse_args()
    #file_check(options.filename)
    main(options)
