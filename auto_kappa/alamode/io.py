# -*- coding: utf-8 -*-
#from optparse import OptionParser
import numpy as np
import pandas as pd

def read_fcs(filename):
    
    ## il at the beginning
    lines = open(filename, 'r').readlines()
    ilines = {}
    order = None
    for il, line in enumerate(lines):
        data = line.split()
        if len(data) == 0:
            continue
        if "*FC" in line:
            if "**" in line:
                break
            order = int(data[0][-1].replace("FC", "")) - 1
            ilines[order] = []
            ilines[order].append(il + 1)
    
    ## il at the end
    for order in ilines:
        il0 = ilines[order][0]
        for il in range(il0, len(lines)):
            data = lines[il].split()
            if len(data) == 0:
                ilines[order].append(il-1)
                break
    
    ### Read all FCs
    dfs = {}
    for order in ilines:
        ndat = 6 + order
        il0 = ilines[order][0]
        il1 = ilines[order][1]
        all_data = []
        for j in range(ndat):
            all_data.append([])
        
        ## set columns: FC[Ry/a0^(order+1)], Distance[Bohr]
        columns = ['global_index', 'local_index', 'fc', 'multiplicity']
        types = ["int", "int", "float", "int"]
        for j in range(order+1):
            columns.append("iatom%d" % (j+1))
            types.append("str")
        columns.append("distance")
        types.append("float")
         
        ## read all data for FC*
        for il in range(il0, il1+1):
            data = lines[il].split()
            for j in range(ndat):
                if types[j] == "int":
                    all_data[j].append(int(data[j]))
                elif types[j] == "float":
                    all_data[j].append(float(data[j]))
                else:
                    all_data[j].append(data[j])
        
        ## convert to DataFrame
        df_each = pd.DataFrame()
        for j in range(ndat):
            df_each[columns[j]] = all_data[j]
        
        dfs[order] = df_each
    
    return dfs

