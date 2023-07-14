# -*- coding: utf-8 -*-
"""
Created on Sun May 29 21:13:45 2022

@author: Shinji
"""

import pandas as pd
import numpy as np
#from numpy.lib.stride_tricks import sliding_window_view
import os
#import codecs
#import re
import sys

list_dir = sys.argv[1]
cpg_file = sys.argv[2]

def running_mean(x, N):
    out = np.zeros_like(x, dtype=np.float64)
    dim_len = x.shape[0]
    for i in range(dim_len):
        if N%2 == 0:
            a, b = i - (N-1)//2, i + (N-1)//2 + 2
        else:
            a, b = i - (N-1)//2, i + (N-1)//2 + 1

        #cap indices to min and max indices
        a = max(0, a)
        b = min(dim_len, b)
        out[i] = np.mean(x[a:b])
    return out

result_dt = pd.DataFrame(pd.read_table(cpg_file, sep='\t', names=["pos"]))
result_dt.set_index('pos')


def reading_smoothing(directory, filename):
    if filename.endswith(".txt"):
            temp = pd.DataFrame(pd.read_table(directory+filename, sep=' ', names=["pos", "dm"]), columns=["pos", "dm"])
            temp["pos"] = temp["pos"] + 1
            print(filename)
            print(temp[0:3])
        temp["smooth"] = running_mean(temp["dm"], 25)
        del temp["dm"]
        named_org = filename.rpartition('_')[0]
        temp = dict(zip(temp.pos, temp.smooth))
        print(filename + " smoothed")
        result_dt[named_org] = result_dt['pos'].map(temp)
    elif filename.endswith(".bed"):
        temp = pd.DataFrame(pd.read_table(directory+filename, sep='\t', names=["chr", "pos", "pos1", "meth"]), 
                            columns=["chr", "pos", "pos1", "meth"])
        del temp["pos1"]
        del temp["chr"]
        print(filename)
        print(temp[0:3])
        named_org = filename.rpartition('_')[0]
        result_dt[named_org] = result_dt['pos'].map(temp)
    else:
        return 0
    #return result_dt

[reading_smoothing(list_dir, filename) for filename in os.listdir(list_dir)]

columnNames = df.columns.values.tolist()
np.savetxt(list_dir+'smoothed.txt', result_dt,
           fmt = "%d\t%1.7f\t%1.7f\t%1.7f\t%1.7f\t%1.7f\t%1.7f\t%1.7f\t%1.7f\t%1.7f\t%1.7f\t%1.7f\t%1.7f\t%1.7f",
           delimiter = "\t", header='\t'.join(columnNames))
