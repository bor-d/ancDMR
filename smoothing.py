# -*- coding: utf-8 -*-
"""
Created on Sun May 29 21:13:45 2022

@author: Shinji
"""

import pandas as pd
import numpy as np
#import sys

#directory with files for exact chrom
list_dir = sys.argv[1]

#ref file is generated by perl script cgpositionFinder.pl and has a name ref/chr*.CpG.positions.txt
ref_file = sys.argv[2]

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

result_dt = pd.DataFrame(pd.read_csv(ref_file, names = ["chr", "pos"], sep = '\t')['pos'])
#result_dt.set_index('pos')
sample_sheet = pd.DataFrame(pd.read_csv(list_dir+"sample_file.txt", sep='\t', 
                                          names=["name", "org", "cov"]))
#result_dt.set_index('name')

def reading_smoothing(directory, sample):
    filename = str(directory)+str(sample)+'.txt'
    
    #deciding sliding window based on metadata
    coverage = sample_sheet.loc[sample_sheet["name"] == sample, 'cov'].item()
    if coverage > 29:
        window = 25
    elif coverage > 10:
        window = 35
    else:
        window = 50
        
    
    try:
        if (sample == "den") or (sample == "nea") or (sample == 'vi33'):
            temp = pd.DataFrame(pd.read_csv(filename, sep=' ', names=["pos", "dm", "gkh"]), columns=["pos", "dm", "gkh"])
            del temp["gkh"]
            temp.dm = temp.dm/100
            print(sample)
            print(temp[0:3])
            temp["smooth"] = running_mean(pd.to_numeric(temp["dm"]), window)
            del temp["dm"]
            temp = dict(zip(temp.pos-1, temp.smooth))
        elif sample.endswith("modern") or sample.startswith('lean') or sample.startswith('obese'):
            temp = pd.DataFrame(pd.read_csv(filename, sep='\t', names=["chr", "pos", "pos1", "meth"]), columns=["chr", "pos", "pos1", "meth"])
            del temp["chr"], temp["pos1"]
            if sample.startswith("osteo"):
                temp = dict(zip(temp.pos, temp.meth))
            elif sample.startswith('lean') or sample.startswith('obese'):
                temp = dict(zip(temp.pos, temp.meth/100))
            else:
                temp.pos = temp.pos-1
                temp = dict(zip(temp.pos, temp.meth/100))
        else:
            temp = pd.DataFrame(pd.read_csv(filename, sep=' ', names=["pos", "dm"]), columns=["pos", "dm"])
            temp["pos"] =  pd.to_numeric(temp["pos"])
            print(sample)
            print(temp[0:3])
            temp["smooth"] = running_mean(temp["dm"], window)
            del temp["dm"]
            temp = dict(zip(temp.pos, temp.smooth))
        result_dt[sample] = result_dt['pos'].map(temp)
        print(sample + " smoothed")
    except FileNotFoundError:
        print(sample+" not found in folder")
    

[reading_smoothing(list_dir, sample) for sample in sample_sheet['name']]

np.savetxt('shit.txt', result_dt,
           fmt = "%d\t%1.7f\t%1.7f\t%1.7f\t%1.7f\t%1.7f\t%1.7f\t%1.7f\t%1.7f\t%1.7f\t%1.7f\t%1.7f\t%1.7f\t%1.7f\t%1.7f\t%1.7f\t%1.7f\t%1.7f\t%1.7f\t%1.7f\t%1.7f\t%1.7f\t%1.7f\t%1.7f\t%1.7f\t%1.7f\t%1.7f\t%1.7f\t%1.7f\t%1.7f\t%1.7f\t%1.7f\t%1.7f\t%1.7f",
           delimiter = "\t")
#format is optional, delimeter '\t' is a must for flawless work of other scripts
