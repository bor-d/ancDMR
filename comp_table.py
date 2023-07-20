# -*- coding: utf-8 -*-



import pandas as pd
import numpy as np
import os
import sys

list_dir = sys.argv[1]
cpg_file = sys.argv[2]

result_dt = pd.DataFrame(pd.read_table(cpg_file, sep='\t'), columns = ["pos"])
result_dt.set_index('pos')


def reading(directory, filename):
	temp = pd.DataFrame(pd.read_table(directory+filename, sep='\t', names=["chr", "pos", "pos1", "meth"]), columns=["chr", "pos", "pos1", "meth"])
	# add column names as standard bed files don't use header
	
            del temp["pos1"]
            del temp["chr"]
            print("Adding:" + str(filename) # making sure everything works
            
            # here I adjusted script for bed files of different format
            if filename.startswith("mod"): 
                temp = dict(zip((temp.pos+1), temp.meth))
                named_org = filename.rpartition('.')[0]
            elif filename.startswith("osteo"):
                temp = dict(zip((temp.pos+1), (temp.meth/100)))
                named_org = filename.rpartition('.')[0]
            elif filename.startswith("bone"):
                temp = dict(zip(temp.pos, (temp.meth/100)))
                named_org = filename.rpartition('.')[0]
                
            # adding new column to final table    
            result_dt[named_org] = result_dt['pos'].map(temp)
        else:
            return 0

[reading(list_dir, filename) for filename in os.listdir(list_dir)]

result_dt.to_csv(list_dir + 'modern_compare.txt', sep='\t', float_format='%1.7f')

