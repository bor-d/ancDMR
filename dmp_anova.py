# -*- coding: utf-8 -*-
"""
Created on Wed Jun  1 20:44:41 2022

@author: Shinji
"""

import pandas as pd
#import statsmodels.api as sm
from scipy import stats as st
import numpy as np
#import scipy as sp
import sys

list_dir = sys.argv[1]
group_file = sys.argv[2]


datafile = list_dir+"smoothed.txt"
data = pd.read_csv(datafile, sep = '\t', header=0)
age_df = pd.read_csv(group_file, sep = '\t', header=0)

f = open(list_dir+"anova.txt", 'w')
f.write(str("pos"+"\t"+"anova_stat"+"\t"+"pval" +'\t'+'significant?'+"\n"))
f.close()

def anova_dmp(*args):
	temp_df = pd.DataFrame({"sample": data.columns[1:], "methylation": args[1:]})
	merged_df = pd.merge(temp_df, age_df, on="sample", how="left")
    
    F, p = st.f_oneway(pd.to_numeric(temp_df[temp_df["age"] == 'ANC']['value'], errors='coerce'), 
                       pd.to_numeric(temp_df[temp_df["age"] == 'HG']['value'], errors='coerce'), 
                       pd.to_numeric(temp_df[temp_df["age"] == 'MODERN']['value'], errors='coerce'))
    f = open(list_dir+"anova.txt", 'a')
    if p <= 0.05: f.write(str(pos)+"\t"+str(F)+"\t"+str(p)+'\t'+'+'+"\n")
    if p > 0.05: f.write(str(pos)+"\t"+str(F)+"\t"+str(p)+'\t'+'-'+"\n")
    f.close()
    
    
    
[anova_dmp(*row) for row in df.itertuples(index=False)]
