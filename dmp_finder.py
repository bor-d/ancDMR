# -*- coding: utf-8 -*-
"""
Created on Wed Jun  1 20:44:41 2022

@author: Shinji
"""

import pandas as pd
#import numpy as np
import scipy.stats as st
#import sys

#directory with files for exact chrom
list_dir = sys.argv[1]
chrom = sys.argv[2]

sample_sheet = pd.DataFrame(pd.read_csv("sample_file.txt", sep='\t', names = ['name', 'org', 'cov']))
header = pd.Series('pos')
header = header.append(sample_sheet['name'])
data = pd.DataFrame(pd.read_csv('shit.txt', sep = '\t', header=None, names=header, nrows = 150)).set_index('pos')

#reading file for mean calculating
datafile = list_dir+"dmp_table.txt"
data = pd.read_csv(datafile, sep = '\t')

#generating output with header
write_ao = "pval_anc_obese.bed"
f = open(write_ao, 'w')
f.write(str("Chrom"+"\t"+"Start"+"\t"+"Stop"+"\t"+"pval_anc_obese"+"\n"))
f.close()

write_lo = "pval_lean_obese.bed"
f = open(write_lo, 'w')
f.write(str("Chrom"+"\t"+"Start"+"\t"+"Stop"+"\t"+"pval_lean_obese"+"\n"))
f.close()

write_al = "pval_anc_lean.bed"
f = open(write_al, 'w')
f.write(str("Chrom"+"\t"+"Start"+"\t"+"Stop"+"\t"+"pval_anc_lean"+"\n"))
f.close()

#main function
def dmp_finder(index, row):
    #try:
        #making df for one position given as arguments
        temp_df = pd.DataFrame(row)
        temp_df = temp_df.reset_index()
        temp_df.rename(columns={index : 'value'}, inplace=True)
        temp_df = pd.concat([temp_df, sample_sheet['org']], axis = 1)
        temp_df = temp_df.set_index(temp_df['org'])
        print(temp_df)
        
        
        #checking distribution
        hg_dist = st.shapiro(temp_df.loc['HG', 'value'])[1] > 0.05
        lean_dist = st.shapiro(temp_df.loc['LEAN', 'value'])[1] > 0.05
        obese_dist = st.shapiro(temp_df.loc['OBESE', 'value'])[1] > 0.05
        
        #checking variances
        if hg_dist and lean_dist: 
            al_var = st.levene(temp_df.loc['HG', 'value'], temp_df.loc['LEAN', 'value'], center='mean')[1] > 0.05
        else:
            al_var = st.levene(temp_df.loc['HG', 'value'], temp_df.loc['LEAN', 'value'])[1] > 0.05
        if lean_dist and obese_dist: 
            ol_var = st.levene(temp_df.loc['LEAN', 'value'], temp_df.loc['OBESE', 'value'], center='mean')[1] > 0.05
        else:
            ol_var = st.levene(temp_df.loc['LEAN', 'value'], temp_df.loc['OBESE', 'value'])[1] > 0.05
        if obese_dist and hg_dist: 
            ao_var = st.levene(temp_df.loc['OBESE', 'value'], temp_df.loc['HG', 'value'], center='mean')[1] > 0.05
        else:
            ao_var = st.levene(temp_df.loc['OBESE', 'value'], temp_df.loc['HG', 'value'])[1] > 0.05
        
        anc_obese = hg_dist*obese_dist*ao_var
        anc_lean = hg_dist*lean_dist*al_var
        lean_obese = obese_dist*lean_dist*ol_var
        
        #if normal then t-test
        if anc_obese:
            pval_ao = st.ttest_ind(temp_df.loc['HG', 'value'], temp_df.loc['OBESE', 'value'])[1]
        #else Mann-Whitney
        else:
            pval_ao = st.mannwhitneyu(temp_df.loc['HG', 'value'], temp_df.loc['OBESE', 'value'], nan_policy = 'omit')[1]
        
        if anc_lean:
            pval_al = st.ttest_ind(temp_df.loc['HG', 'value'], temp_df.loc['LEAN', 'value'])[1]
        #else Mann-Whitney
        else:
            pval_al = st.mannwhitneyu(temp_df.loc['HG', 'value'], temp_df.loc['LEAN', 'value'], nan_policy = 'omit')[1]
        
        if lean_obese:
            pval_lo = st.ttest_ind(temp_df.loc['LEAN', 'value'], temp_df.loc['OBESE', 'value'])[1]
        #else Mann-Whitney
        else:
            pval_lo = st.mannwhitneyu(temp_df.loc['LEAN', 'value'], temp_df.loc['OBESE', 'value'], nan_policy = 'omit')[1]
        
        #writing to file
        f = open(write_ao, 'a')
        f.write(str(chrom)+"\t"+'value'+"\t"+str(index+1)+"\t"+str(pval_ao)+"\n")
        f.close()
       
        f = open(write_lo, 'a')
        f.write(str(chrom)+"\t"+'value'+"\t"+str(index+1)+"\t"+str(pval_lo)+"\n")
        f.close()
        
        f = open(write_al, 'a')
        f.write(str(chrom)+"\t"+'value'+"\t"+str(index+1)+"\t"+str(pval_al)+"\n")
        f.close()
        
    #except (ValueError, TypeError, KeyError):
     #   f = open("fail.txt", 'a')
      #  f.write('value'+"\t"+"failed calculating stats"+"\n")
       # f.close()
 
    
#list comprehension for faster calc    
[dmp_finder(index, row) for index, row in data.iterrows()]

    
