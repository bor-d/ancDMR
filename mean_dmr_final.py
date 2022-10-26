# -*- coding: utf-8 -*-
"""
Created on Wed Jun  1 20:44:41 2022

@author: Shinji
"""

import pandas as pd
import numpy as np
#import sys
import warnings
#lots of warnings when no meth level available (RRBS data for 2 groups)
#suggest ignore this NA warnings
#warnings.filterwarnings("ignore")

list_dir = sys.argv[1]
chrom = sys.argv[2]

sample_sheet = pd.DataFrame(pd.read_csv(list_dir+"sample_file.txt", sep='\t', 
                                          names=["name", "org", "cov"]))
header = pd.Series('pos')
header = header.append(sample_sheet['name'])

datafile = "shit.txt"
data = pd.read_csv(datafile, sep = '\t', header = None, names = header)

reg_file_al = list_dir + "DMR_targeted.txt"
reg_df_al = pd.DataFrame(pd.read_csv(reg_file_al, sep = '\t', header = None, 
                                  names = ['chrom', 'start', 'end']), 
                      columns = ['chrom', 'start', 'end'])
#del reg_df_al['min_p'], reg_df_al['z_p']
#reg_df_al.rename(columns = {'z_sidak_p': 'p_adj'}, inplace=True)
#reg_df_al.drop(reg_df_al.loc[reg_df_al['p_adj'] > 0.05].index, inplace=True)

#reg_file_ao = list_dir + "filtered_ao_pval.txt"
#reg_df_ao = pd.DataFrame(pd.read_csv(reg_file_ao, sep = '\t', header = None, 
#                                  names = ['chrom', 'start', 'end', 'min_p', 'n_probes', 'z_p', 'z_sidak_p']), 
#                      columns = ['chrom', 'start', 'end', 'min_p', 'n_probes', 'z_p', 'z_sidak_p'])
#del reg_df_ao['min_p'], reg_df_ao['z_p']
#reg_df_ao.rename(columns = {'z_sidak_p': 'p_adj'}, inplace=True)
#reg_df_ao.drop(reg_df_ao.loc[reg_df_ao['p_adj'] > 0.05].index, inplace=True)

write_al = "means.txt"
f = open(write_al, 'w')
f.write("Region_start"+"\t"+"Region_end"+"\t"+"anc_mean"+"\t"+"lean_mean"+'\t'+"Obese_mean"+'\t'+"delta_meth"+'\t'+"change"+"\n")
f.close()

write_mean = "mean_pos.txt"
f = open(write_mean, 'w')
f.write("pos"+"\t"+"anc"+"\t"+"lean"+"\t"+"obese"+"\n")
f.close()

def mean_adding(index, row):
        try:
            temp_df = pd.DataFrame(row)
            temp_df = temp_df.reset_index()
            temp_df.rename(columns={index : 'value'}, inplace=True)
        
            pos = temp_df.loc[temp_df['index'] == 'pos', 'value']
            temp_df.drop(temp_df.loc[temp_df['index']=='pos'].index, inplace=True)
            temp_df = temp_df.reset_index()
            del temp_df['level_0']
                
            temp_df = pd.concat([temp_df, sample_sheet['org']], axis = 1)
            temp_df = temp_df.set_index(temp_df['org'])
            
            obese = np.nanmean(temp_df.loc['OBESE', 'value'])
            lean = np.nanmean(temp_df.loc['LEAN', 'value'])
            anc = np.nanmean(temp_df.loc['HG', 'value'])
        
        #print (temp_df)
        
            f = open(write_mean, 'a')
            f.write(str(pos[0])+"\t"+str(anc)+"\t"+str(lean)+"\t"+str(obese)+"\n")
            f.close()
        except KeyError:
            print(str(index)+' too many nans')

#firstly calc means for each group at a position
[mean_adding(index, row) for index, row in data.iterrows()]

means_file = "mean_pos.txt"
data_means = pd.DataFrame(pd.read_csv(means_file, sep = '\t', header = 0), 
                      columns = ['pos', 'anc', 'lean', 'obese']).dropna(axis='index')

def reg_mean(reg_start, reg_end):
    
    positions = []
    
    for index, row in data_means.iterrows():
        #print(row['pos'])
        if row['pos'] >= reg_start and row['pos'] <= reg_end:
            positions.append(row['pos'])
            #print(positions, row)
    
    anc_reg_mean = np.nanmean(data_means.query('pos == @positions')['anc'])
    lean_reg_mean = np.nanmean(data_means.query('pos == @positions')['lean'])
    obese_reg_mean = np.nanmean(data_means.query('pos == @positions')['obese'])
    
    #if flag == 'al':
    delta = abs(anc_reg_mean - obese_reg_mean)
    
    if anc_reg_mean >= obese_reg_mean:
        meth_state = 'hypo'
    else:
        meth_state = 'hyper'
    
    f = open(write_al, 'a')
    f.write(str(reg_start)+"\t"+str(reg_end)+"\t"+str(anc_reg_mean)+"\t"+str(lean_reg_mean)+'\t'+str(obese_reg_mean)+'\t'+str(delta)+'\t'+str(meth_state)+"\n")
    f.close()
    #elif flag == 'ao':
    #    delta = abs(anc_reg_mean - obese_reg_mean)
    #    f = open(write_ao, 'a')
    #   f.write(str(reg_start)+"\t"+str(reg_end)+"\t"+str(anc_reg_mean)+"\t"+str(obese_reg_mean)+'\t'+str(delta)+'\t'+str(reg_pval)+'\t'+str(n_probes)+"\n")
    #   f.close()

#sort comb-p output
reg_df_al.sort_values(by = 'start', inplace=True, ignore_index=True)
#reg_df_ao.sort_values(by = 'start', inplace=True, ignore_index=True)

#calculate region mean from position means

[reg_mean(reg_start, reg_end) for reg_start, reg_end in zip(reg_df_al['start'], reg_df_al['end']) ]
                                                                                                   # reg_df_al['p_adj'], reg_df_al['n_probes'])]

#[reg_mean(reg_start, reg_end, reg_pval, n_probes, 'ao') for reg_start, reg_end, reg_pval, n_probes in zip(reg_df_ao['start'], reg_df_ao['end'], 
#                                                                                                    reg_df_ao['p_adj'], reg_df_ao['n_probes'])]
