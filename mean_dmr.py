# -*- coding: utf-8 -*-



import pandas as pd
import numpy as np
import sys
import warnings

list_dir = sys.argv[1]
chrom = sys.argv[3]
age_file = sys.argv[2]

sample_sheet = pd.DataFrame(pd.read_csv(age_file, sep='\t', header = 0))

datafile = list_dir + "smoothed.txt"
data = pd.read_csv(datafile, sep = '\t', header = 0)
data.set_index('pos', inplace=True)

reg_file = list_dir + "dmr_targeted.txt"
reg_df = pd.DataFrame(pd.read_csv(reg_file, sep = '\t', header = None, 
                                  names = ['chrom', 'start', 'end']), 
                      columns = ['chrom', 'start', 'end'])

dmr_means_file = list_dir + "dmr_means.txt"
f = open(dmr_means_file, 'w')
f.write('Chrom'+'\t'+"Region_start"+"\t"+"Region_end"+"\t"+"anc_mean"+"\t"+"hg_mean"+'\t'+"modern_mean"+'\t'+"delta_meth"+'\t'+"change"+"\n")
f.close()

cpg_means_file = list_dir + "cpg_means.txt"
fl = open(cpg_means_file, 'w')
fl.write("pos"+"\t"+"anc"+"\t"+"hg"+"\t"+"modern"+"\n")
fl.close()

def mean_calc(*args):
        try:
            temp_df = pd.DataFrame({"sample": data.columns[1:], "methylation": args[1:]})
            
            pos = args[0]
                            
            temp_df = pd.merge(temp_df, age_df, on="sample", how="left")
            temp_df = temp_df.set_index(temp_df['age'])
            
            modern = np.nanmean(temp_df.loc['MODERN', 'value'])
            hg = np.nanmean(temp_df.loc['HG', 'value'])
            anc = np.nanmean(temp_df.loc['ANC', 'value'])
        
            if not np.isnan(modern) or not np.isnan(hg):
                with open(cpg_means_file, 'a') as fl:
                	fl.write(str(pos)+"\t"+str(anc)+"\t"+str(hg)+"\t"+str(modern)+"\n")
                
        except KeyError:
            print(str(index)+' too many nans')


[mean_calc(*row) for row in data.itertuples(index=False)]

data_means = pd.DataFrame(pd.read_csv(cpg_means_file, sep = '\t', header = 0))

def reg_mean(reg_start, reg_end):
    
    positions = []
    
    for index, row in data_means.iterrows():
        if row['pos'] >= reg_start and row['pos'] <= reg_end:
            positions.append(row['pos'])
    
    anc_reg_mean = np.nanmean(data_means.query('pos == @positions')['anc'])
    hg_reg_mean = np.nanmean(data_means.query('pos == @positions')['hg'])
    modern_reg_mean = np.nanmean(data_means.query('pos == @positions')['modern'])
    
    delta = abs(hg_reg_mean - modern_reg_mean)
    
    if hg_reg_mean >= modern_reg_mean:
        meth_state = 'hypo'
    else:
        meth_state = 'hyper'
    
    with open(dmr_means_file, 'a') as f:
    	f.write(str(chrom)+'\t'+str(reg_start)+"\t"+str(reg_end)+"\t"+str(anc_reg_mean)+"\t"+str(hg_reg_mean)+'\t'+str(modern_reg_mean)+'\t'+str(delta)+'\t'+str(meth_state)+"\n")


reg_df.sort_values(by = 'start', inplace=True, ignore_index=True)

[reg_mean(reg_start, reg_end) for reg_start, reg_end in zip(reg_df['start'], reg_df['end']) ]


