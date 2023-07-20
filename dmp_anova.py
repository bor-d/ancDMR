# -*- coding: utf-8 -*-



import pandas as pd
from scipy import stats as st
import numpy as np
import sys

list_dir = sys.argv[1]
group_file = sys.argv[2]


datafile = list_dir+"filtered_bones.txt"
data = pd.read_csv(datafile, sep = '\t', header=0)
age_df = pd.read_csv(group_file, sep = '\t', header=0)

f = open(list_dir+"anova.txt", 'w')
f.write(str("pos"+"\t"+"anova_stat"+"\t"+"pval" +'\t'+'significant'+"\n"))
f.close()

def anova_dmp(*args):
	temp_df = pd.DataFrame({"sample": data.columns[1:], "methylation": args[1:]})
	merged_df = pd.merge(temp_df, age_df, on="sample", how="left")
    
    F, p = st.f_oneway(pd.to_numeric(temp_df[temp_df["age"] == 'ANC']['value'], errors='coerce'), 
                       pd.to_numeric(temp_df[temp_df["age"] == 'HG']['value'], errors='coerce'), 
                       pd.to_numeric(temp_df[temp_df["age"] == 'MODERN']['value'], errors='coerce'))
    f = open(list_dir+"anova.txt", 'a')
    if p <= 0.05: f.write(str(args[0])+"\t"+str(F)+"\t"+str(p)+'\t'+'yes'+"\n")
    if p > 0.05: f.write(str(args[0])+"\t"+str(F)+"\t"+str(p)+'\t'+'no'+"\n")
    f.close()
    
    
    
[anova_dmp(*row) for row in data.itertuples(index=False)]

# filtering positions for Tukey post hoc test

aov_file = list_dir+"filtered_bones.txt"
aov_df = pd.read_csv(datafile, sep = '\t', header = 0)

aov_df = aov_df[aov_df['significant'] == 'yes']

filtered = pd.merge(data, aov_df, on='pos', how='inner')

# perform Tukey post hoc test

f = open(list_dir+"tukey_res.txt", 'w')
f.write(str("pos"+"\t"+"diff_hg_anc"+"\t"+"pval_hg_anc"+"\t"+"diff_hg_mod"+"\t"+"pval_hg_mod"+'\t'+"diff_anc_mod"+"\t"+'pval_anc_mod'+"\t"+"significant"+"\n"))
f.close()

def tukey_testing(*args):
    try:
        temp_df = pd.DataFrame({'sample': data.columns[1:], 'methylation': args[1:]})
	merged_df = pd.merge(temp_df, age_df, on='sample', how="left")
        temp_df['methylation'] = pd.to_numeric(temp_df['methylation'], errors='coerce')
    
    	# we will use bioinfokit (v1.0.3 or later) for performing tukey HSD test
    	# check documentation here https://github.com/reneshbedre/bioinfokit

    	# perform multiple pairwise comparison (Tukey's HSD)
    	# unequal sample size data, tukey_hsd uses Tukey-Kramer test
        res = stat()
        res.tukey_hsd(df = temp_df, res_var = 'methylation', xfac_var='age', anova_model='methylation ~ C(age)')
    	#print(res.tukey_summary)
    
        pval1, pval2, pval3 = [1, 1, 1]
        diff1, diff2, diff3 = [0, 0, 0]
    
    
        if np.array(res.tukey_summary['p-value']).size > 2:
            pval1, pval2, pval3 = res.tukey_summary['p-value']
            diff1, diff2, diff3 = res.tukey_summary['Diff']
        #print(pval1, pval2, pval3, diff1, diff2, diff3)
        elif np.array(res.tukey_summary['p-value']).size > 1:
            pval1, pval2, pval3 = res.tukey_summary['p-value']
            diff1, diff2, diff3 = res.tukey_summary['Diff']
        
    
        f = open(list_dir+"tukey_res.txt", 'a')
        if (pval1 > 0.05) & (pval2 > 0.05) & (pval3 > 0.05):
            f.write(str(args[0])+"\t"+str(diff1)+"\t"+str(pval1)+"\t"+str(diff2)+"\t"+str(pval2)+'\t'+str(diff3)+"\t"+str(pval3)+"\t"+"none"+"\n")
        elif (pval1 > 0.05) & (pval2 > 0.05) & (pval3 <= 0.05):
            f.write(str(args[0])+"\t"+str(diff1)+"\t"+str(pval1)+"\t"+str(diff2)+"\t"+str(pval2)+'\t'+str(diff3)+"\t"+str(pval3)+"\t"+"anc_vs_mod"+"\n")
        elif (pval1 > 0.05) & (pval2 >= 0.05) & (pval3 > 0.05):
            f.write(str(args[0])+"\t"+str(diff1)+"\t"+str(pval1)+"\t"+str(diff2)+"\t"+str(pval2)+'\t'+str(diff3)+"\t"+str(pval3)+"\t"+"hg_vs_mod"+"\n")
        elif (pval1 <= 0.05) & (pval2 > 0.05) & (pval3 > 0.05):
            f.write(str(args[0])+"\t"+str(diff1)+"\t"+str(pval1)+"\t"+str(diff2)+"\t"+str(pval2)+'\t'+str(diff3)+"\t"+str(pval3)+"\t"+"hg_vs_anc"+"\n")
        elif (pval1 > 0.05) & (pval2 <= 0.05) & (pval3 <= 0.05):
            f.write(str(args[0])+"\t"+str(diff1)+"\t"+str(pval1)+"\t"+str(diff2)+"\t"+str(pval2)+'\t'+str(diff3)+"\t"+str(pval3)+"\t"+"mod_vs_all"+"\n")
        elif (pval1 <= 0.05) & (pval2 > 0.05) & (pval3 <= 0.05):
            f.write(str(args[0])+"\t"+str(diff1)+"\t"+str(pval1)+"\t"+str(diff2)+"\t"+str(pval2)+'\t'+str(diff3)+"\t"+str(pval3)+"\t"+"anc_vs_all"+"\n")
        elif (pval1 <= 0.05) & (pval2 <= 0.05) & (pval3 > 0.05):
            f.write(str(args[0])+"\t"+str(diff1)+"\t"+str(pval1)+"\t"+str(diff2)+"\t"+str(pval2)+'\t'+str(diff3)+"\t"+str(pval3)+"\t"+"hg_vs_all"+"\n")
        elif (pval1 <= 0.05) & (pval2 <= 0.05) & (pval3 <= 0.05):
            f.write(str(args[0])+"\t"+str(diff1)+"\t"+str(pval1)+"\t"+str(diff2)+"\t"+str(pval2)+'\t'+str(diff3)+"\t"+str(pval3)+"\t"+"all_different"+"\n")
        else: f.close()
    
        f.close()
    except (ValueError, TypeError, KeyError):
        f = open(list_dir+"tukey_fail.txt", 'a')
        f.write(str(args[0])+"\t"+"failed"+"\n")
        f.close()
    
    
    
[tukey_testing(*row) for row in filtered.itertuples(index=False)]



    
