# -*- coding: utf-8 -*-




import pandas as pd
#import statsmodels.api as sm
from scipy import stats as st
import numpy as np
#import scipy as sp
import sys

list_dir = sys.argv[1]
group_file = sys.argv[2]

datafile = list_dir+"modern_compare.txt"
data = pd.read_csv(datafile, sep = '\t', header = 0)
                          
tissue_df = pd.read_csv(group_file, sep = '\t', header=0)

f = open(list_dir+"aov_tissuespec.txt", 'w')
f.write(str("pos"+"\t"+"anova"+"\t"+"pval" +'\t'+'significant'+"\n"))
f.close()

def anova(*args):
    temp_df = pd.DataFrame({"sample": data.columns[1:], "methylation": args[1:]})
    merged_df = temp_df.set_index('sample').combine_first(tissue_df.set_index('sample'))
    
    F, p = st.f_oneway(*[pd.to_numeric(merged_df[merged_df["tissue"] == t]['methylation'], errors='coerce') for t in ['M', 'F', 'B']])
    
    return [args[0], F, p, 'yes' if p < 0.05 else 'no']

anova_results = [anova(*row) for row in data.itertuples(index=False)]

with open(list_dir + "aov_tissuespec.txt", 'a') as f:
    for result in anova_results:
        f.write("\t".join(str(val) for val in result) + "\n")


# filtering non-tissuespecific regions

anc_file = list_dir+"smoothed.txt"
anc_data = pd.read_csv(anc_file, sep = '\t', header = 0)

spec_file = list_dir+"aov_tissuespec.txt"
spec_data = pd.read_csv(spec_file, sep = '\t', header = 0)

spec_data = spec_data[spec_data['significant'] == 'yes']

filtered = pd.merge(anc_data, spec_data, on='pos', how='inner')

final_file = list_dir+"filtered_bones.txt"
    
filtered.to_csv(final_file, sep='\t', index=False, float_format='%1.7f')
           
           
