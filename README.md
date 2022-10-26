These scripts were developed for ancient DNA epigenetic analysis  
Main goal: smooth methylation signal in samples and perform statistical analysis to find differentially methylated positions and regions

Scripts are meant to work with single chromosome at a time so genome needs to be splitted first
Computing methylation profiles is meant to be done with DamMet (Kristian Hanghøj, 2019, https://gitlab.com/KHanghoj/DamMet)

Bash script for retrieval methylation levels from DamMet output is meant to be added later

smoothing.py takes methylation levels calculated by DamMet and other mathylation levels that one wants to use (% of meth) and outputs table containing positions as a first column and methylation levels of samples smoothed by calculating running mean within specified sliding window

dmp_finder.py uses this table to calculate t-test or Mann-Whitney statistics for each pair of groups at every position and outputs table with positions and pvalues

After that using combined p-values is recommended 
https://github.com/brentp/combined-pvalues

mean_dmr_final.py uses region_p bed file and output table of smoothing.py to calculate mean methylation within group in discovered regions; it prints table with chrom, region start and end as first three columns and then columns for each group of samples
