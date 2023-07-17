These scripts were developed for ancient DNA epigenetic analysis  
Main goal: smooth methylation signal in samples and perform statistical analysis to find differentially methylated positions and regions


Scripts are designed to work with single chromosome at a time so genome needs to be splitted first.
Computing methylation profiles is to be done with DamMet (Kristian Hanghøj, 2019, https://gitlab.com/KHanghoj/DamMet)


Step 1. Ancient profiles generation: 

bash (or sbatch) dammet_run_per_sample.sh <path-to-single-BAM-per-sample> <path-to-excluded-vcf> <output-directory-and-sample-prefix>
CpG positions file can be generated by perl script cgpositionFinder.pl and has a name like chr*.CpG.positions.txt

extract_dm.sh <path-to-F-with-name-prefix> <path-to-directory-with-CpG-reference> <out-prefix>

smoothing_dm_adding_modern.py takes methylation levels calculated by DamMet (tab delimeted txt files produced by extract_dm.sh) and other methylation levels (4 columns .bed where column 4 is methylation %) and outputs table containing positions as a first column and methylation levels of samples smoothed by calculating running mean within specified sliding window.
This script requiers two arguments: first is path to directory where all ancient .txt (result of extract_dm.sh) and modern .bed files are stored and second is path to file containing CpG positions of your reference genome.
CpG positions file can be generated by perl script cgpositionFinder.pl and has a name like chr*.CpG.positions.txt


Step 2. Filtering non-tissuespecific CpG-sites:

bash (or sbatch) run_tissue_comparison.sh

Script addresses 2 Python scripts:
- comp_table.py generates a tab delimeted file storing all modern samples methylation with adjustment to several .bed formats
- anova_tissues.py performs ANOVA on the generated table and filters ancient methylation table from step 1. Resulting table contains only nontissue-specific positions

Step 3. DMR search:

bash dmr_search.sh

Script addresses 3 Python scripts and performs calculation of combined p-values (https://github.com/brentp/combined-pvalues) to aggregate DMPs into DMRs:
- dmp_anova.py uses filtered table to calculate ANOVA statistics for three groups (ancient, HGs and modern people) at every CpG site
- dmp_tukey.py performs post hoc tests on CPG sites that significantly changed methylation. Only CpGs with  methylation change between HGs/ancient and modern people are considered DMPs
- comb-p section searches for DMR based on calculated DMPs
- discovered DMRs are intersected with target regions (please, provide the target file)
- mean_dmr.py calculates mean methylation per DMR and outputs calculated means for each group and methylation change direction

Discovered DMRs can be annotated with GREAT or other tools
https://github.com/brentp/combined-pvalues

mean_dmr.py uses region_p bed file and output table of smoothing_dm_adding_modern.py to calculate mean methylation within group in discovered regions; it prints bed format table with chromosome, region start and end as first three columns and then columns for each group of samples
