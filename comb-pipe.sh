#!/bin/bash
#SBATCH -p comet
#SBATCH -c 5

for i in {1..23}
do
	if (($i == 23)); then i=X; fi
	
	# filtering Tukey results (leaving hg/modern p < 0.05)	
	awk 'BEGIN {FS = "\t"; OFS = "\t"} NR == 1 || ($5 < 0.05)' chr$i/tukey_res.txt | cut -f1,5 > chr$i/dmps.txt
	
	# adding chr$i to each line to generate bed like format
	awk -v i="$i" 'BEGIN {FS = "\t"; OFS = "\t"} {$1 = $1 "\t" $1; print "chr" i "\t" $0}'
	
	comb-p pipeline -c 4 --seed 5e-2 --dist 200 -p ./chr$i/dmr_chr$i --region-filter-p 0.1 --anno hg19 chr$i/dmps.txt

	awk 'int($5) >= 5' chr$i/dmr_chr$i.regions-p.bed > chr$i/base_dmr.txt
	# filtering regions longer than 5bp
	
	awk '$7 <= 0.05' chr$i/base_dmr.txt > chr$i/long_dmr.txt
	# keeping only significantly changed regions
	
	wc -l chr$i/long_dmr.txt chr$i/dmr_cutoff.txt
	# counting DMRs before and after pval cutoff
	
	bedtools intersect -wo -a chr$i/dmr_cutoff.txt -b chr$i/target_reg.txt > chr$i/dmr_targeted.txt
	# please split targets into chr dirs first

done

# comb-p options

# -c sets column where pval stored
# --seed sets pval to start region
# --dist region to search in for next pval > seed
# -p prefix
# --region-filter post-analysis filter FDR
# --anno reference genome to annot with
