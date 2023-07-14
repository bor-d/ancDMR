#!/bin/bash
#SBATCH -p comet

f_prefix=$1
ref_dir=$2
out_prefix=$3

for j in {1..23}
do
	if (($j == 23)); then j=X; fi
	awk 'FNR==NR{exists[$2]++; meth[$2]=$NF/100;next} FNR==1{next};{
	n3=split($NF, sites, ","); 
	d=0;counter=1; for (i=1;i<=length(sites);i++){
  		if(exists[sites[i]]>0){
        	d+=meth[sites[i]]; counter++ }else {print $4, sites[i] > "missing.txt"}}; n4=split($1, abc, ":");
	print abc[2], $6}' ${ref_dir}/chr$j.CpG.positions.txt ${f_prefix}.$j.F > ${out_prefix}_chr$j.txt
	
done 
