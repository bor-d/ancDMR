#!/usr/bin/env bash
#SBATCH -p comet

BAM=$1
# N positions to estimate deamination rates from each end (possibly read length)
MAX_POS=60
# path to reference genome
FASTA=/mnt/storage/lab3/Borodko/ref/human_g1k_v37.fasta
# chromosome to analyze
CHROM=1
# minimum mapping quality
minMapQ=20
# minimum base quality
minBaseQ=20
# global prior on the fraction of methylated CpG's
M=0.75
# CpGs to exclude due to low coverage
EXCLUDE=$2
# N cpgs to include per window for estimating f (the more - the smoother is profile, better coverage - less N)
NCPG=25
#OUT PREFIX
OUT_PREFIX=$3
# Using precalculated Deamination rates as the provided bam only covers a window of 30kb bp.
#PRECALCDEAM=test.deamrates

for i in $(seq 1 1 23)
do
	if (($i==23)) then i=X fi
	CHROM=$i
	OUT=${OUT_PREFIX}_chr$i
	../DM estDEAM -b ${BAM} -r ${FASTA} -c ${CHROM} -q ${minMapQ} -Q ${minBaseQ} -E $EXCLUDE -P ${MAX_POS} -M ${M} -O ${OUT} -skip -verbose
	../DM estF -b ${BAM} -r ${FASTA} -c ${CHROM} -q ${minMapQ} -Q ${minBaseQ} -E $EXCLUDE -P ${MAX_POS} -M ${M} -O ${OUT} -N ${NCPG} -verbose -skip_empty_cpg

done

# launching methylation extractor script - command for slurm cluster
#sbatch compare_dm.sh
bash extract_dm.sh
