#!/bin/bash

# Set the working directory
cd variant_calling/platypus/

# DEPENDENCY: refgenome/Gmax_508_v4.0_mit_chlp.fasta
refgenome=../../refgenome/Gmax_508_v4.0_mit_chlp.fasta

# Necessary to avoid premature termination (allows more open files)
ulimit -S -n 8192

# DEPENDENCY: BAM files aligned with BWA and merged with bamaddrg (illumina_data/merged_bams/ILLUMINA_BAM_MERGING)
# Creating the list of bam_files for platypus to run
ls ../../illumina_data/merged_bams/*bam > all_bam_files.txt

# DEPENDENCY: variant_calling/platypus/
Platypus.py callVariants --bamFiles=all_bam_files.txt --nCPU=20 --refFile=${refgenome} --skipRegionsFile=chloroplast,mitochondrion \
	--logFileName=platypus_all.log --bufferSize=10000 --maxReads=10000000 --minReads=10 --maxReadLength=150 \
	--maxSize=500 --minMapQual=20 --minBaseQual=20 --maxVariants=15 --filterReadPairsWithSmallInserts=0 \
	--minVarFreq=0.01 --output=platypus_all.vcf

