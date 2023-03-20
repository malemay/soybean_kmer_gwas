#!/bin/bash

#Creating a variable for the reference genome location
# DEPENDENCY: refgenome/Gmax_508_v4.0_mit_chlp.fasta
refgenome=../../refgenome/Gmax_508_v4.0_mit_chlp.fasta

# Calling the genotypes in parallel with smoove call
mkdir -p results-smoove
# DEPENDENCY: utilities/srr_id_correspondence.txt
# DEPENDENCY: illumina_data/merged_bams/ILLUMINA_BAM_MERGING
cut -d " " -f1 ../../utilities/srr_id_correspondence.txt | parallel -j 20 "
    smoove call --outdir results-smoove/ --excludechroms mitochondrion,chloroplast --name {} --fasta $refgenome -p 1 --genotype ../../illumina_data/merged_bams/{}_merged.bam" 

# Merging the results with smoove merge
smoove merge --name merged -f $refgenome --outdir ./results-smoove results-smoove/*.genotyped.vcf.gz

# Calling the genotypes with smoove genotype
mkdir -p results-genotyped
cut -d " " -f1 ../../utilities/srr_id_correspondence.txt | parallel -j 20 "
    smoove genotype -d -x -p 1 --name {}_joint --outdir results-genotyped/ --fasta $refgenome --vcf results-smoove/merged.sites.vcf.gz ../../illumina_data/merged_bams/{}_merged.bam"

# Pasting the results of all samples in a single vcf file
smoove paste --name all_samples results-genotyped/*.vcf.gz

