#!/bin/bash

# This script uses bcftools call and bcftools mpileup to compute genotypes from WGS
# data and then uses bcftools gtcheck to compare the genotypes with SoySNP50K calls

# Setting the working directory
cd illumina_data/soysnp50k_genotyping/

# Creating a variable for the reference genome
# DEPENDENCY : refgenome/Gmax_508_v4.0_mit_chlp.fasta
refgenome=../../refgenome/Gmax_508_v4.0_mit_chlp.fasta

# Getting the name of the sample
# DEPENDENCY: illumina_data/soysnp50k_genotyping/pi_ids.txt
for sample in $(cut -f1 pi_ids.txt)
do

	# Creating a string of options to make coding clearer
	opt="-G read_groups.txt -T alleles.tsv.gz --no-BAQ --min-BQ 15 --min-MQ 20 --skip-indels"

	# Calling the SoySNP50K SNPs in all WGS samples
	# DEPENDENCY: mapped Illumina reads for all samples
	bcftools mpileup -f $refgenome $opt -Ou ../merged_bams/${sample}_merged.bam \
		| bcftools call -m -Ov -T alleles.tsv.gz -C alleles --insert-missed --keep-alts - > calls/${sample}_calls.vcf

	# Creating a directory for the output
	mkdir -p gtcheck
	mkdir -p calls

	bgzip -c calls/${sample}_calls.vcf > calls/${sample}_calls.vcf.gz
	tabix calls/${sample}_calls.vcf.gz

	# Running bcftools gtcheck to compare the results of the WGS data and the SoySNP50K data
	# DEPENDENCY: illumina_data/soysnp50k_genotyping/soysnp50K_gmax_v4.vcf
	bbgzip soysnp50K_gmax_v4.vcf
	bcftools gtcheck -G 1 -g soysnp50K_gmax_v4.vcf.gz calls/${sample}_calls.vcf.gz > gtcheck/${sample}_gtcheck.txt

done

touch GTCHECK
