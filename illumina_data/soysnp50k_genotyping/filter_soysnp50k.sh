#!/bin/bash

# This script extracts the lines used in our study from the VCF file of the SoySNP50K data
# and filters the SoySNP50K subset such that only SNPs with a MAF > 0.1 are kept

# Setting the working directory
cd illumina_data/soysnp50k_genotyping/

# Getting the names of the accessions in our phenotypic dataset
# DEPENDENCY: phenotypic_data/phenotypic_data.csv
awk -F ";" '{print $5}' ../../phenotypic_data/phenotypic_data.csv | \
	sed 's/"//g' | \
	grep -v "^ACCESSION" | \
	grep -v "^NA$" > phenotypic_accessions.txt

# Using that file to extract the columns of interest from the SoySNP50K VCF file
# DEPENDENCY: external_data/soysnp50k_wm82.a2_41317.vcf.gz 
vcftools --gzvcf ../../external_data/soysnp50k_wm82.a2_41317.vcf.gz --keep phenotypic_accessions.txt --recode --stdout > soysnp50K_subset.vcf

# Using vcftools to filter based on minor allele frequency
vcftools --vcf soysnp50K_subset.vcf --maf 0.1 --recode --stdout | awk '/^#/ {print} /^Chr/ {print}' > soysnp50K_maf10.vcf

