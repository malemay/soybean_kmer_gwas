#!/bin/bash

# Setting the working directory
cd illumina_data/soysnp50k_genotyping/

# Creating a directory for the output files
mkdir -p calls

# Creating a -T file suitable for input to bcftools pileup & call so they can constraint calling to those alleles
bcftools sort -Ou soysnp50K_gmax_v4.vcf | \
	bcftools query -f '%CHROM\t%POS\t%REF,%ALT\n' - | \
	bgzip -c > alleles.tsv.gz

tabix -s1 -b2 -e2 alleles.tsv.gz

# Creating a bgzip-compressed and indexed version of the SoySNP50K vcf
bcftools sort -Ov soysnp50K_gmax_v4.vcf | bgzip -c > soysnp50K_gmax_v4.vcf.gz
tabix soysnp50K_gmax_v4.vcf.gz

# Creating a file with read groups
> read_groups_tmp.txt
> pi_ids_tmp.txt

# DEPENDENCY: phenotypic_data/phenotypic_data.csv
for i in $(seq 2 $(wc -l ../../phenotypic_data/phenotypic_data.csv | cut -d " " -f1))
do
	line=$(head -n $i ../../phenotypic_data/phenotypic_data.csv | tail -n1)
	pi=$(echo $line | cut -d ";" -f5 | sed 's/"//g')
	bayer=$(echo $line | cut -d ";" -f2 | sed 's/"//g')

	printf "$bayer\t$pi\n" >> pi_ids_tmp.txt
	printf '*\t../merged_bams/'${bayer}'_merged.bam\t'${pi}'\n' >> read_groups_tmp.txt
done

awk '$2 !~ /NA/ {print}' pi_ids_tmp.txt > pi_ids.txt
awk '$3 !~ /NA/ {print}' read_groups_tmp.txt > read_groups.txt

rm pi_ids_tmp.txt read_groups_tmp.txt

