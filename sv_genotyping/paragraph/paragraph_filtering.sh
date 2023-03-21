#!/bin/bash

# This code uses bcftools to:
# 1- Filter the population-scale SV vcf by setting to missing the genotype calls with DP < 2
# 2- Filter out the variants at MAF < 0.02
# 3- Filter out the variants with more than 50% missing data

# Setting the working directory
cd sv_genotyping/paragraph/

# Sorting and indexing all samples prior to merging
# DEPENDENCY: utilities/srr_id_correspondence.txt
for sample in $(cut -d " " -f1 ../../utilities/srr_id_correspondence.txt)
do
	# Sorting the vcf files
	# DEPENDENCY: sv_genotyping/paragraph/PARAGRAPH_GENOTYPING
	bcftools sort -m 2000M -Oz ${sample}_results/genotypes.vcf.gz > ${sample}_results/sorted_genotypes.vcf.gz

	# Indexing the vcf files
	tabix ${sample}_results/sorted_genotypes.vcf.gz
done

# Merging the sorted genotype calls of the 389 samples genotyped with Paragraph using bcftools
bcftools merge -Oz $(ls *_results/sorted_genotypes.vcf.gz) > paragraph_svs.vcf.gz
tabix paragraph_svs.vcf.gz

# Applying the filter and only then computing the tags to add
bcftools filter --exclude "FORMAT/DP < 2" --set-GTs . -Ou paragraph_svs.vcf.gz | \
	bcftools filter --exclude "F_MISSING > 0.5" -Ou |
	bcftools view --min-af 0.02:minor -Ov > paragraph_filtered.vcf

# Creating a bgzip-compressed version of this file in filtered_variants/paragraph/ to retrieve the variants from their IDs
bgzip -c paragraph_filtered.vcf > ../../filtered_variants/paragraph/filtered_variants.vcf.gz
tabix ../../filtered_variants/paragraph/filtered_variants.vcf.gz
	
# Using gawk to recode the alleles
awk 'BEGIN {OFS = "\t"} 
	/^#/ {print} 
	!/^#/ {
		$4 = "A" 
		$5 = "T" 
		gsub(/Gm0?/,"", $1)
		print}' \
		paragraph_filtered.vcf > paragraph_formatted.vcf

# Using TASSEL to convert to diploid Hapmap format
./run_pipeline.pl -Xmx5000m -vcf paragraph_formatted.vcf -sortPositions \
	-export paragraph_formatted.hmp.txt -exportType HapmapDiploid

