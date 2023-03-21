#!/bin/bash

# Setting the current directory
cd variant_calling/assemblies/

# Creating a variable for the path to the reference genome
# DEPENDENCY: refgenome/Gmax_508_v4.0_mit_chlp.fasta
refgenome=../../refgenome/Gmax_508_v4.0_mit_chlp.fasta

# Converting the output of smvu to VCF using an R script
# DEPENDENCY: variant_calling/assemblies/svmu_to_vcf.R
# DEPENDENCY: variant_calling/assemblies/SVMU_CALLING
# DEPENDENCY: genome assemblies by Liu et al. (2020) in external_data/genome_assemblies/
Rscript svmu_to_vcf.R

# Getting the sample information
# DEPENDENCY: variant_calling/assemblies/assembly_samples.txt
for sample_id in $(cut -d " " -f2 assembly_samples.txt)
do
	# First we excluded all records that have N in either the REF or ALT columns
	echo "Original VCF contains $(grep -cv '^#' ${sample_id}/${sample_id}.vcf) SVs" >&2
	awk 'BEGIN {OFS="\t"} /^#/ {print} $0 !~ /^#/ && $4 !~ /N/ && $5 !~ /N/ {print}' ${sample_id}/${sample_id}.vcf > ${sample_id}/${sample_id}_filtered.vcf
	echo "Filtered VCF contains $(grep -cv '^#' ${sample_id}/${sample_id}_filtered.vcf) SVs" >&2

	# Now we normalize the variant representation using bcftools norm
	bcftools norm -f $refgenome -Ov ${sample_id}/${sample_id}_filtered.vcf > ${sample_id}/${sample_id}_normalized.vcf
done

touch ASSEMBLY_FILTERING
