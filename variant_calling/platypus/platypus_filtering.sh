#!/bin/bash

# Setting the working directory
cd variant_calling/platypus/

# DEPENDENCY: refgenome/Gmax_508_v4.0_mit_chlp.fasta.fai
fai=../../refgenome/Gmax_508_v4.0_mit_chlp.fasta.fai

# Removing sites with a FILTER tag other than PASS or .
# DEPENDENCY: variant_calling/platypus/platypus_all.vcf
bcftools view --apply-filters .,PASS platypus_all.vcf -Ov > platypus_PASS.vcf

# Removing alleles not found in any genotype calls
# and filtering for a minimum MAF of 0.005
# We first need to add the header on the fly to avoid errors
bcftools reheader --fai $fai platypus_PASS.vcf | bcftools view --trim-alt-alleles -Ou | bcftools view --min-af 0.005 -Ov > platypus_maf.vcf

# Removing variants with more than 80% missing data
bcftools filter -e "F_MISSING > 0.8" -Ov platypus_maf.vcf > platypus_missing.vcf

# Removing variants with more than 20% heterozygosity
bcftools filter -e '(N_PASS(GT="het") / (N_SAMPLES - N_MISSING)) > 0.2' -Ov platypus_missing.vcf > platypus_hetfilter.vcf

# This code uses bcftools to:
# 1- Remove alleles that are not found in any genotype calls
# 2- Split multiallelic sites into biallelic sites
# 3- Filter out the variants with more than 50% missing data
# 4- Filter out the variants at MAF < 0.02
# 5- Set the IDs of the markers for use in pruning

# Applying the filters with bcftools
bcftools view --trim-alt-alleles -Ou platypus_hetfilter.vcf | \
	bcftools norm --multiallelics -any -Ou | \
	bcftools filter --exclude "F_MISSING > 0.5" -Ou | \
	bcftools view --min-af 0.02:minor -Ov | \
	bcftools annotate --set-id '%CHROM\_%POS\_' -Ov | \
	awk '/^#/ {x+=1 ; print}
	     !/^#/ {$3 = $3 "" NR - x; print}' OFS='\t' > platypus_maf02_missing50.vcf

# Using plink to prune the set of SNPs with the following parameters:
# - Processing 1,000 markers at a time
# - Sliding window of 100 markers
# - r2 threshold of 0.9
plink --vcf platypus_maf02_missing50.vcf --allow-extra-chr --indep-pairwise 1000 100 0.9

# Using vcftools to extract the pruned SNPs and indels from the full file
input_file=platypus_maf02_missing50.vcf
pruned_snps=plink.prune.in
vcftools --vcf $input_file --snps $pruned_snps --recode --stdout > platypus_pruned.vcf

# Using gawk to recode the alleles
# We also take this opportunity to remove markers that do not match chromosomes of the form Gm[0-9]{2}
awk 'BEGIN {OFS = "\t"} 
	/^#/ {print} 
	!/^#/ && $1 ~ /^Gm[0-9]{2}$/{
		$4 = "A" 
		$5 = "T" 
		gsub(/Gm0?/,"", $1)
		print}' \
		platypus_pruned.vcf > platypus_formatted.vcf

# Using TASSEL to convert to diploid Hapmap format
# OUTPUT: variant_calling/platypus/platypus_formatted.hmp.txt
./run_pipeline.pl -Xmx20000m -vcf platypus_formatted.vcf -sortPositions -export platypus_formatted.hmp.txt -exportType HapmapDiploid

