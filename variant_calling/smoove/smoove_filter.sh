#!/bin/bash

# Setting the current directory
cd variant_calling/smoove/

# Creating a variable for the reference genome
# DEPENDENCY: refgenome/Gmax_508_v4.0_mit_chlp.fasta
refgenome=../../refgenome/Gmax_508_v4.0_mit_chlp.fasta

# Creating a variable for the regions (assembled chromosome pseudomolecules) from which to extract the variants
regions=Gm01,Gm02,Gm03,Gm04,Gm05,Gm06,Gm07,Gm08,Gm09,Gm10,Gm11,Gm12,Gm13,Gm14,Gm15,Gm16,Gm17,Gm18,Gm19,Gm20

# Creating a variable for the file holding smoove variants
# DEPENDENCY: variant_calling/smoove/all_samples.smoove.square.vcf.gz
smoove_file=all_samples.smoove.square.vcf.gz

# We first need to remove some invalid entries whose deletion size is > 0 in the original VCF
bcftools view --exclude "INFO/SVTYPE ~ 'DEL' & INFO/SVLEN > 0" -Oz $smoove_file > all_samples.vcf.gz

# Indexing the vcf file
tabix all_samples.vcf.gz

# Removing variants which are breakends, or larger than +- 500,000 nucleotides, as well as sample fields
bcftools view --regions $regions -G all_samples.vcf.gz -Ov | \
	grep -v "SVTYPE=BND" | \
	bcftools filter --include "INFO/SVLEN>-500000 && INFO/SVLEN<500000 && POS>1" -Ov - > filtered_smoove.vcf

# The header is all mixed-up, so we will put it back in order
# This first command creates a temporary file from which the ##contig header lines have been removed
grep -v "^##contig" filtered_smoove.vcf > tmp_file.vcf

# This second command fixes the header and put the file back into filtered_smoove.vcf
# DEPENDENCY : refgenome/Gmax_508_v4.0_mit_chlp.fasta.fai
bcftools reheader --fai ${refgenome}.fai tmp_file.vcf > filtered_smoove.vcf

# Deleting the temporary file
rm tmp_file.vcf

# Using BayesTyperTools convertAllele to convert the alleles to fully specified sequence 
bayesTyperTools convertAllele -g $refgenome -v filtered_smoove.vcf --keep-imprecise 1 --keep-partial 1 -o converted_smoove

# Using bcftools norm to normalize each of the files
# DEPENDENCY : variant_calling/smoove/ACO_header_line.txt
bcftools norm -f $refgenome -d none -Ou converted_smoove.vcf | bcftools annotate -x "^INFO/SVTYPE" -Ov - | \
	gawk 'BEGIN {OFS = "\t"} /#/ {print $0} !/^#/ {print $0 ";ACO=smoove"}' | \
	bcftools annotate --header-lines ACO_header_line.txt -Ov > normalized_smoove.vcf

# This code extracts the SVs >= 50 nucleotides from the smoove variants
# It also sets the ID to the name of the caller + SV type + line number
# DEPENDENCY : scripts/extract_svs_50.awk
../../scripts/extract_svs_50.awk normalized_smoove.vcf | \
	awk 'BEGIN {OFS="\t"} /^#/ {print} !/^#/ {$3 = NR ; print}' | \
	bcftools annotate --set-id "%INFO/ACO\_%INFO/SVTYPE\_%ID" -Ov - > smoove_filtered.vcf

