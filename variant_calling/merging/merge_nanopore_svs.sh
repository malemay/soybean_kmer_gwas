#!/bin/bash

# Using SVmerge for merging the Illumina SVs together with Oxford Nanopore SVs called with Sniffles

# The edlib-aligner is in the path

# Setting the working directory
cd variant_calling/merging/

# First processing the Illumina input file
# DEPENDENCY: variant_calling/merging/illumina.clustered.vcf
illumina_input=illumina.clustered.vcf

# This code will prepare the input Illumina SV vcf
# so that it has proper IDs and INFO values (SVTYPE + ACO) to
# be used as input to SVmerge

# We keep only the ACO, SVTYPE and ClusterIDs INFO ; CLusterIDs is renamed to initial IDs because SVmerge will overwrite it
# We also change the ID to illumina_{row number} to have unique IDs for merging with the Oxford Nanopore dataset
# It will also make it easier to identify which variants come from the Illumina dataset, and which come from the Nanopore dataset
bcftools annotate -x "^INFO/ACO,INFO/SVTYPE,INFO/ClusterIDs" -Ov $illumina_input | sed 's/ClusterIDs/InitialIDs/g' | \
	awk 'BEGIN {OFS="\t"} /^#/ {print} !/^#/ {$3 = "illumina_" NR; print}' > illumina_svs.vcf

# Running SVmerge itself
# DEPENDENCY : refgenome/Gmax_508_v4.0_mit_chlp.fasta
# DEPENDENCY : variant_calling/merging/nanopore_merging_files.txt
# DEPENDENCY : external_data/nanopore_svs/nanopore_svs.vcf
ln -s ../../external_data/nanopore_svs/nanopore_svs.vcf
SVmerge -ref ../../refgenome/Gmax_508_v4.0_mit_chlp.fasta \
	-fof nanopore_merging_files.txt \
	-prefix nanopore_svmerged \
	-maxdist 15 \
	-reldist 0.2 \
	-relsizediff 0.1 \
	-relshift 0.1 

# Running the R code that will output a VCF file where Illumina variants are favoured
# relative to Oxford Nanopore variants whenever a cluster of variants includes variants
# from both sources
# DEPENDENCY: variant_calling/merging/select_svs.R
Rscript select_svs.R

# We sort the vcfs and remove all the annotations that we do not need anymore
bcftools sort -m 4000.0M -Ou illumina_merged.vcf | \
	bcftools annotate -x "^INFO/ACO,INFO/InitialIDs,INFO/SVTYPE,INFO/ClusterIDs" -Ov - > illumina_merged_sorted.vcf

# Using awk to remove the records where the ALT allele contains any N
# The reasoning behind this is that genotyping cannot be properly done
# if the ALT allele is not entirely defined

# I have not removed REF alleles containing N as these may contain
# legitimate deletions of regions that include ambiguous (N) nucleotides
awk 'BEGIN {OFS = "\t"} /^#/ {print} $0 !~ /^#/ && $5 !~ /N/ {print}' illumina_merged_sorted.vcf > usda_svs.vcf

