#!/bin/bash

# Using SVmerge for merging the SVs from all realigned vcf files

# I have used rather stringent parameters so as not to lose information
# from the set of variants; alleles must be sufficiently precise to
# enable genotyping from Illumina data, so we won't merge variants
# that are too different

# Setting the working directory
cd external_data/nanopore_svs/

# Running the command
# DEPENDENCY: refgenome/Gmax_508_v4.0_mit_chlp.fasta
# DEPENDENCY: external_data/nanopore_svs/files.txt
# DEPENDENCY: all normalized Oxford Nanopore SVs from figshare
SVmerge -ref ../../refgenome/Gmax_508_v4.0_mit_chlp.fasta \
	-fof files.txt \
	-prefix svmerged_preliminary \
	-maxdist 15 \
	-reldist 0.2 \
	-relsizediff 0.1 \
	-relshift 0.1 \
	-seqspecific

# Running a script that will preferentially merge realigned variants
# DEPENDENCY: external_data/nanopore_svs/merge_realigned.R
Rscript merge_realigned.R

# Sort vcf entries
bcftools sort -m 500.0M -Ov svmerged.clustered.vcf > svmerged_clustered_sorted.vcf

# This code will prepare the input Oxford Nanopore SV vcf
# so that they are given proper IDs and INFO values (SVTYPE + ACO) to
# be used as input to SVmerge

# Then processing the Oxford Nanopore variants
nanopore_input=svmerged_clustered_sorted.vcf

# For these, we first begin by removing all annotations except SVTYPE
# Then over a single pass with awk we move the original ID to the INFO field
# Then we add the annotation header line with bcftools annotate
# DEPENDENCY: external_data/nanopore_svs/header_lines.txt
bcftools annotate -x "^INFO/SVTYPE" -Ov $nanopore_input | \
	awk 'BEGIN {OFS="\t"} /^#/ {print} !/^#/ {$0=$0 ";ACO=sniffles;InitialIDs=" $3; $3="nanopore_" NR; print}' | \
	bcftools annotate --header-lines header_lines.txt > nanopore_svs.vcf

