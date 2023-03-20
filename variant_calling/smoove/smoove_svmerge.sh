#!/bin/bash

# Using SVmerge for merging the smoove SVs used as input to Paragraph
# This first merging pass is performed with reduced distance
# between the SVs (-maxdist 5 vs the standard -maxdist 15)
# so as to reduce the number of pairwise comparisons at the next step
# when SVs from all 4 SV calling programs will be combined

# Setting the current directory
cd variant_calling/smoove

# Running the command
# DEPENDENCY: refgenome/Gmax_508_v4.0_mit_chlp.fasta
# DEPENDENCY: variant_calling/smoove/smoove_filtered.vcf
$SVmerge -ref ../../refgenome/Gmax_508_v4.0_mit_chlp.fasta \
	-variants smoove_filtered.vcf \
	-prefix smoove_svmerged \
	-maxdist 5 \
	-reldist 0.2 \
	-relsizediff 0.1 \
	-relshift 0.1 

