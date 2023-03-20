#!/bin/bash

# Using SVmerge for merging the SVs identified from Illumina and nanopore data (usda_svs.vcf)
# and those identified from the analysis of the assemblies of Liu et al. (2020) (assembly_svs.vcf)
# with svmu and svmutools

# The edlib-aligner is in the path

# Setting the working directory
cd variant_calling/merging/

# DEPENDENCY : refgenome/Gmax_v4/Gmax_508_v4.0_mit_chlp.fasta
# DEPENDENCY: variant_calling/merging/assembly_merging_files.txt
# DEPENDENCY: variant_calling/merging/usda_svs.vcf
# DEPENDENCY: variant_calling/merging/assembly_svs.vcf
SVmerge -ref ../../refgenome/Gmax_508_v4.0_mit_chlp.fasta \
	-fof assembly_merging_files.txt \
	-prefix assembly_svmerged \
	-maxdist 15 \
	-reldist 0.2 \
	-relsizediff 0.1 \
	-relshift 0.1 

# Changing the file name
mv assembly_svmerged.clustered.vcf candidate_svs.vcf

