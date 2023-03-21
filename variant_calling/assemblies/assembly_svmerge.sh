#!/bin/bash

# Using SVmerge for merging the SVs from all VCF files derived from Liu et al. 2020 genomes + Lee assembly obtained from SoyBase

# I have used rather stringent parameters so as not to lose information
# from the set of variants; alleles must be sufficiently precise to
# enable genotyping from Illumina data, so we won't merge variants
# that are too different

# Setting the current directory
cd variant_calling/assemblies/

# Running the command
# DEPENDENCY: refgenome/Gmax_508_v4.0_mit_chlp.fasta
# DEPENDENCY: variant_calling/assemblies/merging_files.txt
# DEPENDENCY: variant_calling/assemblies/ASSEMBLY_FILTERING
$SVmerge -ref ../../refgenome/Gmax_508_v4.0_mit_chlp.fasta \
	-fof merging_files.txt \
	-prefix svmerged \
	-maxdist 15 \
	-reldist 0.2 \
	-relsizediff 0.1 \
	-relshift 0.1 \
	-seqspecific

mv svmerged.clustered.vcf assembly_svs.vcf

