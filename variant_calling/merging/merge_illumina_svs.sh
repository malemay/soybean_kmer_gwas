#!/bin/bash

# Setting the working directory
cd variant_calling/merging/

# DEPENDENCY: refgenome/Gmax_508_v4.0_mit_chlp.fasta
# DEPENDENCY: variant_calling/merging/illumina_merging_files.txt
# DEPENDENCY: variant_calling/asmvar/asmvar_svmerged.clustered.vcf
# DEPENDENCY: variant_calling/manta/manta_svmerged.clustered.vcf
# DEPENDENCY: variant_calling/smoove/smoove_svmerged.clustered.vcf
# DEPENDENCY: variant_calling/svaba/svaba_svmerged.clustered.vcf

# Create symlinks
ln -s variant_calling/asmvar/asmvar_svmerged.clustered.vcf
ln -s variant_calling/manta/manta_svmerged.clustered.vcf
ln -s variant_calling/smoove/smoove_svmerged.clustered.vcf
ln -s variant_calling/svaba/svaba_svmerged.clustered.vcf

SVmerge -ref ../../refgenome/Gmax_508_v4.0_mit_chlp.fasta \
	-fof illumina_merging_files.txt \
	-prefix illumina \
	-maxdist 15 \
	-reldist 0.2 \
	-relsizediff 0.1 \
	-relshift 0.1 

