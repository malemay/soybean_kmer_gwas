#!/bin/bash

# Setting the current directory
cd variant_calling/assemblies/

# Getting the sample information
for i in $(seq 1 29)
do
	# DEPENDENCY: variant_calling/assemblies/assembly_samples.txt
	sample_info=$(head -n $i assembly_samples.txt | tail -n1)
	sample_path=$(echo $sample_info | cut -d " " -f1)
	sample_id=$(echo $sample_info | cut -d " " -f2)

	# Creating a directory for the sample and moving to it
	mkdir -p $sample_id
	cd $sample_id

	# Filtering the delta file with 1-to-1 alignment
	# DEPENDENCY: variant_calling/assemblies/MUMMER_ALIGNMENT
	delta-filter -1 ../alignment/${sample_id}.delta > ${sample_id}_filtered.delta

	# The version of svmu used here has been modified for efficiency and lower memory footprint
	# It has been compiled from the code in branch true-dspr in my local svmu git directory
	# DEPENDENCY: refgenome/Gmax_508_v4.0_mit_chlp.fasta
	# DEPENDENCY: genome assemblies by Liu et al. (2020) in external_data/genome_assemblies/
	svmu ${sample_id}_filtered.delta ../../../refgenome/Gmax_508_v4.0_mit_chlp.fasta ../../../external_data/genome_assemblies/${sample_path} 100 h

	# Moving back to the parent directory
	cd ..
done

touch SVMU_CALLING
