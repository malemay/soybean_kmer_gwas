#!/bin/bash

# Setting the current directory
cd kmers_table/

for i in $(seq 1 389)
do
	# DEPENDENCY: utilities/srr_id_correspondence.txt
	sample=$(head -n $i ../utilities/srr_id_correspondence.txt | tail -n1 | cut -d " " -f1)

	# Creating a directory for the sample and moving to it
	mkdir -p $sample
	cd $sample

	# Listing the files corresponding to that sample
	srr_ids=$(head -n $i ../../utilities/srr_id_correspondence.txt | tail -n1 | cut -d " " -f2)
	IFS=';' read -r -a srr_array <<< "$srr_ids"

	# Making sure that fastq_files.txt is empty to begin with
	> fastq_files.txt

	# Completing the command with all samples
	# DEPENDENCY: illumina_data/BBDUK_TRIMMING
	for id in ${srr_array[@]}
	do
		echo "../../illumina_data/trimmed_fastq/${id}/${id}_R1_trimmed.fastq.gz" >> fastq_files.txt
		echo "../../illumina_data/trimmed_fastq/${id}/${id}_R2_trimmed.fastq.gz" >> fastq_files.txt
		echo "../../illumina_data/trimmed_fastq/${id}/${id}_sing_trimmed.fastq.gz" >> fastq_files.txt
	done

	# Getting the canonized k-mers for that sample
	kmc -t5 -k31 -ci2 @fastq_files.txt ${sample}_kmc_canon ./ 1> kmc_canon.1 2> kmc_canon.2

	# Getting all the k-mers for that sample
	kmc -t5 -k31 -ci0 -b @fastq_files.txt ${sample}_kmc_all ./ 1> kmc_all.1 2> kmc_all.2

	kmers_add_strand_information -c ${sample}_kmc_canon -n ${sample}_kmc_all -k 31 -o ${sample}_kmers_with_strand

	# Going back to the parent directory
	cd ..
done

# Creating a tab-separated file with a list of paths to kmers_with_strand files and sample name
>kmers_list_paths.txt

for i in $(cut -d " " -f 1 ../utilities/srr_id_correspondence.txt)
do
	printf "${i}/${i}_kmers_with_strand\t${i}\n" >> kmers_list_paths.txt
done

# Then generating the list of k-mers to use for the analysis using the parameters suggested in the manual
list_kmers_found_in_multiple_samples -l kmers_list_paths.txt -k 31 --mac 5 -p 0.2 -o kmers_to_use

# Generating the kmers presence/absence table across all individuals
build_kmers_table -l kmers_list_paths.txt -k 31 -a kmers_to_use -o kmers_table

