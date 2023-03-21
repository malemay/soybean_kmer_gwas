#!/bin/bash

# Setting the current directory
cd gwas_results/kmers/

# Getting the phenotype (trait) analyzed
phenotype=$1

# First we extract the presence-absence patterns of the k-mers identified as significantly associated at 5% threshold
# Creating a directory to store the data related to the reads matching the kmers
kmerdir=${phenotype}/katcher_results
mkdir -p ${kmerdir}

# First we extract only the k-mer sequence from the file of significantly associated k-mers,
# sorted in decreasing order of significance
# DEPENDENCY: gwas_results/kmers/KMER_GWAS
tail -n+2 ${phenotype}/kmers/pass_threshold_5per | sort -gk9 > ${kmerdir}/pass_threshold_5per_sorted.txt
egrep -o "[A-Z]{31}" ${kmerdir}/pass_threshold_5per_sorted.txt > ${kmerdir}/kmers_list.txt

# Getting the presence/absence table of these kmers
# DEPENDENCY: kmers_table/kmers_table.names
# DEPENDENCY: kmers_table/kmers_table.table
filter_kmers -t ../../kmers_table/kmers_table -k ${kmerdir}/kmers_list.txt -o ${kmerdir}/significant_pav_table.txt

# Looping over all samples to extract the reads containing significant k-mers with katcher
# DEPENDENCY: utilities/srr_id_correspondence.txt
for sample in $(cut -d " " -f1 ../../utilities/srr_id_correspondence.txt)
do
	# A variable for the path to the input bam file
	# DEPENDENCY: illumina_data/merged_bams/ILLUMINA_BAM_MERGING
	input_bam=../../illumina_data/merged_bams/${sample}_merged.bam

	# We check whether that sample was used for that phenotype
	# DEPENDENCY: gwas_results/kmers/KMER_GWAS
	sample_in_pheno=$(cut -f1 ${phenotype}.pheno | tail -n +2 | grep -c "^${sample}\$")

	# If so we extract the k-mers from that sample
	if [ $sample_in_pheno -eq 1 ]
	then

		echo "Sample ${sample} included in analysis of ${phenotype}. Extracting k-mers." 2>&1
	
		# First creating a directory for the output and moving to it
		mkdir -p ${phenotype}/katcher_results/${sample}
		cd ${phenotype}/katcher_results/${sample}

		# Extracting the list of k-mers from the PAV table
		list_kmers ../significant_pav_table.txt $sample > ${sample}_kmers.txt

		# Launching katcher using that list of k-mers
		katcher -i ../../../${input_bam} -k ${sample}_kmers.txt -o ${sample}_kmers.bam -t 8 -b 1600000 -m 300

		# Adding the p-values to the output bam file from katcher
		add_pvalues ${sample}_kmers.bam ../../kmers/pass_threshold_5per > ${sample}_pvalues.bam

		# Sorting and indexing the .bam file with the p-values
		samtools sort ${sample}_pvalues.bam > ${sample}_pvalues_sorted.bam
		samtools index ${sample}_pvalues_sorted.bam

		# Moving back to the starting directory
		cd ../../..

	else
		# Otherwise we skip it as it is not of interest
		echo "Sample ${sample} not included in analysis of ${phenotype}. Skipping." 2>&1
	fi
done

touch ${phenotype}/katcher_results/KATCHER

