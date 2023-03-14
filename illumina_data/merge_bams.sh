#!/bin/bash

# Building a list of commands that will run bamaddrg on all samples with appropriate sample names and read groups
# Sample names are actual cultivar names whereas read groups are the SRR ID of the corresponding bam file

# Setting the working directory
cd illumina_data/

# Making sure that the command list is newly created
> command_list.sh

# DEPENDENCY: utilities/srr_id_correspondence.txt
cat ../utilities/srr_id_correspondence.txt | while read i
do
	# Assigning the sample names and SRR IDs to variables
	bayer_id=$(echo $i | cut -d " " -f1)
	srr_ids=$(echo $i | cut -d " " -f2)
	IFS=';' read -r -a srr_array <<< "$srr_ids"

	# Outputting the mkdir and bamaddrg commands to the file
	printf "bamaddrg " >> command_list.sh

	# Completing the command with all samples
	for id in ${srr_array[@]}
	do
		# DEPENDENCY: reads mapped with bwa
		printf -- "-b aligned_reads/${id}/${id}_all.sort.bam -s $bayer_id -r $id " >> command_list.sh
	done

	# Output redirection and adding a newline character before the next sample
	printf "> merged_bams/${bayer_id}_merged.bam\n" >> command_list.sh
done

cat command_list.sh | parallel -j20 {}

