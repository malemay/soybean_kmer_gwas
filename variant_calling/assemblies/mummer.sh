#!/bin/bash

# Setting the directories
cd variant_calling/assemblies/
mkdir alignment

for i in $(seq 1 29)
do
	# DEPENDENCY: variant_calling/assemblies/assembly_samples.txt
	sample_info=$(head -n $i assembly_samples.txt | tail -n1)
	sample_path=$(echo $sample_info | cut -d " " -f1)
	sample_id=$(echo $sample_info | cut -d " " -f2)

	# DEPENDENCY: refgenome/Gmax_508_v4.0_mit_chlp.fasta
	# DEPENDENCY: genome assemblies by Liu et al. (2020) in external_data/genome_assemblies/
	nucmer -c 1000 --prefix alignment/${sample_id} --threads 4 \
		../../refgenome/Gmax_508_v4.0_mit_chlp.fasta \
		../../external_data/genome_assemblies/${sample_path}
done

touch MUMMER_ALIGNMENT
