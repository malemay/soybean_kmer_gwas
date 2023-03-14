#!/bin/bash

# Setting the working directory
cd illumina_data/

mkdir -p aligned_reads

# DEPENDENCY : refgenome/Gmax_508_v4.0_mit_chlp.fasta
# DEPENDENCY: utilities/correct_sra_metadata.csv
# DEPENDENCY: trimmed FASTQ Illumina reads
cut -d, -f1 utilities/correct_sra_metadata.csv | tail -n+2 | sed 's/"//g' | parallel -j10 '
	# Creating the output directory for sample $i
	mkdir -p aligned_reads/{}
	output_dir=aligned_reads/{}
	refgenome=../refgenome/Gmax_508_v4.0_mit_chlp.fasta

	# Aligning the paired reads
	bwa mem -t 4 $refgenome trimmed_fastq/{}/{}_R1_trimmed.fastq.gz trimmed_fastq/{}/{}_R2_trimmed.fastq.gz | \
		samtools view -bSh - > ${output_dir}/{}_paired.bam

	# Aligning the reads that were unpaired after trimming
	bwa mem -t 4 $refgenome trimmed_fastq/{}/{}_sing_trimmed.fastq.gz | samtools view -bSh - > ${output_dir}/{}_unpaired.bam

	# Sorting the reads
	samtools sort -@ 4 -o ${output_dir}/{}_paired.sort.bam  ${output_dir}/{}_paired.bam
	samtools sort -@ 4 -o ${output_dir}/{}_unpaired.sort.bam  ${output_dir}/{}_unpaired.bam

	# Merging the paired and unpaired reads
	samtools merge ${output_dir}/{}_all.sort.bam ${output_dir}/{}_paired.sort.bam ${output_dir}/{}_unpaired.sort.bam

	# Indexing the reads
	samtools index ${output_dir}/{}_all.sort.bam 

	# Removing the files that are no longer needed
	rm ${output_dir}/{}_paired.bam ${output_dir}/{}_unpaired.bam ${output_dir}/{}_paired.sort.bam ${output_dir}/{}_unpaired.sort.bam'

