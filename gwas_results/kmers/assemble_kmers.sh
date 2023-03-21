#!/bin/bash

# Setting the current directory
cd gwas_results/kmers/

# DEPENDENCY: gwas_results/kmers/assembly_params.txt
# Setting the trait and getting the ID of the sample to analyze
for i in $(seq 1 $(wc -l assembly_params | cut -d " " -f1))
do
	params=$(head -n $i assembly_params.txt | tail -n 1)
	trait=$(echo $params | cut -d " " -f1)
	sample=$(echo $params | cut -d " " -f2)

	# Path to the executable and reference genome
	# DEPENDENCY: refgenome/Gmax_508_v4.0_mit_chlp.fasta
	reference=../../refgenome/Gmax_508_v4.0_mit_chlp.fasta

	# Creating the output directory if it does not already exist
	outdir=${trait}/assemblies/${sample}
	mkdir -p ${outdir}

	# Extracting the aligned read pairs that matched a significant k-mer
	# DEPENDENCY: katcher output for all traits
	# DEPENDENCY: illumina_data/merged_bams/ILLUMINA_BAM_MERGING
	extract_qname ${trait}/katcher_results/${sample}/${sample}_pvalues_sorted.bam \
		../../illumina_data/merged_bams/${sample}/${sample}_merged.bam > ${outdir}/${sample}_${trait}_significant_pairs.bam

	# Extracting those reads in fastq format with 1 read for the first read in pair and another for the second; singleton reads will not be used
	samtools sort -n ${outdir}/${sample}_${trait}_significant_pairs.bam | \
		samtools fastq -1 ${outdir}/${sample}_${trait}_R1.fastq.gz -2 ${outdir}/${sample}_${trait}_R2.fastq.gz -s /dev/null -0 /dev/null

	# Using those reads for assembly with SPAdes
	spades.py -t1 -m15 --careful \
		-1 ${outdir}/${sample}_${trait}_R1.fastq.gz \
		-2 ${outdir}/${sample}_${trait}_R2.fastq.gz \
		-o ${outdir}/${sample}_${trait}_assembly

	# Aligning with bwa
	bwa mem $reference ${outdir}/${sample}_${trait}_assembly/scaffolds.fasta | \
		samtools view -b | samtools sort > ${outdir}/${sample}_${trait}_bwa.bam

	samtools index ${outdir}/${sample}_${trait}_bwa.bam
done

touch KMER_ASSEMBLIES
