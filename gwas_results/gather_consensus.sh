#!/bin/bash

# This script gathers the consensus sequences for assemblies
# based on significant k-mers for a given genomic region

# Loading the samtools module
module load samtools/1.15

# Getting the ID of the locus to plot from the command line
locus=$1

# Using this locus to extract the relevant data
# DEPENDENCY: utilities/kmer_plot_ranges.txt
locus_data=$(grep $locus utilities/kmer_plot_ranges.txt)
trait=$(echo $locus_data | cut -d " " -f2)
region=$(echo $locus_data | cut -d " " -f3)
chrom=$(echo $region | cut -d ":" -f1)

# DEPENDENCY: refgenome/Gmax_508_v4.0_mit_chlp.fasta
reference=refgenome/Gmax_508_v4.0_mit_chlp.fasta

# Creating a fasta file with the consensus sequences based on assemblies
samtools faidx -n 100000 $reference $region | sed 's/>.*/>Williams82/' > gwas_results/kmer_consensus/${locus}_sequences.fa
# Another fasta file with the consensus sequences based on alignment of reads with bwa in cases where assembly fails
samtools faidx -n 100000 $reference $region | sed 's/>.*/>Williams82/' > gwas_results/kmer_consensus/${locus}_bwa_sequences.fa

# DEPENDENCY: aligned assembled sequences for that locus
for i in $(ls gwas_results/kmer_data/${trait}/assemblies/)
do
	if [ -s gwas_results/kmer_data/${trait}/assemblies/${i}/${i}_${trait}_bwa.bam ]
	then
		samtools consensus -l 100000 --show-del no -m simple -r $region \
			gwas_results/kmer_data/${trait}/assemblies/${i}/${i}_${trait}_bwa.bam |
			sed "s/${chrom}/${i}/" >> gwas_results/kmer_consensus/${locus}_sequences.fa
	fi

	samtools consensus -l 100000 --show-del no -m simple -r $region \
		illumina_data/merged_bams/${i}_merged.bam |
		sed "s/${chrom}/${i}/" >> gwas_results/kmer_consensus/${locus}_bwa_sequences.fa
done

