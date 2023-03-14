#!/bin/bash

# Setting the working directory
cd illumina_data/

# DEPENDENCY: external_data/adapters.fa
# DEPENDENCY: illumina sequencing raw data
# Launching the command on all samples using parallel
cut -d, -f1 ../utilities/correct_sra_metadata.csv | tail -n+2 | sed 's/"//g' | parallel -j8 "
  mkdir -p trimmed_fastq/{}
  $bbduk in1=raw_fastq/{}_1.fastq.gz in2=raw_fastq/{}_2.fastq.gz \
	out1=trimmed_fastq/{}/{}_R1_trimmed.fastq.gz out2=trimmed_fastq/{}/{}_R2_trimmed.fastq.gz outs=trimmed_fastq/{}/{}_sing_trimmed.fastq.gz \
	ref=$adapters literal=AGATCGGAAGAGCACACGTCTGAACTCCAGTCA,AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT \
	stats=trimmed_fastq/{}/{}_stats_trimmed.txt \
	qout=33 ktrim=r k=23 mink=9 hdist=1 tpe tbo \
	threads=4 \
	qtrim=r trimq=10 minlength=35 minavgquality=15"

