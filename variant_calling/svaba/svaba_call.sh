#!/bin/bash

# Creating variables for the svaba executable and the location of the reference genome
# DEPENDENCY: refgenome/Gmax_508_v4.0_mit_chlp.fasta
refgenome=../../refgenome/Gmax_508_v4.0_mit_chlp.fasta

# Running 5 svaba tasks of 4 threads each in parallel
# DEPENDENCY: utilities/srr_id_correspondence.txt
# DEPENDENCY: illumina_data/merged_bams/ILLUMINA_BAM_MERGING
cut -d " " -f1 ../../utilities/srr_id_correspondence.txt | parallel -j5 "svaba run -t ../../illumina_data/merged_bams/{}_merged.bam -p 4 -a {} -G $refgenome --germline -I -L 6 --verbose 1"

touch SVABA_CALLING
