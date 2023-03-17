#!/bin/bash

# Setting the working directory
cd illumina_data/

# Processing all samples on 10 cores using parallel
# DEPENDENCY: illumina_data/merged_bams/ILLUMINA_BAM_MERGING
# DEPENDENCY: utilities/srr_id_correspondence.txt
cut -d " " -f1 ../utilities/srr_id_correspondence.txt | parallel -j10 'samtools stats merged_bams/{}_merged.bam > merged_bams/{}_stats.txt'

touch SAMTOOLS_STATS
