#!/bin/bash

# Setting the working directory
cd illumina_data/

# DEPENDENCY: illumina_data/merged_bams/ILLUMINA_BAM_MERGING
# DEPENDENCY: utilities/srr_id_correspondence.txt
# Processing all samples on 18 cores using parallel
cut -d " " -f1 ../utilities/srr_id_correspondence.txt | \
	parallel -j18 'samtools coverage merged_bams/{}_merged.bam > merged_bams/{}_coverage.txt'

touch SAMTOOLS_COVERAGE
