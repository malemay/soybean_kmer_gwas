#!/bin/bash
cut -d " " -f1 utilities/srr_id_correspondence.txt | parallel -j20 "samtools index illumina_data/merged_bams/{}_merged.bam"

