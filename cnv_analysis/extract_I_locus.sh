#!/bin/bash

# Extracting the 27 kb of the I locus according to data in Tuteja and Vodkin (2008), p. S-55
# Cluster A of the CHS genes goes from 27,026 to 37,937.
# Cluster B goes from 43,806 to 54,716
# In between is a spacer sequence that we also extract

# Loading the required modules
module load bwa/0.7.17
module load samtools/1.15

# The I locus sequence will then be mapped to the Williams82 reference version 4
# to identify its limits
# DEPENDENCY: external_data/BAC77G7-a.fasta
# DEPENDENCY: external_data/BAC77G7-a.fasta.fai
samtools faidx external_data/BAC77G7-a.fasta BAC77G7a:27026-54716 > cnv_analysis/i_locus.fa

# Using bwa to align this part of the contig to the Williams82 assembly version 4
# DEPENDENCY: refgenome/Gmax_508_v4.0_mit_chlp.fasta
bwa mem refgenome/Gmax_508_v4.0_mit_chlp.fasta cnv_analysis/i_locus.fa | samtools view -b > cnv_analysis/i_locus.bam
samtools index cnv_analysis/i_locus.bam

