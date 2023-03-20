#!/bin/bash

# This script creates the configuration files for all Manta batches
# and then launches those batches

# Setting the current directory
cd variant_calling/manta/

# Running the script the creates the configurations for all batches
# DEPENDENCY: refgenome/Gmax_508_v4.0_mit_chlp.fasta
# DEPENDENCY: variant_calling/manta/config_commands.sh
# DEPENDENCY: variant_calling/manta/manta_config.txt
# DEPENDENCY: illumina_data/merged_bams/ILLUMINA_BAM_MERGING
./config_commands.sh

# Launching the analysis in parallel on all 78 batches
seq 1 78 | parallel -j30 "
cd manta_batch_{}
./runWorkflow.py -j 1 -g 10"

