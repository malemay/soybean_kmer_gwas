#!/bin/bash

# Setting the current directory
cd kmers_table/

# Generating the kmers presence/absence table across all individuals
# DEPENDENCY: kmers_table/kmers_table.table
# DEPENDENCY: kmers_table/kmers_table.names
emma_kinship_kmers -k 31 -t kmers_table --maf 0.02 > kmers_table.kinship

