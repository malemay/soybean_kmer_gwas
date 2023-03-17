#!/bin/bash

# Setting the working directory

# Compressing and indexing the file, otherwise the next command won't work
# DEPENDENCY: variant_calling/merging/candidate_svs.vcf
# DEPENDENCY: sv_genotyping/paragraph/truncate_variants.awk
./truncate_variants.awk ../../variant_calling/merging/candidate_svs.vcf > all_svs_truncated.vcf

# Padding the variants as required by Paragraph
# DEPENDENCY: refgenome/Gmax_508_v4.0_mit_chlp.fasta
# DEPENDENCY: sv_genotyping/paragraph/addMissingPaddingGmax4.py
python addMissingPaddingGmax4.py all_svs_truncated.vcf > all_svs_padded.vcf 

