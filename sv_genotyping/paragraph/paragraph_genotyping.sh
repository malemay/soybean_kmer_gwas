#!/bin/bash

# Setting the working directory
cd sv_genotyping/paragraph/

# DEPENDENCY: refgenome/Gmax_508_v4.0_mit_chlp.fasta
refgenome=../../refgenome/Gmax_508_v4.0_mit_chlp.fasta

# First I do the whole thing on the first sample to generate the json files
for i in $(cut -d " " -f1 ../../ utilities/srr_id_correspondence.txt)
do
	# DEPENDENCY: sv_genotyping/paragraph/MANIFEST_FILES
	# DEPENDENCY: sv_genotyping/paragraph/all_svs_padded.vcf
	# DEPENDENCY: 
	m_opt=$(grep $i manifest_files/${i}_manifest.txt | awk '{print int($3 * 20)}')
	python3 multigrmpy -t20 -M $m_opt -i all_svs_padded.vcf -m manifest_files/${i}_manifest.txt -r $refgenome \
		--scratch-dir tmpdir/ -o ${i}_results
done

