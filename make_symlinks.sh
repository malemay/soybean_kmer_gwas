#!/bin/bash

for i in $(tail -n+2 utilities/loci.txt | cut -d "," -f1)
do
	trait=$(grep ^$i utilities/loci.txt | cut -d "," -f2)
	
	for j in platypus vg paragraph
	do
		cd gwas_results/${j}
		ln -s ${trait}_gwas.csv ${i}_gwas.csv
		ln -s ${trait}_threshold_5per.txt ${i}_threshold_5per.txt
		cd ../..

	done

	cd gwas_results/kmers
	ln -s ${trait}_kmer_positions.rds ${i}_kmer_positions.rds
	ln -s ${trait}_threshold_5per ${i}_threshold_5per
	cd ../..

done

for i in $(seq 1 $(wc -l utilities/signal_ids.txt | awk '{print $1}'))
do
	line=$(head -n $i utilities/signal_ids.txt | tail -n1)
	id=$(echo $line | cut -d "," -f1)
	trait=$(echo $line | cut -d "," -f2)
	
	for j in platypus vg paragraph
	do
		cd gwas_results/${j}
		ln -s ${trait}_gwas.csv ${id}_gwas.csv
		ln -s ${trait}_threshold_5per.txt ${id}_threshold_5per.txt
		cd ../..

	done

	cd gwas_results/kmers
	ln -s ${trait}_kmer_positions.rds ${id}_kmer_positions.rds
	ln -s ${trait}_threshold_5per ${id}_threshold_5per
	cd ../..

done

