#!/bin/bash

# Creating the symlinks from the main results in ~/sv_gwas/usda_lines/gwas to the local gwas_results directory

# First for the three genotyping programs
for i in platypus vg paragraph
do
	mkdir -p gwas_results/${i}
	cd gwas_results/${i}

	for trait in $(cat ../../utilities/trait_names.txt)
	do
		ln -s ~/sv_gwas/usda_lines/gwas/${i}/${trait}/GAPIT.MLM.${trait}.GWAS.Results.csv ${trait}_gwas.csv
		ln -s ~/sv_gwas/usda_lines/gwas/${i}/${trait}/${trait}_threshold_5per.txt
	done

	cd ../..
done

# And then for the kmers approach
mkdir -p gwas_results/kmers
cd gwas_results/kmers

for trait in $(cat ../../utilities/trait_names.txt)
do
	ln -s ~/sv_gwas/usda_lines/gwas/kmers/${trait}/katcher_results/${trait}_kmer_positions.rds
	ln -s ~/sv_gwas/usda_lines/gwas/kmers/${trait}/kmers/threshold_5per ${trait}_threshold_5per
done

