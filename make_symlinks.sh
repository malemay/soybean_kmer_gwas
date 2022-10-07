#!/bin/bash

# Creating the symlinks from the main results in ~/sv_gwas/usda_lines/gwas to the local gwas_results directory

# First for the three genotyping programs
for i in platypus vg paragraph
do
	mkdir -p gwas_results/${i}
	cd gwas_results/${i}

	for trait in $(cat ../../utilities/trait_names.txt)
	do
		ln -s ~/sv_gwas/usda_lines/gwas/filtered_gwas/${i}/${trait}/GAPIT.MLM.${trait}.GWAS.Results.csv ${trait}_gwas.csv
		ln -s ~/sv_gwas/usda_lines/gwas/filtered_gwas/${i}/${trait}/${trait}_threshold_5per.txt
	done

	cd ../..
done

# And then for the kmers approach
mkdir -p gwas_results/kmers
cd gwas_results/kmers

for trait in $(cat ../../utilities/trait_names.txt)
do
	# Creating symlinks to the processed results of the k-mer analysis
	ln -s ~/sv_gwas/usda_lines/gwas/filtered_gwas/kmers/${trait}/katcher_results/${trait}_kmer_positions.rds
	ln -s ~/sv_gwas/usda_lines/gwas/filtered_gwas/kmers/${trait}/kmers/threshold_5per ${trait}_threshold_5per
done

cd ../..

# Creating symlinks to the main directories of the k-mer analysis
mkdir -p gwas_results/kmer_data
cd gwas_results/kmer_data

for trait in $(cat ../../utilities/trait_names.txt)
do
	ln -s ~/sv_gwas/usda_lines/gwas/filtered_gwas/kmers/${trait}
done

cd ../..

# Creating symlinks for the filtered VCF files and their indexes
ln -s ~/sv_gwas/usda_lines/gwas/platypus/platypus_pruned.vcf.gz filtered_variants/platypus/filtered_variants.vcf.gz
ln -s ~/sv_gwas/usda_lines/gwas/platypus/platypus_pruned.vcf.gz.tbi filtered_variants/platypus/filtered_variants.vcf.gz.tbi
ln -s ~/sv_gwas/usda_lines/gwas/paragraph/paragraph_filtered.vcf.gz filtered_variants/paragraph/filtered_variants.vcf.gz
ln -s ~/sv_gwas/usda_lines/gwas/paragraph/paragraph_filtered.vcf.gz.tbi filtered_variants/paragraph/filtered_variants.vcf.gz.tbi
ln -s ~/sv_gwas/usda_lines/gwas/vg/vg_filtered.vcf.gz filtered_variants/vg/filtered_variants.vcf.gz
ln -s ~/sv_gwas/usda_lines/gwas/vg/vg_filtered.vcf.gz.tbi filtered_variants/vg/filtered_variants.vcf.gz.tbi

