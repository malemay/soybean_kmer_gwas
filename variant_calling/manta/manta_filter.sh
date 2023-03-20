#!/bin/bash

# Setting the current directory
cd variant_calling/manta/

# 1- Filtering out undesired variants

# Candidate variants will be filtered here to remove breakends (MantaBND) 
# as well as variants on chloroplast, mitochondrion or scaffolds, or larger than 500 000 nucleotides

# This is a list of the reference chromosomes used to filter the variants
regions=Gm01,Gm02,Gm03,Gm04,Gm05,Gm06,Gm07,Gm08,Gm09,Gm10,Gm11,Gm12,Gm13,Gm14,Gm15,Gm16,Gm17,Gm18,Gm19,Gm20

# DEPENDENCY: The candidateSV.vcf.gz files of the 78 batches used to call Manta
# DEPENDENCY: variant_calling/manta/MANTA_CALLING
# Looping over the 78 batches
for i in $(seq 1 78)
do
	manta_file=manta_batch_${i}/results/variants/candidateSV.vcf.gz

	bcftools view --regions $regions -Ou $manta_file | \
		bcftools filter --include "INFO/SVLEN>-500000 && INFO/SVLEN<500000" -Ov - | \
		grep -v "SVTYPE=MantaBND" > manta${i}_SV_filtered.vcf
done


# 2- Using bayesTyperTools to convert the alleles to explicit sequence

# DEPENDENCY : refgenome/Gmax_508_v4.0_mit_chlp.fasta
for i in $(seq 1 78)
do
	bayesTyperTools convertAllele -g ../../refgenome/Gmax_508_v4.0_mit_chlp.fasta -v manta${i}_SV_filtered.vcf \
		--keep-imprecise 1 --keep-partial 1 -o manta${i}_SV_converted
done

# 3- Normalizing the allele sequence representation using bcftools norm

# Using bcftools norm to normalize each of the files, and index them with tabix
for i in $(seq 1 78)
do
	bcftools norm -Oz -f ../../refgenome/Gmax_508_v4.0_mit_chlp.fasta manta${i}_SV_converted.vcf > manta${i}_SV_norm.vcf.gz
	tabix manta${i}_SV_norm.vcf.gz
done

# 4- Merging all Manta normalized files using bcftools merge, and also adjusting the annotations

# DEPENDENCY : variant_calling/manta/ACO_header_line.txt
bcftools merge -m none -Ou $(ls *norm.vcf.gz) | bcftools norm -d none -Ou - | bcftools annotate -x "^INFO/SVTYPE" -Ov - | \
	gawk 'BEGIN {OFS = "\t"} /#/ {print $0} !/^#/ {print $0 ";ACO=manta"}' | \
	bcftools annotate --header-lines ACO_header_line.txt -Ov > manta_merged.vcf

# 5- Extracting SVs > 50 bp from the merged file to generate the final manta SV dataset

# This code extracts the SVs >= 50 nucleotides from the Manta variants
# It also sets the ID to the name of the caller + SV type + line number
# DEPENDENCY : scripts/extract_svs_50.awk
../../scripts/extract_svs_50.awk manta_merged.vcf | \
	awk 'BEGIN {OFS="\t"} /^#/ {print} !/^#/ {$3 = NR ; print}' | \
	bcftools annotate --set-id "%INFO/ACO\_%INFO/SVTYPE\_%ID" -Ov - > manta_filtered.vcf

