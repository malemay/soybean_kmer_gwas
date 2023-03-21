#!/bin/bash

# Setting the current directory
cd variant_calling/svaba/

# Creating a variable for the reference fasta
# DEPENDENCY: refgenome/Gmax_508_v4.0_mit_chlp.fasta
refgenome=../../refgenome/Gmax_508_v4.0_mit_chlp.fasta

# Creating a variable for the regions to keep from the indel files
regions=Gm01,Gm02,Gm03,Gm04,Gm05,Gm06,Gm07,Gm08,Gm09,Gm10,Gm11,Gm12,Gm13,Gm14,Gm15,Gm16,Gm17,Gm18,Gm19,Gm20

# Running the R script that classifies the SvABA variants and sets their sequences
# DEPENDENCY: variant_calling/svaba/convert_svaba.R
# DEPENDENCY: scripts/svaba_process.R
# DEPENDENCY: variant_calling/svaba/SVABA_CALLING
Rscript convert_svaba.R

# Looping over all the samples
# DEPENDENCY: utilities/srr_id_correspondence.txt
for i in $(cut -d " " -f1 ../../utilities/srr_id_correspondence.txt)
do
	# Removing the contig headers from each of the files and changing the Number for PL in the header from . to G
	grep -v '^##contig' ${i}.converted.vcf  | sed 's/##FORMAT=<ID=PL,Number=\./##FORMAT=<ID=PL,Number=G/' > ${i}_sv_tmp.vcf

	# For the indels, we also add the annotation for the SVTYPE
	# DEPENDENCY : variant_calling/svaba/annotate_svtype.awk
	grep -v '^##contig' ${i}.svaba.indel.vcf | sed 's/##FORMAT=<ID=PL,Number=\./##FORMAT=<ID=PL,Number=G/' | ./annotate_svtype.awk > ${i}_indel_tmp.vcf

	# Adding the headers
	bcftools reheader --fai ${refgenome}.fai ${i}_sv_tmp.vcf > ${i}_sv_header.vcf
	bcftools reheader --fai ${refgenome}.fai ${i}_indel_tmp.vcf > ${i}_indel_header.vcf

	# Compressing indel_header.vcf so I can filter out the scaffolds and organellar genomes
	bgzip ${i}_indel_header.vcf
	tabix ${i}_indel_header.vcf.gz

	# Processing the result through bcftools norm, and keeping only genotypes and the SV annotation
	bcftools norm -f $refgenome -cwx -Ou ${i}_sv_header.vcf | bcftools view -G -Ou - | bcftools annotate -x '^INFO/SVTYPE' -Ov - > ${i}_sv_norm.vcf
	bcftools view -G --regions $regions -Ou ${i}_indel_header.vcf.gz | bcftools norm -f $refgenome -Ou - | bcftools annotate -x '^INFO/SVTYPE' -Ov - > ${i}_indel_norm.vcf

	# Removing the temporary files
	rm ${i}_sv_tmp.vcf ${i}_indel_tmp.vcf ${i}_sv_header.vcf ${i}_indel_header.vcf.gz ${i}_indel_header.vcf.gz.tbi
done

# Zipping and indexing all the files with bgzip and tabix so we can use bcftools merge on them
ls *indel_norm.vcf | parallel -j2 "bgzip {}"
ls *sv_norm.vcf | parallel -j2 "bgzip {}"

ls *indel_norm.vcf.gz | parallel -j2 "tabix {}"
ls *sv_norm.vcf.gz | parallel -j2 "tabix {}"

# Merging all SvABA normalized files using bcftools merge, and also adjusting the annotations
# DEPENDENCY : variant_calling/svaba/ACO_header_line.txt
bcftools merge -m none -Ou $(ls *sv_norm.vcf.gz) $(ls *indel_norm.vcf.gz) | bcftools norm -d none -Ov - | \
	gawk 'BEGIN {OFS = "\t"} /^#/ {print $0} !/^#/ {print $0 ";ACO=svaba"}' | \
	bcftools annotate --header-lines ACO_header_line.txt -Ov - > svaba_merged.vcf

# This code extracts the SVs >= 50 nucleotides from the svaba variants
# It also sets the ID to the name of the caller + SV type + line number
# DEPENDENCY : scripts/extract_svs_50.awk
../../scripts/extract_svs_50.awk svaba_merged.vcf | \
	awk 'BEGIN {OFS="\t"} /^#/ {print} !/^#/ {$3 = NR ; print}' | \
	bcftools annotate --set-id "%INFO/ACO\_%INFO/SVTYPE\_%ID" -Ov - > svaba_filtered.vcf

