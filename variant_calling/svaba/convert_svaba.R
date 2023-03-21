# Loading the GenomicRanges and parallel packages
library(GenomicRanges)
library(parallel)

# Sourcing the functions used to process svaba output
# DEPENDENCY : scripts/svaba_process.R
source("../../scripts/svaba_process.R")

# Creating a list of IDs over which to iterate
# DEPENDENCY : ../../../srr_id_correspondence.txt
ids <- read.table("../../srr_id_correspondence.txt", header = FALSE, stringsAsFactors = FALSE, col.names = c("cultivar", "srr_ids"))
ids <- ids$cultivar

# Using parallel::mclapply to run on 6 cores in parallel
# DEPENDENCY : SvABA SV vcf files
mclapply(ids, 
	 FUN = function(i) svaba_classifier(vcf_file = paste0(i, ".svaba.sv.vcf"),
					    output_file = paste0(i, ".classified.vcf")),
	 mc.cores = 6)

# DEPENDENCY : refgenome/Gmax_508_v4.0_mit_chlp.fasta
mclapply(ids,
	 FUN = function(i) svaba_converter(vcf_file = paste0(i, ".classified.vcf"),
					   output_file = paste0(i, ".converted.vcf"),
					   refgenome = "../../refgenome/Gmax_508_v4.0_mit_chlp.fasta",
					   contig_pattern = "^Gm[0-9]{2}$",
					   max_span = 500000),
	 mc.cores = 6)

