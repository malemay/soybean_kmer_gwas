# This script combines the information from all reads of all samples
# for a given trait into a GRanges object that can be used for
# downstream applications such as generating Manhattan plots and computing signal ranges

# Loading the required libraries
library(gwask)
library(GenomicRanges)
library(parallel)

phenotype <- commandArgs(trailingOnly = TRUE)[1]

# Setting the working directory
setwd("gwas_results/kmers/")

# Getting the set of bam files to read
# DEPENDENCY: output of katcher for a given trait (gwas_results/kmers/%/katcher_results/KATCHER)
bam_files <- dir(path = paste0(phenotype, "/katcher_results"),
		 pattern = ".*_pvalues_sorted\\.bam$",
		 full.names = TRUE,
		 recursive = TRUE)

if(!length(bam_files)) {
	warning("No input bam files available for phenotype ", phenotype)

	reads <- GRanges()

	reads$qname <- character()
	reads$flag <- integer()
	reads$pos <- integer()
	reads$qwidth <- integer()
	reads$mapq <- integer()
	reads$cigar <- character()
	reads$mrnm <-  character()
	reads$mpos <- integer()
	reads$isize <- integer()
	reads$pvalue <- numeric()
	reads$kmer_canon <- character()
	reads$kmer_pos_count <- integer()
	reads$log10p <- numeric()

	# Setting the seqlevels and seqlengths of the output GRanges from the .fai index
	# DEPENDENCY: refgenome/Gmax_508_v4.0_mit_chlp.fasta
	# DEPENDENCY: refgenome/Gmax_508_v4.0_mit_chlp.fas.fai
	fai_info <- Rsamtools::scanFaIndex("../../refgenome/Gmax_508_v4.0_mit_chlp.fasta")
	GenomeInfoDb::seqlevels(reads) <- GenomeInfoDb::seqlevels(fai_info)                                                                                                     
	GenomeInfoDb::seqlengths(reads) <- GenomeInfoDb::seqlengths(fai_info)
} else {

	message("Processing the following bam files:\n")
	writeLines(bam_files, con = stderr())

	sample_names <- sub("_pvalues_sorted.bam", "", basename(bam_files))

	# Reading the results for several samples at once using parallel::mclapply
	reads <- mclapply(bam_files, 
			  FUN = format_kmer_gwas,
			  ref_fasta = "../../refgenome/Gmax_508_v4.0_mit_chlp.fasta",
			  pattern = NULL,
			  min_mapq = 20,
			  min_count = 10,
			  mc.cores = 10)

	# Once we have the reads we assign a sample column to each of them
	for(i in 1:length(reads)) {
		if(length(reads[[i]])) reads[[i]]$sample <- sample_names[i]
	}

	# Joining all the reads in a common data.frame
	reads <- do.call("c", reads)

	# Then we keep the most common position for each k-mer
	reads <- reads[order(reads$kmer_canon, reads$kmer_pos_count, decreasing = TRUE)]
	reads <- reads[!duplicated(reads$kmer_canon)]
	reads <- sort(reads, ignore.strand = TRUE)

	# There is a special situation if a p-value is so small that -log10(p-value) yields Inf
	# In this case we set the value to the -log10 of the smallest representable value
	reads[reads$pvalue == 0]$log10p <- -log10(.Machine$double.xmin)

}

message("Total of ", length(reads), " unique k-mers found for phenotype ", phenotype)

# Saving the combined reads to file
saveRDS(reads, file = paste0(phenotype, "/katcher_results/", phenotype, "_kmer_positions.rds"))

# Creating a copy in the local directory
file.copy(paste0(phenotype, "/katcher_results/", phenotype, "_kmer_positions.rds"), paste0(phenotype, "_kmer_positions.rds"))

