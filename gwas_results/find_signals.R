# Loading the required libraries
suppressMessages(library(GenomicRanges))
suppressMessages(library(gwastools))

# Getting the name of the trait and of the program that called the genotypes from the command line
trait <- commandArgs(trailingOnly = TRUE)[1]
program <- commandArgs(trailingOnly = TRUE)[2]

# A parameter that sets the distance between markers for merging signals with extract_signals
signal_distance <- 10^5

# Reading the GWAS results and threshold
# DEPENDENCY: GWAS association results
# DEPENDENCY: GWAS thresholds
gwas_results <- readRDS(paste0("gwas_results/", program, "/", trait, "_gwas.rds"))
threshold <- -log10(as.numeric(readLines(paste0("gwas_results/", program, "/", trait, "_threshold_5per.txt"))))

# Extracting the signals found for that trait
gwas_signals <- extract_signals(gwas_results,
				threshold = threshold,
				distance = signal_distance)

# Saving the signals to an RDS file for retrieval later on
saveRDS(gwas_signals,
	file = paste0("gwas_results/", program, "/", trait, "_signal.rds"),
	compress = FALSE)

# Also creating the appropriate symlinks such that the analyses on each locus can relate to the right trait
signal_ids <- read.table("utilities/signal_ids.txt", header = FALSE, sep = ",")
signal_ids <- signal_ids[signal_ids[[2]] == trait, ]

if(nrow(signal_ids)) {
	for(i in 1:nrow(signal_ids)) {
		signal_file1 <- paste0(signal_ids[i, 2], "_signal.rds")
		signal_file2 <- paste0("gwas_results/", program, "/", signal_ids[i, 1], "_locus_signal.rds")
		if(file.exists(signal_file2)) unlink(signal_file2)
		file.symlink(signal_file1, signal_file2)
	}
}

