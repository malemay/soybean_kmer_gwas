# Loading the required libraries
suppressMessages(library(grid))
suppressMessages(library(gwastools))
suppressMessages(library(GenomicRanges))

# Getting the name of the trait-locus combination and of the program that called the genotypes from the command line
id <- commandArgs(trailingOnly = TRUE)[1]
program <- commandArgs(trailingOnly = TRUE)[2]

# A parameter that sets the padding of the x-scale when the signal consists of only one marker
xscale_padding <- 10

# Getting the ID of the gene associated with the locus
# DEPENDENCY: utilities/all_signals.rds
target_signal <- readRDS("utilities/all_signals.rds")
target_signal <- target_signal[id]

stopifnot(length(target_signal) == 1)

gene_name <- target_signal$gene_name_v4
### ----------

# Getting the set of reference genes
# DEPENDENCY: refgenome/*.rds
genes <- readRDS("refgenome/gmax_v4_genes.rds")

# Handling the special case where the gene_name value contains more than one gene or no gene at all
if(!is.na(gene_name) && grepl(";", gene_name)) {
	gene_names <- strsplit(gene_name, ";")[[1]]
	gene <- genes[gene_names]
	gene <- reduce(gene, min.gapwidth = 10^6, ignore.strand = TRUE)
} else if (!is.na(gene_name)) {
	gene <- genes[gene_name]
} else {
	gene <- NULL
}

# Reading the GWAS results, signals, and region of top markers from file
gwas_results <- readRDS(paste0("gwas_results/", program, "/", id, "_gwas_locus.rds"))
gwas_signals <- readRDS(paste0("gwas_results/", program, "/", id, "_signal_locus.rds"))
top_markers <- readRDS(paste0("gwas_results/", program, "/", id, "_top_markers.rds"))

# Extracting the signal for the locus of interest from those just read
signal <- subsetByOverlaps(gwas_signals, target_signal, ignore.strand = TRUE)

if(!length(signal)) {
	warning("No signal found by ", program, " for signal ID ", id)
	ptx_plot <- grid::gTree(children = gList(textGrob(paste0("No signal found by ", program, " for signal ID ", id))))
} else {
	# Setting the x-scale while handling special cases where there is no signal or the signal is only 1 bp wide
	if(length(signal) == 1 && width(signal) == 1) {
		start(signal) <- start(signal) - xscale_padding
		end(signal)   <- end(signal) + xscale_padding
	}

	# Creating a grob representing the p-value plot and transcript
	ptx_plot <- pvalueGrob(gwas_results = gwas_results,
			       interval = signal,
			       feature = gene,
			       col = "blue",
			       pruned_col = if(program == "platypus") "firebrick2" else NULL,
			       shading = top_markers,
			       yexpand = c(0.1, 0.1))
}

# Saving the grob to an RDS file for retrieval later on
saveRDS(ptx_plot,
	file = paste0("figures/grobs/", program, "_", id, "_signal.rds"),
	compress = FALSE)

