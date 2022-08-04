# Loading the required libraries
suppressMessages(library(grid))
suppressMessages(library(gwastools))
suppressMessages(library(GenomicRanges))

# Getting the name of the trait-locus combination and of the program that called the genotypes from the command line
id <- commandArgs(trailingOnly = TRUE)[1]
program <- commandArgs(trailingOnly = TRUE)[2]

# A parameter that sets the padding of the x-scale when the signal consists of only one marker
xscale_padding <- 500
# A parameter that sets the distance between markers for merging signals with extract_signals
signal_distance <- 250000

# Getting the ID of the gene associated with the locus
# DEPENDENCY: utilities/all_signals.rds
target_signal <- readRDS("utilities/all_signals.rds")
target_signal <- target_signal[id]

stopifnot(length(target_signal) == 1)

gene_name <- target_signal$gene_name_v4

# A vector of Glycine max chromosome names
chromosomes <- paste0("Gm", ifelse(1:20 < 10, "0", ""), 1:20)
### ----------

# Getting the set of genes, transcripts, exons and coding sequences
# DEPENDENCY: refgenome/*.rds
genes <- readRDS("refgenome/gmax_v4_genes.rds")
transcripts <- readRDS("refgenome/gmax_v4_transcripts.rds")
exons <- readRDS("refgenome/gmax_v4_exons.rds")
cds <- readRDS("refgenome/gmax_v4_cds.rds")

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

# Reading the data.frame of marker associations; treatment is different for k-mers than for other approaches
# DEPENDENCY: GWAS association results
# DEPENDENCY: GWAS thresholds
threshold <- -log10(as.numeric(readLines(paste0("gwas_results/", program, "/", id, "_locus_threshold_5per.txt"))))

# Reading the GWAS results from the RDS file
gwas_results <- readRDS(paste0("gwas_results/", program, "/", id, "_locus_gwas.rds"))

# Extracting the signal at the location of the gene
gwas_signals <- extract_signals(gwas_results, threshold = threshold, distance = signal_distance)
signal <- subsetByOverlaps(gwas_signals, target_signal)

if(!length(signal)) {
	warning("No signal found by ", program, " for signal ID ", id)
	ptx_plot <- grid::gTree(children = gList(textGrob(paste0("No signal found by ", program, " for signal ID ", id))))

} else if (length(signal) > 1) {
	warning(length(signal), " signals found by ", program, " for signal ID ", id)
	ptx_plot <- grid::gTree(children = gList(textGrob(paste0(length(signal), " signals found by ", program, " for signal ID ", id))))

} else {
	# Setting the x-scale while handling special cases where there is no signal or the signal is only 1 bp wide
	if(width(signal) == 1) {
		start(signal) <- start(signal) - xscale_padding
		end(signal)   <- end(signal) + xscale_padding
	}

	# Creating a grob representing the p-value plot and transcript
	ptx_plot <- pvalue_tx_grob(gwas_results = gwas_results,
				   feature = gene,
				   xscale = signal,
				   first_tx_only = TRUE,
				   xexpand = c(0.05, 0.05),
				   yexpand = c(0.1, 0.1),
				   genes = genes,
				   transcripts = transcripts,
				   exons = exons,
				   cds = cds,
				   transcript_margins = c(0, 3.6, 0, 1.1),
				   pvalue_margins = c(5.1, 3.6, 0.5, 1.1),
				   pvalue_fraction = 0.95)
}

# Saving the grob to an RDS file for retrieval later on
saveRDS(ptx_plot,
	file = paste0("figures/grobs/", program, "_", id, "_signal.rds"),
	compress = FALSE)

