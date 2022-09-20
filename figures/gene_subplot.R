# Loading the required libraries
suppressMessages(library(grid))
suppressMessages(library(gwastools))
suppressMessages(library(GenomicRanges))

# Setting some analysis parameters
xextend <- 0.1 # The expansion factor on either side of the gene for the transcript plot x-scale

# Getting the name of the trait-locus combination and of the program that called the genotypes from the command line
id <- commandArgs(trailingOnly = TRUE)[1]
program <- commandArgs(trailingOnly = TRUE)[2]

# Getting the ID of the gene associated with the trait and the name of the trait from the all_signals.rds file
# DEPENDENCY: utilities/all_signals.rds
target_signal <- readRDS("utilities/all_signals.rds")
target_signal <- target_signal[id]

gene_name <- target_signal$gene_name_v4
### ----------

# Getting the set of reference genes
# DEPENDENCY: refgenome/gmax_v4_genes.rds
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

if(is.null(gene)) stop("No gene found where there should be one")

# Extending the window for the causal gene
window_size <- GenomicRanges::width(gene)
start(gene) <- start(gene) - xextend * window_size
end(gene) <- end(gene) + xextend * window_size

# Reading the data.frame of marker associations; treatment is different for k-mers than for other approaches
# DEPENDENCY: GWAS association results
# DEPENDENCY: GWAS thresholds

# Loading the results and signals from the rds file
gwas_results <- readRDS(paste0("gwas_results/", program, "/", id, "_gwas_locus.rds"))
gwas_signals <- readRDS(paste0("gwas_results/", program, "/", id, "_signal_locus.rds"))

# Getting the name of the trait associated with the signal, and from it the significance threshold
trait <- target_signal$trait
# DEPENDENCY: significance thresholds (should not need to be included in Makefile, because the dependence is implicit)
threshold <- -log10(as.numeric(readLines(paste0("gwas_results/", program, "/", trait, "_threshold_5per.txt"))))

# Getting the one signal associated with this gene
signal <- subsetByOverlaps(gwas_signals, gene)

if(!length(signal) == 1) warning("There is no signal overlapping the gene ", gene, " for locus ", id)

# Using the right viewport to focus on the p-values over the length of the gene
ptx_plot <- pvalueGrob(gwas_results = gwas_results,
		       interval = gene,
		       feature = NULL,
		       threshold = threshold,
		       col = "blue",
		       pruned_col = if(program == "platypus") "firebrick2" else NULL,
		       yexpand = c(0.1, 0.1))

# Saving the grob to an RDS file for retrieval later on
saveRDS(ptx_plot,
	file = paste0("figures/grobs/", program, "_", id, "_gene.rds"),
	compress = FALSE)

