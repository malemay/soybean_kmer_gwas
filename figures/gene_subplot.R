# Loading the required libraries
suppressMessages(library(grid))
suppressMessages(library(gwask))
suppressMessages(library(GenomicRanges))

# Getting the name of the trait-locus combination and of the program that called the genotypes from the command line
id <- commandArgs(trailingOnly = TRUE)[1]
program <- commandArgs(trailingOnly = TRUE)[2]

# Setting some analysis parameters
xextend <- 0.1 # The expansion factor on either side of the gene for the transcript plot x-scale
# The R locus is a special case because I want to see its promoter region
# if(id == "hilum_color_blackbrown_R") xextend <- 2

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

# Reading the CNV range (if applicable)
# DEPENDENCY: utilities/cnv_genes.txt
cnv_genes <- read.table("utilities/cnv_genes.txt")
if(id %in% cnv_genes[[1]]) {
	gene <- cnv_range <- readRDS(cnv_genes[cnv_genes[[1]] == id, 2])
} else {
	cnv_range <- NULL
}

# The Ps locus is a special case because we want to set the plotting window to include the CNV
if(id == "pubescence_density_Ps") gene <- cnv_range

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

# Setting the feature optionally (just for genes for which we want to highlight specific variants)
# DEPENDENCY: utilities/variant_ranges.txt
variant_ranges <- read.table("utilities/variant_ranges.txt", head = FALSE, sep = "\t", stringsAsFactors = FALSE)
key <- paste0(program, "_", id)

if(length(index <- which(key == variant_ranges[[1]]))) {
	feature <- as(variant_ranges[index, 2], "GRanges")
} else {
	feature <- NULL
}

# Using the right viewport to focus on the p-values over the length of the gene
ptx_plot <- pvalueGrob(gwas_results = gwas_results,
		       interval = gene,
		       shading = cnv_range,
		       feature = feature,
		       threshold = threshold,
		       col = "blue",
		       pruned_col = if(program == "platypus") "firebrick2" else NULL,
		       cex.points = 0.7,
		       yexpand = c(0.1, 0.1))

# Saving the grob to an RDS file for retrieval later on
saveRDS(ptx_plot,
	file = paste0("figures/grobs/", program, "_", id, "_gene.rds"),
	compress = FALSE)

