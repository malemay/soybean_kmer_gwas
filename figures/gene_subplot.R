# Loading the required libraries
suppressMessages(library(gwastools))
suppressMessages(library(grid))
suppressMessages(library(GenomicRanges))

# Getting the name of the trait-locus combination and of the program that called the genotypes from the command line
id <- commandArgs(trailingOnly = TRUE)[1]
program <- commandArgs(trailingOnly = TRUE)[2]

# Getting the ID of the gene associated with the trait and the name of the trait from the all_signals.rds file
# DEPENDENCY: utilities/all_signals.rds
target_signal <- readRDS("utilities/all_signals.rds")
target_signal <- target_signal[id]

gene_name <- target_signal$gene_name_v4
trait <- target_signal$trait

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

if(is.null(gene)) stop("No gene found where there should be one")

# Reading the data.frame of marker associations; treatment is different for k-mers than for other approaches
# DEPENDENCY: GWAS association results
# DEPENDENCY: GWAS thresholds

# Loading the results and signals from the rds file
gwas_results <- readRDS(paste0("gwas_results/", program, "/", id, "_locus_gwas.rds"))
gwas_signals <- readRDS(paste0("gwas_results/", program, "/", id, "_locus_signal.rds"))

# Getting the one signal associated with this gene
signal <- subsetByOverlaps(gwas_signals, gene)

if(!length(signal) == 1) warning("There is no signal overlapping the gene ", gene, " for locus ", id)

# Using the right viewport to focus on the p-values over the length of the gene
ptx_plot <- pvalue_tx_grob(gwas_results = gwas_results,
			       xscale = gene,
			       xexpand = c(0.1, 0.1),
			       yexpand = c(0.1, 0.1),
			       genes = genes,
			       transcripts = transcripts,
			       exons = exons,
			       cds = cds,
			       transcript_margins = c(0, 3.6, 0, 1.1),
			       pvalue_margins = c(5.1, 3.6, 0.5, 1.1))


# Saving the grob to an RDS file for retrieval later on
saveRDS(ptx_plot,
	file = paste0("figures/grobs/", program, "_", id, "_gene.rds"),
	compress = FALSE)

