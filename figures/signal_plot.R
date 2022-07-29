# Loading the required libraries
library(gwastools)
library(grid)
library(GenomicRanges)

### ---------- ANALYSIS PARAMETERS
# Getting the name of the trait-locus combination from the command line
id <- commandArgs(trailingOnly = TRUE)[1]

loci <- read.table("utilities/loci.txt", header = TRUE, sep = ",", stringsAsFactors = FALSE)
gene_v1 <- loci[loci$id == id, "gene_v1"]
trait <- loci[loci$id == id, "trait"]
### ----------

# Loading the lookup table from Gmax v1 gene names to Gmax v4 gene names
# DEPENDENCY: refgenome/lookup_v1_to_v4.rds
lookup_v1_to_v4 <- readRDS("refgenome/lookup_v1_to_v4.rds")

# Getting the set of genes, transcripts, exons and coding sequences
# DEPENDENCY: refgenome/*.rds
genes <- readRDS("refgenome/gmax_v4_genes.rds")
transcripts <- readRDS("refgenome/gmax_v4_transcripts.rds")
exons <- readRDS("refgenome/gmax_v4_exons.rds")
cds <- readRDS("refgenome/gmax_v4_cds.rds")

# Reading the data.frame of significantly associated k-mers and formatting it
# DEPENDENCY : refgenome/Gmax_508_v4.0_mit_chlp.fasta.fai
# DEPENDENCY: kmer positions
kmer_results <- readRDS(paste0("gwas/kmers/", trait, "/katcher_results/", trait, "_kmer_positions.rds"))

# Extracting the signal at the location of the gene
kmer_signals <- extract_signals(kmer_results,
				threshold = as.numeric(readLines(paste0("gwas/kmers/", trait, "/kmers/threshold_5per"))),
				distance = 10^4)

# The next steps depend on whether there is a known gene model associated with the phenotype
if(!is.na(gene_v1)) {

	# Looking up the name of the gene in version 4
	gene_v4 <- lookup_v1_to_v4[gene_v1]

	# Getting the one signal associated with this gene
	signal <- subsetByOverlaps(kmer_signals, genes[gene_v4])

	if(!length(signal) == 1) {
		stop("There is no signal overlapping the gene ", gene_v4, " for phenotype ", trait)
	}

	# Plotting using grid functions and the gwastools::pvalue_tx_plot function
	png(paste0("figures/", id, "_signal.png"), width = 8, height = 10, units = "in", res = 300)

	grid.newpage()

	# Creating a horizontal viewport layout
	pushViewport(viewport(layout = grid.layout(ncol = 3,
						   widths = unit(c(0.48, 1, 0.48), c("npc", "null", "npc")))))

	# Using the left viewport to plot the p-values over the whole GWAS signal
	pushViewport(viewport(layout.pos.col = 1))
	pvalue_tx_plot(gwas_results = kmer_results,
		       feature = genes[gene_v4],
		       xscale = signal,
		       xexpand = c(0.05, 0.05),
		       yexpand = c(0.1, 0.1),
		       genes = genes,
		       transcripts = transcripts,
		       exons = exons,
		       cds = cds,
		       transcript_margins = c(0, 3.6, 0, 1.1),
		       pvalue_margins = c(5.1, 3.6, 0.5, 1.1))
	upViewport()

	# Using the right viewport to focus on the p-values over the length of the gene
	pushViewport(viewport(layout.pos.col = 3))
	pvalue_tx_plot(gwas_results = kmer_results,
		       xscale = genes[gene_v4],
		       xexpand = c(0.1, 0.1),
		       yexpand = c(0.1, 0.1),
		       genes = genes,
		       transcripts = transcripts,
		       exons = exons,
		       cds = cds,
		       transcript_margins = c(0, 3.6, 0, 1.1),
		       pvalue_margins = c(5.1, 3.6, 0.5, 1.1))
	upViewport()
	dev.off()
} else {
	# Otherwise we get the signal with the greatest number of markers
	signal <- kmer_signals[which.max(kmer_signals$n_markers)]

	png(paste0("figures/", id, "_signal.png"), width = 8, height = 10, units = "in", res = 300)

	pvalue_tx_plot(gwas_results = kmer_results,
		       xscale = signal,
		       xexpand = c(0.05, 0.05),
		       yexpand = c(0.1, 0.1),
		       genes = genes,
		       transcripts = transcripts,
		       exons = exons,
		       cds = cds,
		       transcript_margins = c(0, 3.6, 0, 1.1),
		       pvalue_margins = c(5.1, 3.6, 0.5, 1.1))

	dev.off()
}

