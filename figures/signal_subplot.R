# Loading the required libraries
library(gwastools)
library(grid)
library(GenomicRanges)

# Getting the name of the trait-locus combination and of the program that called the genotypes from the command line
id <- commandArgs(trailingOnly = TRUE)[1]
program <- commandArgs(trailingOnly = TRUE)[2]

# Getting the ID of the gene associatde with the trait and the name of the trait from the utilities/loci.txt file
# DEPENDENCY: utilities/loci.txt
loci <- read.table("utilities/loci.txt", header = TRUE, sep = ",", stringsAsFactors = FALSE)
gene_v1 <- loci[loci$id == id, "gene_v1"]
trait <- loci[loci$id == id, "trait"]

# The reference genome is needed for formatting the GWAS results from programs other than k-mers
# DEPENDENCY: refgenome/Gmax_508_v4.0_mit_chlp.fasta
refgenome <- "refgenome/Gmax_508_v4.0_mit_chlp.fasta"

# A vector of Glycine max chromosome names
chromosomes <- paste0("Gm", ifelse(1:20 < 10, "0", ""), 1:20)
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

# Reading the data.frame of marker associations; treatment is different for k-mers than for other approaches
# DEPENDENCY: GWAS association results
# DEPENDENCY: GWAS thresholds
if(program == "kmers") {
	gwas_results <- readRDS(paste0("gwas_results/kmers/", trait, "_kmer_positions.rds"))
	threshold <- as.numeric(readLines(paste0("gwas_results/kmers/", trait, "_threshold_5per")))
} else if(program %in% c("platypus", "paragraph", "vg")) {
	# Loading the results from the CSV file and formatting them as a GRanges object
	gwas_results <- format_gapit_gwas(filename = paste0("gwas_results/", program, "/", trait, "_gwas.csv"),
					  ref_fasta = refgenome,
					  chromosomes = chromosomes,
					  pattern = "^Gm[0-9]{2}$")
	threshold <- -log10(as.numeric(readLines(paste0("gwas_results/", program, "/", trait, "_threshold_5per.txt"))))
}

# Extracting the signal at the location of the gene
gwas_signals <- extract_signals(gwas_results, threshold = threshold, distance = 10^5)

# The next steps depend on whether there is a known gene model associated with the phenotype
if(!is.na(gene_v1) && !is.na(gene_v4 <- lookup_v1_to_v4[gene_v1])) {

	# Getting the one signal associated with this gene
	signal <- subsetByOverlaps(gwas_signals, genes[gene_v4])

	if(!length(signal) == 1) warning("There is no signal overlapping the gene ", gene_v4, " for phenotype ", trait)

	# Setting the x-scale while handling special cases where there is no signal or the signal is only 1 bp wide
	if(length(signal)) {
		if(width(signal) > 1) {
			xscale <- signal
		} else {
			start(signal) <- start(signal) - 500
			end(signal)   <- end(signal) + 500
			xscale <- signal
		}
	} else {
		xscale <- genes[gene_v4]
	}

	# Creating a grob representing the p-value plot and transcript
	ptx_plot <- pvalue_tx_grob(gwas_results = gwas_results,
				   feature = genes[gene_v4],
				   xscale = xscale,
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

} else {
	# Otherwise we get the signal with the greatest number of markers
	signal <- gwas_signals[which.max(gwas_signals$n_markers)]

	if(length(signal)) {
		if(width(signal) > 1) {
			xscale <- signal
		} else {
			start(signal) <- start(signal) - 500
			end(signal)   <- end(signal) + 500
			xscale <- signal
		}

		ptx_plot <- pvalue_tx_grob(gwas_results = gwas_results,
					   xscale = xscale,
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

	} else {
		warning("No signal found")
		ptx_plot <- gTree()
	} 
}

# Saving the grob to an RDS file for retrieval later on
saveRDS(ptx_plot,
	file = paste0("figures/ggplots/", program, "_", id, "_signal.rds"),
	compress = FALSE)

