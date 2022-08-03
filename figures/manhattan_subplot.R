# Loading the required packages
suppressMessages(library(grid))
suppressMessages(library(gwastools))
suppressMessages(library(GenomicRanges))


# Reading in the command line arguments; the first one is the trait, the second is the genotyping program
trait <- commandArgs(trailingOnly = TRUE)[1]
program <- commandArgs(trailingOnly = TRUE)[2]

# The location of the reference genome fasta
# DEPENDENCY: refgenome/Gmax_508_v4.0_mit_chlp.fasta
refgenome <- "refgenome/Gmax_508_v4.0_mit_chlp.fasta"

# A vector of Glycine max chromosome names
chromosomes <- paste0("Gm", ifelse(1:20 < 10, "0", ""), 1:20)

# It is a special case if we are analyzing k-mers
if(program == "kmers") {
	gwas_results <- readRDS(paste0("gwas_results/kmers/", trait, "_kmer_positions.rds"))

	gwas_plot <- manhattanGrob(gwas_results,
				   threshold = as.numeric(readLines(paste0("gwas_results/kmers/", trait, "_threshold_5per"))),
				   numeric_chrom = TRUE,
				   margins = c(5.1, 3.1, 0.1, 0.1))

} else if(program %in% c("platypus", "paragraph", "vg")) {
	# Loading the results from the CSV file and formatting them as a GRanges object
	gwas_results <- format_gapit_gwas(filename = paste0("gwas_results/", program, "/", trait, "_gwas.csv"),
					  ref_fasta = refgenome,
					  chromosomes = chromosomes,
					  pattern = "^Gm[0-9]{2}$")

	# Plotting the results using the gwastools::manhattan_plot function
	gwas_plot <- manhattanGrob(gwas_results,
				   threshold = -log10(as.numeric(readLines(paste0("gwas_results/", program, "/", trait, "_threshold_5per.txt")))),
				   numeric_chrom = TRUE,
				   margins = c(5.1, 3.1, 0.1, 0.1))
} else {
	stop("Unrecognized program option")
}

# Outputting to an rds file; actual plotting will be done in a different script
saveRDS(gwas_plot, file = paste0("figures/ggplots/", program, "_", trait, "_manhattan.rds"), compress = FALSE)

