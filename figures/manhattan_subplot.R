# Loading the required packages
library(ggplot2)
library(gwastools)
library(GenomicRanges)


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

	gwas_plot <- manhattan_plot(formatted_data = gwas_results,
				    gwas_type = "kmer",
				    threshold = as.numeric(readLines(paste0("gwas_results/kmers/", trait, "_threshold_5per")))) +
	ggplot2::ggtitle(paste0("kmer-", trait)) +
	theme(text = element_text(size = 8))

} else if(program %in% c("platypus", "paragraph", "vg")) {
	# Loading the results from the CSV file and formatting them as a GRanges object
	gwas_results <- format_gapit_gwas(filename = paste0("gwas_results/", program, "/", trait, "_gwas.csv"),
					  ref_fasta = refgenome,
					  chromosomes = chromosomes,
					  pattern = "^Gm[0-9]{2}$")

	# Plotting the results using the gwastools::manhattan_plot function
	gwas_plot <- manhattan_plot(gwas_results,
				    gwas_type = "gapit",
				    threshold = -log10(as.numeric(readLines(paste0("gwas_results/", program, "/",
										   trait, "_threshold_5per.txt"))))) +
	ggplot2::ggtitle(paste0(program, "-", trait)) +
	theme(text = element_text(size = 8))
} else {
	stop("Unrecognized program option")
}

# Outputting to an rds file; actual plotting will be done in a different script
saveRDS(gwas_plot, file = paste0("figures/ggplots/", program, "_", trait, "_manhattan.rds"), compress = FALSE)

