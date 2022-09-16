# Loading the required packages
suppressMessages(library(grid))
suppressMessages(library(gwastools))
suppressMessages(library(GenomeInfoDb))
suppressMessages(library(GenomicRanges))

# Reading in the command line arguments; the first one is the trait, the second is the genotyping program
trait <- commandArgs(trailingOnly = TRUE)[1]
program <- commandArgs(trailingOnly = TRUE)[2]

# Reading the threshold from file
# DEPENDENCY: significance threshold
threshold <- -log10(as.numeric(readLines(paste0("gwas_results/", program, "/", trait, "_threshold_5per.txt"))))

# Loading the results from the rds file
# DEPENDENCY: GWAS results
gwas_results <- readRDS(paste0("gwas_results/", program, "/", trait, "_gwas.rds"))

# We filter out variants other than unachored scaffolds
GenomeInfoDb::seqlevels(gwas_results, pruning.mode = "coarse") <- grep("^Gm[0-9]{2}$", GenomeInfoDb::seqlevels(gwas_results), value = TRUE, invert = TRUE)
# And we keep only the reference sequences that are present
GenomeInfoDb::seqlevels(gwas_results, pruning.mode = "coarse") <- seqlevels(gwas_results)[seqlevels(gwas_results) %in% seqnames(gwas_results)]

# We first check if there are any markers left
if(!length(gwas_results)) {
	gwas_plot <- textGrob(paste0("No significant associations found for ", trait, " using ", program))
} else {
	# Plotting the results using the gwastools::manhattan_plot function
	gwas_plot <- manhattanGrob(gwas_results,
				   threshold = threshold,
				   ref_signals = NULL,
				   new_signals = NULL,
				   numeric_chrom = FALSE,
				   margins = c(5.1, 4.1, 0.1, 0.1))
}

# Plotting directly as we are only using a single program (the k-mers)
# Outputting to a PNG file
png(paste0("figures/", trait, "_scaffolds_manhattan.png"), width = 9, height = 3, units = "in", res = 100)
grid.newpage()
grid.draw(gwas_plot)
dev.off()

