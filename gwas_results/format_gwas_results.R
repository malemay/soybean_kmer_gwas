# Loading the required packages
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
	# For k-mers we simply copy the kmer_positions file to a new file
	gwas_results <- readRDS(paste0("gwas_results/kmers/", trait, "_kmer_positions.rds"))

} else if(program %in% c("platypus", "paragraph", "vg")) {
	# Loading the results from the CSV file and formatting them as a GRanges object
	gwas_results <- format_gapit_gwas(filename = paste0("gwas_results/", program, "/", trait, "_gwas.csv"),
					  ref_fasta = refgenome,
					  chromosomes = chromosomes,
					  pattern = "^Gm[0-9]{2}$")
} else {
	stop("Unrecognized program option")
}

# Outputting to an rds file for retrieval in downstream analyses
saveRDS(gwas_results, file = paste0("gwas_results/", program, "/", trait, "_gwas.rds"), compress = FALSE)

# Also creating the appropriate symlinks such that signals can relate to their underlying results file
signal_ids <- read.table("utilities/signal_ids.txt", header = FALSE, sep = ",")
signal_ids <- signal_ids[signal_ids[[2]] == trait, ]

if(nrow(signal_ids)) {

	for(i in 1:nrow(signal_ids)) {
		gwas_file1 <- paste0(signal_ids[i, 2], "_gwas.rds")
		gwas_file2 <- paste0("gwas_results/", program, "/", signal_ids[i, 1], "_locus_gwas.rds")
		if(file.exists(gwas_file2)) unlink(gwas_file2)
		file.symlink(gwas_file1, gwas_file2)

		threshold_file1 <- paste0(signal_ids[i, 2], "_threshold_5per.txt")
		threshold_file2 <- paste0("gwas_results/", program, "/", signal_ids[i, 1], "_locus_threshold_5per.txt")

		if(file.exists(threshold_file2)) unlink(threshold_file2)
		file.symlink(threshold_file1, threshold_file2)
	}

}
