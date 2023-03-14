# Loading the required packages
suppressMessages(library(gwask))
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

} else if(program %in% c("platypus", "paragraph")) {
	# Loading the results from the CSV file and formatting them as a GRanges object
	gwas_results <- format_gapit_gwas(filename = paste0("gwas_results/", program, "/", trait, "_gwas.csv"),
					  ref_fasta = refgenome,
					  chromosomes = chromosomes,
					  vcf_file = paste0("filtered_variants/", program, "/filtered_variants.vcf.gz"),
					  pattern = "^Gm[0-9]{2}$")
} else {
	stop("Unrecognized program option")
}

# Outputting to an rds file for retrieval in downstream analyses
saveRDS(gwas_results, file = paste0("gwas_results/", program, "/", trait, "_gwas.rds"), compress = FALSE)

