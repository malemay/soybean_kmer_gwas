# This script provides functionality for computing and plotting
# the linkage disequilibrium (LD) for a set of k-mers using
# their presence/absence pattern.

# Loading the required libraries
suppressMessages(library(GenomicRanges))
suppressMessages(library(gwask))

# Setting the trait to analyze and loading the k-mer results for
trait <- commandArgs(trailingOnly = TRUE)[1]

# Creating variables for the paths to external programs
# WARNING: should consider setting these paths in the Makefile
filter_kmers <- "~/programs/kmers_gwas/bin/filter_kmers "
kmers_ld <- "~/scripts/kmers_ld/kmers_ld "

# Creating variables for the file paths to improve readability
# DEPENDENCY: kmers_table/kmers_table.table
# DEPENDENCY: kmers_table/kmers_table.names
kmers_table <- "kmers_table/kmers_table"

kmers_file <- paste0("gwas_results/kmers/", trait, "_significant_kmers.txt")
pav_file <- paste0("gwas_results/kmers/", trait, "_pav_table.txt")
nodup_pav_file <- paste0("gwas_results/kmers/", trait, "_pav_table_nodup.txt")
clustered_ld_file <- paste0("gwas_results/kmers/", trait, "_clustered_ld.txt")

# Reading the aligned significant k-mers for this trait
# DEPENDENCY: positions of the significant k-mers for trait
kmer_positions <- readRDS(paste0("gwas_results/kmers/", trait, "_gwas.rds"))
names(kmer_positions) <- kmer_positions$kmer_canon

# Optionally doing some processing to subsample only some k-mers if there are too many
if(length(kmer_positions) > 1500) {
	kmer_positions <- subsample_kmers(kmer_positions, nkmers = 1500, npvalue = 500)
}

# Writing the significant k-mers to file
significant_kmers <- kmer_positions$kmer_canon
writeLines(significant_kmers, con = kmers_file)

# Extracting the presence/absence values for those k-mers using a k-mers GWAS program
command <- paste0(filter_kmers, " -t ", kmers_table, " -k ", kmers_file, " -o ", pav_file)
system(command)

# Reading the PAV table from file and removing variants that have exactly the same PAV profile (removed this: now all variants are considered)
pav_table <- read.table(pav_file, header = TRUE)
#pav_table <- pav_table[!duplicated(pav_table[, -1]), ]

# Writing this de-duplicated PAV table to file
write.table(pav_table, file = nodup_pav_file, sep = "\t",
	    col.names = TRUE, row.names = FALSE, quote = FALSE)

# Computing the LD matrix using a C program
system(paste0(kmers_ld, nodup_pav_file, " > ", clustered_ld_file))

