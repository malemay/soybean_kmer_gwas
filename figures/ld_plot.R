# This scripts plots the pairwise LD between a set of k-mers based
# on a pre-computed LD matrix

# Loading the required packages
suppressMessages(library(gwastools))
suppressMessages(library(grid))
suppressMessages(library(GenomicRanges))

# Getting the trait to analyze from the command line
trait <- commandArgs(trailingOnly = TRUE)[1]

# Reading the LD matrix from file
# DEPENDENCY: clustered LD matrix for trait considered
clustered_ld_file <- paste0("gwas_results/kmers/", trait, "_clustered_ld.txt")
clustered_ld <- as.matrix(read.table(clustered_ld_file, header = FALSE, row.names = 1, sep = "\t"))
colnames(clustered_ld) <- rownames(clustered_ld)

# Reading the aligned significant k-mers for this trait
# DEPENDENCY: positions of the significant k-mers for trait
kmer_positions <- readRDS(paste0("gwas_results/kmers/", trait, "_gwas.rds"))
names(kmer_positions) <- kmer_positions$kmer_canon

# Sorting the k-mers in the LD matrix based on their genomic position
gsorted_ld <- ld_sort(clustered_ld, kmer_positions, sort_param = "position")

# Plotting the results using markers sorted by genomic position
png(paste0("figures/", trait, "_ld.png"), width = 9, height = 10, units = "in", res = 400)

grid.newpage()
ld_plot(gsorted_ld, kmer_positions, top_legend = FALSE, ylabels = TRUE, ylabel_pattern = "^Gm[0-9]{2}$")

dev.off()


# Getting a GRanges objct describing the whole genome
# DEPENDENCY: refgenome/Gmax_508_v4.0_mit_chlp.fasta.fai
#refgenome <- Rsamtools::scanFaIndex("refgenome/Gmax_508_v4.0_mit_chlp.fasta")
#names(refgenome) <- as.character(seqnames(refgenome))

# Plotting the results using a subset of markers located at certain genomic positions
#pubescence_form_all_subset <- GRanges(seqnames ="Gm12", ranges = IRanges(start = 38700000, end = 38800000))
#pubescence_form_all_subset <- c(pubescence_form_all_subset, refgenome[c("Gm04", "Gm13", "Gm15")])

#pubescence_color_all_subset <- GRanges(seqnames = "Gm06", ranges = IRanges(start = 18100000, end = 18800000))
#pubescence_color_all_subset <- c(pubescence_color_all_subset, refgenome[c("Gm06_scaffold_183", "Gm16", "Gm01")])

#grid.newpage()
#ld_plot(ld_sort(clustered_ld, kmer_positions, sort_param = "position", positions = pubescence_color_all_subset),
	#kmer_positions)

