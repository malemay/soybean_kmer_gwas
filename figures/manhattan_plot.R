# Loading the required packages
library(grid)
library(ggplot2)
library(gwastools)
library(GenomicRanges)

# The location of the reference genome fasta
# DEPENDENCY: refgenome/Gmax_508_v4.0_mit_chlp.fasta
refgenome <- "refgenome/Gmax_508_v4.0_mit_chlp.fasta"

# A vector of Glycine max chromosome names
chromosomes <- paste0("Gm", ifelse(1:20 < 10, "0", ""), 1:20)

i <- commandArgs(trailingOnly = TRUE)[1]

#-------------------------------------------------- 
message("Reading in the GWAS results from vg on ", i)

# Loading the results from the CSV file and formatting them as a GRanges object
vg_results <- format_gapit_gwas(filename = paste0("gwas/vg_gwas/", i, "/GAPIT.MLM.", i, ".GWAS.Results.csv"),
				ref_fasta = refgenome,
				chromosomes = chromosomes,
				pattern = "^Gm[0-9]{2}$")

# Plotting the results using the gwastools::manhattan_plot function
vg_plot <- manhattan_plot(vg_results,
			  gwas_type = "gapit",
			  threshold = -log10(as.numeric(readLines(paste0("gwas/vg_gwas/",
									 i, "/", i,
									 "_threshold_5per.txt"))))) +
	ggplot2::ggtitle(paste0("vg-", i)) +
	theme(text = element_text(size = 8))

#-------------------------------------------------- 
message("Reading in the GWAS results from platypus on ", i)

# Loading the results from the CSV file and formatting them as a GRanges object
platypus_results <- format_gapit_gwas(filename = paste0("gwas/platypus/", i, "/GAPIT.MLM.", i, ".GWAS.Results.csv"),
				      ref_fasta = refgenome,
				      chromosomes = chromosomes,
				      pattern = "^Gm[0-9]{2}$")

# Plotting the results using the gwastools::manhattan_plot function
platypus_plot <- manhattan_plot(platypus_results,
				gwas_type = "gapit",
				threshold = -log10(as.numeric(readLines(paste0("gwas/platypus/",
									       i, "/", i,
									       "_threshold_5per.txt"))))) +
	ggplot2::ggtitle(paste0("platypus-", i)) +
	theme(text = element_text(size = 8))

# -------------------------------------------------- 
message("Reading in the GWAS results from Paragraph on ", i)

# Loading the results from the CSV file and formatting them as a GRanges object
paragraph_results <- format_gapit_gwas(filename = paste0("gwas/paragraph/", i, "/GAPIT.MLM.", i, ".GWAS.Results.csv"),
				       ref_fasta = refgenome,
				       chromosomes = chromosomes,
				       pattern = "^Gm[0-9]{2}$")

#Plotting the results using the gwastools::manhattan_plot function
paragraph_plot <- manhattan_plot(paragraph_results,
				 gwas_type = "gapit",
				 threshold = -log10(as.numeric(readLines(paste0("gwas/paragraph/",
										i, "/", i,
										"_threshold_5per.txt"))))) +
	ggplot2::ggtitle(paste0("paragraph-", i)) +
	theme(text = element_text(size = 8))

#-------------------------------------------------- 
message("Reading in the k-mer GWAS results on ", i)

kmer_results <- readRDS(paste0("gwas/kmers/", i, "/katcher_results/", i, "_kmer_positions.rds"))

kmer_plot <- manhattan_plot(formatted_data = kmer_results,
			    gwas_type = "kmer",
			    threshold = as.numeric(readLines(paste0("gwas/kmers/", i, "/kmers/threshold_5per")))) +
	ggplot2::ggtitle(paste0("kmer-", i)) +
	theme(text = element_text(size = 8))

# Outputting to a PNG file
png(paste0("figures/", i, "_manhattan.png"), width = 9, height = 10, units = "in", res = 400)
# Resetting the viewport
grid.newpage()

# Setting a viewport with a top panel for vg results and a bottom panel for Paragrpah results
pushViewport(viewport(layout = grid.layout(nrow = 4)))
print(platypus_plot, vp = viewport(layout.pos.row = 1))
print(vg_plot, vp = viewport(layout.pos.row = 2))
print(paragraph_plot, vp = viewport(layout.pos.row = 3))
print(kmer_plot, vp = viewport(layout.pos.row = 4))

dev.off()

