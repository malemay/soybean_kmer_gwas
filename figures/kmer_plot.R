# This script generates a k-mer plot from the consensus sequences of a group of samples at a given locus

# Loading the required libraries
library(grid)
library(Biostrings)
library(IRanges)
library(stringr)
library(gwastools)

# Setting some analysis parameters
max_kmers <- 10000 # the maximum number of k-mers to map onto the haplotype sequences
kmer_length <- 31 # the length of the k-mers used in the analysis
min_sequences <- 80 # min number of samples that must have a consensus sequence from assembly; otherwise obtained from bwa
min_frequency <- 5 # min number of times that a haplotype must occur to be used for plotting

# Reading the locus from the command line
locus <- commandArgs(trailingOnly = TRUE)[1]

# Getting the trait that corresponds to that locus
# DEPENDENCY: utilities/kmer_plot_ranges.txt
locus_data <- read.table("utilities/kmer_plot_ranges.txt")
trait <- locus_data[locus_data[[1]] == locus, 2]

# Getting the lookup table that links the trait name used in the USDA database to the one used for our analyses
# DEPENDENCY: phenotypic_data/trait_names.rds
trait_names <- readRDS("phenotypic_data/trait_names.rds")
usda_trait <- names(trait_names)[sapply(trait_names, function(x) grepl(paste0("^", x), trait))]

# Getting the phenotypic data
# DEPENDENCY: phenotypic_data/phenotypic_data.csv
phenotypes <- read.table("phenotypic_data/phenotypic_data.csv", sep = ";", header = TRUE,)

# Reading the k-mers associated with that analysis and their p-values
# DEPENDENCY: k-mer p-values
kmer_pvalues <- read_kmer_pvalues(kmer_file = paste0("gwas_results/kmer_data/", trait,
						     "/katcher_results/pass_threshold_5per_sorted.txt"),
				  max_kmers = max_kmers,
				  kmer_length = 31)

# Reading the consensus sequence from the assemblies
# DEPENDENCY: consensus sequence computed from the assembly of significant reads
sequences <- read_consensus(input_fasta = paste0("gwas_results/kmer_consensus/", locus, "_sequences.fa"))

# DEPENDENCY: consensus sequence computed from the alignment using BWA
# Using consensus sequences from bwa in case the assemblies did not yield great results
if(length(sequences) < min_sequences) {
	warning("Only ", length(sequences), "have a consensus sequence from assembly.",
	       	"Using consensus sequences from read alignment by bwa for locus ", locus)
	sequences <- read_consensus(paste0("gwas_results/kmer_consensus/", locus, "_bwa_sequences.fa"))
}

# Get a character vector of the haplotypes found with min_frequency in the dataset
haplotypes <- get_haplotypes(sequences = sequences, min_frequency = min_frequency)

# Creating a data.frame linking haplotypes to phenotypes
haplotype_data <- link_phenotypes(sequences = sequences,
				  haplotypes = haplotypes,
				  phenotypes = phenotypes,
				  id_column = "bayer_id",
				  phenotype_column = usda_trait)

# Get the positions in the haplotypes that overlap any of the significant k-mers
kmer_overlaps <- match_kmers(haplotypes = haplotypes,
			     kmers = kmer_pvalues,
			     kmer_length = kmer_length)


# Generating a data.frame with separate nucleotides and their associated p-value for each haplotype
plotting_data <- format_haplotypes(haplotypes = haplotypes,
				   overlaps = kmer_overlaps)


# Performing multiple alignment of the sequences to find the gaps
alignment <- mafft_align(fasta_path = paste0("gwas_results/kmer_consensus/", locus, "_haplotypes.fa"),
			 haplotypes = haplotypes,
			 mafft_path = "mafft",
			 mafft_options = "--auto")

# Adjusting the positions for plotting depending on the gaps in the alignment
plotting_data <- adjust_gaps(hapdata = plotting_data,
			     alignment = alignment)

# Finding the positions where the aligned nucleotides differ between haplotypes
difflist <- nucdiff(plotting_data)

# Determining the plotting color from a common palette for all haplotypes
max_pvalue <- max(do.call("rbind", plotting_data)$log10p)

plotting_data <- lapply(plotting_data, function(x) {
		       x$color <- map_color(values = x$log10p, min = 0, max = max_pvalue, n.colors = 9, pal = "YlOrRd")
		       x
		       })


# Outputting to a png file
png(paste0("figures/", locus, "_kmers.png"), width = 8, height = 2, units = "in", res = 100)

# Initializing the device
grid.newpage()

# Drawing a box around the plotting region
grid.rect()

# Dividing the viewport into three viewports:
# - The first viewport will be used for showing the haplotype sequences
# - The second viewport will be used for the colour scale of the p-values
# - The third viewport will be used for showing a table of the phenotypes observed per haplotype
pushViewport(viewport(layout = grid.layout(nrow = 3, heights = unit(c(0.3, 0.1, 0.6), "npc"))))


# Moving into the viewport associated with the sequences and drawing them
pushViewport(viewport(layout.pos.row = 1))
grid.haplotypes(plotting_data, difflist, fontsize = 8)
upViewport()

# Moving into the viewport for the color scale
pushViewport(viewport(layout.pos.row = 2))
grid.text("Color scale")
upViewport()

# Moving into the viewport for the table
pushViewport(viewport(layout.pos.row = 3))
grid.phenotable(haplotype_data)
upViewport()

dev.off()

