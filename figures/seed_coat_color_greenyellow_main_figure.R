# Code for the figure on seed coat color to include in the main text
# of the manuscript

# Loading the required libraries
suppressMessages(library(grid))
suppressMessages(library(gwask))
suppressMessages(library(GenomicRanges))

# Setting the trait and locus for this series of plots
locus <- commandArgs(trailingOnly = TRUE)[1]
cropping <- c(10, 15)

# DEPENDENCY: utilities/kmer_plot_ranges.txt
locus_data <- read.table("utilities/kmer_plot_ranges.txt")
trait <- locus_data[locus_data[[1]] == locus, 2]

# Loading the grid objects to plot
# DEPENDENCY: figures/grobs/platypus_seed_coat_color_greenyellow_manhattan.rds
# DEPENDENCY: figures/grobs/kmers_seed_coat_color_greenyellow_manhattan.rds
platypus_gwide_manhattan <- readRDS(paste0("figures/grobs/platypus_", trait, "_manhattan.rds"))
kmers_gwide_manhattan <- readRDS(paste0("figures/grobs/kmers_", trait, "_manhattan.rds"))

# DEPENDENCY: figures/grobs/platypus_seed_coat_color_greenyellow_Td_signal.rds
# DEPENDENCY: figures/grobs/kmers_seed_coat_color_greenyellow_Td_signal.rds
platypus_zoomed_manhattan <- readRDS(paste0("figures/grobs/platypus_", locus, "_signal.rds"))
kmers_zoomed_manhattan <- readRDS(paste0("figures/grobs/kmers_", locus, "_signal.rds"))

# Loading the data associated with the k-mer analysis
# DEPENDENCY: gwas_results/kmer_consensus/%_plotting_data.rds
# DEPENDENCY: gwas_results/kmer_consensus/%_difflist.rds
# DEPENDENCY: gwas_results/kmer_consensus/%_causal_gene.rds
# DEPENDENCY: gwas_results/kmer_consensus/%_phenodata.rds
plotting_data <- readRDS(paste0("gwas_results/kmer_consensus/", locus, "_plotting_data.rds"))
difflist <- readRDS(paste0("gwas_results/kmer_consensus/", locus, "_difflist.rds"))
causal_gene <- readRDS(paste0("gwas_results/kmer_consensus/", locus, "_causal_gene.rds"))
haplotype_data <- readRDS(paste0("gwas_results/kmers/", locus, "_phenodata.rds"))

# Loading the reference gene and transcript data for plotting the transcript
# DEPENDENCY: reference genome annotation
genes <- readRDS("refgenome/gmax_v4_genes.rds")
transcripts <- readRDS("refgenome/gmax_v4_transcripts.rds")
cds <- readRDS("refgenome/gmax_v4_cds.rds")
exons <- readRDS("refgenome/gmax_v4_exons.rds")

# Also extracting the plotting range
plotting_range <- locus_data[locus_data[[1]] == locus, 3]
chrom <- sub(":.*", "", plotting_range)
grange <- as.numeric(strsplit(sub(".*:", "", plotting_range), "-")[[1]])

# Sourcing some functions used for plotting
# DEPENDENCY: figures/main_figure_functions.R
source("figures/main_figure_functions.R")

# Drawing the figure in a PNG device
png(paste0("figures/", locus, "_main_figure.png"), width = 6, height = 10, units = "in", res = 200)

# Resetting the plotting page
grid.newpage()

# Creating the layout for the whole figure
pushViewport(viewport(layout = grid.layout(nrow = 9, ncol = 1, heights = unit(c(3, 1.2, 3, 1.5, 0.8, 0.3, 4, 1.5, 1.6), "null")),
		      name = "main"))

# Drawing the genome-wide Manhattan plots
pushViewport(viewport(layout.pos.row = 1, name = "gwide_manhattan"))
grid.text("(a)", x = 0.02, y = 0.95)
draw_manhattan(list(platypus_gwide_manhattan, kmers_gwide_manhattan), 
	       sigline_regexp = "sigline_1",
	       siglabel_regexp = "siglabel_1",
	       labels = c("SNP/indels", "k-mers"),
	       label_pos = 0.12,
	       fontsize = 10)


# Drawing the zoomed-in Manhattan plots
seekViewport("main")
pushViewport(viewport(layout.pos.row = 3, name = "zoomed_manhattan"))
grid.text("(b)", x = 0.02, y = 0.95)
draw_zoomed(list(platypus_zoomed_manhattan, kmers_zoomed_manhattan),
	    labels = c("SNP/indels", "k-mers"),
	    fontsize = 10)

# Plotting the gene model
seekViewport("main")
pushViewport(viewport(layout.pos.row = 5, name = "transcript"))
grid.text("(c)", x = 0.02, y = 0.95)
grid.draw(transcriptsGrob(genes = genes,
			  transcripts = transcripts,
			  exons = exons,
			  cds = cds,
			  xscale = causal_gene,
			  highlight = GenomicRanges::GRanges(seqnames = chrom,
							     ranges = IRanges::IRanges(start = grange[1] + cropping[1], end = grange[2] - cropping[2])),
			  strand_colors = c("dodgerblue", "dodgerblue"),
			  draw_arrows = TRUE,
			  first_tx_only = TRUE))

# Plotting the k-mer plot
seekViewport("main")
pushViewport(viewport(layout.pos.row = 7, name = "kmers"))
grid.text("(d)", x = 0.02, y = 0.95)
pushViewport(viewport(width = 0.90))
draw_haplotypes(plotting_data = plotting_data, difflist = difflist, plotting_range = plotting_range, cropping = cropping, fontsize = 8)


# Plotting the contingency table
seekViewport("main")
pushViewport(viewport(layout.pos.row = 9, name = "table"))
grid.text("(e)", x = 0.02, y = 0.95)
pushViewport(viewport(width = 0.85))
grid.phenotable(phenodata = haplotype_data)

dev.off()

