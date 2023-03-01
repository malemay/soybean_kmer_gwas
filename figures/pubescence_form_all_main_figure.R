# Code for the figure on pubescence form (Pa1 locus) to include in the main text
# of the manuscript

# Loading the required libraries
suppressMessages(library(grid))
suppressMessages(library(gwask))
suppressMessages(library(GenomicRanges))

# Setting the trait and locus for this series of plots
locus <- commandArgs(trailingOnly = TRUE)[1]

# DEPENDENCY: utilities/kmer_plot_ranges.txt
locus_data <- read.table("utilities/kmer_plot_ranges.txt")
trait <- locus_data[locus_data[[1]] == locus, 2]

# Loading the grid objects to plot
# DEPENDENCY: figures/grobs/platypus_pubescence_form_all_manhattan.rds
# DEPENDENCY: figures/grobs/kmers_pubescence_form_all_manhattan.rds
platypus_gwide_manhattan <- readRDS(paste0("figures/grobs/platypus_", trait, "_manhattan.rds"))
kmers_gwide_manhattan <- readRDS(paste0("figures/grobs/kmers_", trait, "_manhattan.rds"))

# DEPENDENCY: figures/grobs/platypus_pubescence_form_all_Pa1_signal.rds
# DEPENDENCY: figures/grobs/kmers_pubescence_form_all_Pa1_signal.rds
platypus_zoomed_manhattan <- readRDS(paste0("figures/grobs/platypus_", locus, "_signal.rds"))
kmers_zoomed_manhattan <- readRDS(paste0("figures/grobs/kmers_", locus, "_signal.rds"))

# DEPENDENCY: figures/grobs/platypus_pubescence_form_all_Pa1_gene.rds
# DEPENDENCY: figures/grobs/kmers_pubescence_form_all_Pa1_gene.rds
platypus_gene_manhattan <- readRDS(paste0("figures/grobs/platypus_", locus, "_gene.rds"))
kmers_gene_manhattan <- readRDS(paste0("figures/grobs/kmers_", locus, "_gene.rds"))

# Loading the reference gene and transcript data for plotting the transcript
# DEPENDENCY: reference genome annotation
genes <- readRDS("refgenome/gmax_v4_genes.rds")
transcripts <- readRDS("refgenome/gmax_v4_transcripts.rds")
cds <- readRDS("refgenome/gmax_v4_cds.rds")
exons <- readRDS("refgenome/gmax_v4_exons.rds")

# Sourcing some functions used for plotting
# DEPENDENCY: figures/main_figure_functions.R
source("figures/main_figure_functions.R")

# Drawing the figure in a PNG device
png(paste0("figures/", locus, "_main_figure.png"), width = 6, height = 10, units = "in", res = 200)

# Resetting the plotting page
grid.newpage()

# Creating the layout for the whole figure
pushViewport(viewport(layout = grid.layout(nrow = 8, ncol = 1, heights = unit(c(3, 1, 3, 1.2, 1, 0.3, 3, 1), "null")),
		      name = "main"))

# Drawing the genome-wide Manhattan plots
pushViewport(viewport(layout.pos.row = 1, name = "gwide_manhattan"))
grid.text("(a)", x = 0.02, y = 0.95)
draw_manhattan(list(platypus_gwide_manhattan, kmers_gwide_manhattan), 
	       sigline_regexp = "sigline_[23]",
	       siglabel_regexp = "siglabel_[23]",
	       label_pos = 0.01,
	       labels = c("SNPs/indels", "k-mers"),
	       fontsize = 10)


# Drawing the zoomed-in Manhattan plots at the locus
seekViewport("main")
pushViewport(viewport(layout.pos.row = 3, name = "zoomed_manhattan"))
grid.text("(b)", x = 0.02, y = 0.95)
draw_zoomed(list(platypus_zoomed_manhattan, kmers_zoomed_manhattan),
	    label_pos = 0.01,
	    labels = c("SNPs/indels", "k-mers"),
	    fontsize = 10)

# Plotting the gene model along with p-values at this location
seekViewport("main")
pushViewport(viewport(layout.pos.row = 5, name = "transcript"))
grid.text("(c)", x = 0.02, y = 0.95)

stopifnot(identical(platypus_gene_manhattan$vp$xscale, kmers_gene_manhattan$vp$xscale))

# Push a viewport with the same margins as in the gene plots
pushViewport(plotViewport(c(0.5, 4.5, 0.5, 0.5),
			  xscale = platypus_gene_manhattan$vp$xscale))
grid.draw(transcriptsGrob(genes = genes,
			  transcripts = transcripts,
			  exons = exons,
			  cds = cds,
			  xscale = GRanges(seqnames = "Gm12", ranges = IRanges(start = platypus_gene_manhattan$vp$xscale[1],
									       end = platypus_gene_manhattan$vp$xscale[2])),
			  strand_colors = c("dodgerblue", "dodgerblue"),
			  draw_arrows = TRUE,
			  first_tx_only = TRUE))

# Drawing the zoomed-in Manhattan plots at the gene
seekViewport("main")
pushViewport(viewport(layout.pos.row = 7, name = "zoomed_manhattan"))
grid.text("(d)", x = 0.02, y = 0.95)
draw_zoomed(list(platypus_gene_manhattan, kmers_gene_manhattan), c("SNPs/indels", "k-mers"), fontsize = 10)

dev.off()

