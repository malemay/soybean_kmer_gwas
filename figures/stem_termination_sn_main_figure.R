# Code for the figure on pubescence color (Td locus) to include in the main text
# of the manuscript

# Loading the required libraries
suppressMessages(library(grid))
suppressMessages(library(gwask))
suppressMessages(library(GenomicRanges))

# DEPENDENCY: utilities/kmer_plot_ranges.txt
trait <- "stem_termination_sn"

# Loading the grid objects to plot
# DEPENDENCY: figures/grobs/kmers_stem_termination_sn_manhattan.rds
kmers_gwide_manhattan <- readRDS(paste0("figures/grobs/kmers_", trait, "_manhattan.rds"))

# Sourcing some functions used for plotting
# DEPENDENCY: figures/main_figure_functions.R
source("figures/main_figure_functions.R")

# Reading the data for the LD plot

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

# Drawing the figure in a PNG device
png(paste0("figures/", trait, "_main_figure.png"), width = 6, height = 6.5, units = "in", res = 400)

# Resetting the plotting page
grid.newpage()

# Creating the layout for the whole figure
pushViewport(viewport(layout = grid.layout(nrow = 2, ncol = 1, heights = unit(c(1, 5.5), "in")), name = "main"))

# Drawing the genome-wide Manhattan plots
pushViewport(viewport(layout.pos.row = 1, name = "gwide_manhattan"))
grid.text("(a)", x = 0.02, y = 0.94)
draw_manhattan(list(kmers_gwide_manhattan), 
	       sigline_regexp = "sigline_[1356]",
	       siglabel_regexp = "siglabel_[1356]",
	       labels = "", fontsize = 10)
# Adjusting the position of the label for stGm16
#grid.edit("manhattan_siglabel_1", hjust = 1.2)
grid.edit("manhattan_siglabel_1", y = unit(3.5, "native"))
grid.edit("manhattan_siglabel_3", y = unit(8, "native"))
#grid.edit("manhattan_siglabel_5", hjust = 1.2)
grid.edit("manhattan_siglabel_6", y = unit(3.5, "native"))
grid.edit("manhattan_siglabel_[1-6]", grep = TRUE, global = TRUE, gp = gpar(fontsize = 10), hjust = 1.2)


# Drawing the LD plot
seekViewport("main")
pushViewport(viewport(layout.pos.row = 2, name = "ld_plot"))
grid.text("(b)", x = 0.02, y = 0.90)
pushViewport(viewport(width = unit(5.3, "in"), height = unit(5.3, "in")))
ld_plot(gsorted_ld, kmer_positions, top_legend = FALSE, ylabels = TRUE, ylabel_pattern = "^Gm[0-9]{2}$", fontsize = 8)

dev.off()

