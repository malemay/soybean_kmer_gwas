# Loading the required packages
library(grid)
library(ggplot2)
library(gwastools)

# Reading the trait being analyzed from the command line
trait <- commandArgs(trailingOnly = TRUE)[1]

# Outputting to a PNG file
png(paste0("figures/", trait, "_manhattan.png"), width = 9, height = 10, units = "in", res = 300)

# Resetting the viewport
grid.newpage()

# DEPENDENCY: Paragraph manhattan subplot
# DEPENDENCY: vg manhattan subplot
# DEPENDENCY: Platypus manhattan subplot
# DEPENDENCY: kmers manhattan subplot

# Setting a viewport with 4 rows, one for each of the programs
pushViewport(viewport(layout = grid.layout(nrow = 4)))

# Reading and printing each of the subplots in turn in their viewport row
print(readRDS(paste0("figures/ggplots/platypus_", trait, "_manhattan.rds")), vp = viewport(layout.pos.row = 1))
print(readRDS(paste0("figures/ggplots/vg_", trait, "_manhattan.rds")), vp = viewport(layout.pos.row = 2))
print(readRDS(paste0("figures/ggplots/paragraph_", trait, "_manhattan.rds")), vp = viewport(layout.pos.row = 3))
print(readRDS(paste0("figures/ggplots/kmers_", trait, "_manhattan.rds")), vp = viewport(layout.pos.row = 4))

dev.off()

