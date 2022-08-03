# Loading the required packages
suppressMessages(library(grid))
suppressMessages(library(gwastools))

# Reading the trait being analyzed from the command line
trait <- commandArgs(trailingOnly = TRUE)[1]

programs <- c("platypus", "vg", "paragraph", "kmers")

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
for(i in 1:4){
	pushViewport(viewport(layout.pos.row = i))
	grid.draw(readRDS(paste0("figures/ggplots/", programs[i], "_", trait, "_manhattan.rds")))
	grid.text(programs[i], 0.1, 0.9, hjust = 0)
	upViewport()
}

dev.off()

