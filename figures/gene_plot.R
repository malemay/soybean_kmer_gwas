# Loading the required packages
library(grid)

# Reading the ID being analyzed from the command line
id <- commandArgs(trailingOnly = TRUE)[1]

programs <- c("platypus", "vg", "paragraph", "kmers")

# Outputting to a PNG file
png(paste0("figures/", id, "_gene.png"), width = 9, height = 10, units = "in", res = 300)

# Resetting the viewport
grid.newpage()

# DEPENDENCY: Paragraph gene subplot
# DEPENDENCY: vg gene subplot
# DEPENDENCY: Platypus gene subplot
# DEPENDENCY: kmers gene subplot

# Setting a viewport with 4 rows, one for each of the programs
pushViewport(viewport(layout = grid.layout(nrow = 4)))

# Reading and printing each of the subplots in turn in their viewport row
for(i in 1:4){
	pushViewport(viewport(layout.pos.row = i))
	grid.draw(readRDS(paste0("figures/ggplots/", programs[i], "_", id, "_gene.rds")))
	grid.text(programs[i], 0.1, 0.9, hjust = 0)
	upViewport()
}

dev.off()

