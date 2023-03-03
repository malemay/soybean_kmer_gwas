# Loading the required packages
suppressMessages(library(grid))
suppressMessages(library(gwask))

# Reading the trait being analyzed from the command line
trait <- commandArgs(trailingOnly = TRUE)[1]

programs <- c("platypus", "paragraph", "kmers")
program_labels <- c("SNPs/indels", "SVs", "kmers")

# Outputting to a PNG file
png(paste0("figures/", trait, "_manhattan.png"), width = 9, height = 10, units = "in", res = 100)

# Resetting the viewport
grid.newpage()

# DEPENDENCY: Paragraph manhattan subplot
# DEPENDENCY: Platypus manhattan subplot
# DEPENDENCY: kmers manhattan subplot

# Setting a viewport with 3 rows, one for each of the programs
pushViewport(viewport(height = 0.95, y = 0, just = "bottom", layout = grid.layout(nrow = 3)))

# Reading and printing each of the subplots in turn in their viewport row
for(i in 1:3){
	pushViewport(viewport(layout.pos.row = i))
	grid.draw(readRDS(paste0("figures/grobs/", programs[i], "_", trait, "_manhattan.rds")))
	if(programs[i] == "kmers") {
		grid.text(expression(paste("(c) ", italic(k)-mers)), 0.1, 1.08, hjust = 0, gp = gpar(fontsize = 14))
	} else {
		grid.text(paste0("(", letters[i], ") ", program_labels[i]), 0.1, 1.08, hjust = 0, gp = gpar(fontsize = 14))
	}
	upViewport()
}

dev.off()

