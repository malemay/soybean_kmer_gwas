# Loading the required packages
suppressMessages(library(grid))
suppressMessages(library(gwask))
suppressMessages(library(GenomicRanges))

# Reading the trait being analyzed from the command line
trait <- commandArgs(trailingOnly = TRUE)[1]

programs <- c("platypus", "paragraph", "kmers")
program_labels <- c("SNPs/indels", "SVs", "kmers")

# Outputting to a PNG file
png(paste0("figures/", trait, "_manhattan.png"), width = 10, height = 10, units = "in", res = 100)

# Resetting the viewport
grid.newpage()

# DEPENDENCY: Paragraph manhattan subplot
# DEPENDENCY: Platypus manhattan subplot
# DEPENDENCY: kmers manhattan subplot

# Setting a viewport with 3 rows, one for each of the programs
pushViewport(viewport(height = 0.95, y = 0, just = "bottom",
		      layout = grid.layout(nrow = 3, ncol = 2, widths = unit(c(72, 28), "null"))))

# Reading and printing each of the subplots in turn in their viewport row
for(i in 1:3){
	pushViewport(viewport(layout.pos.row = i, layout.pos.col = 1))
	grid.draw(readRDS(paste0("figures/grobs/", programs[i], "_", trait, "_manhattan.rds")))
	if(programs[i] == "kmers") {
		grid.text(expression(paste("(c) ", italic(k)-mers)), 0.1, 1.08, hjust = 0, gp = gpar(fontsize = 14))
	} else {
		grid.text(paste0("(", letters[i], ") ", program_labels[i]), 0.1, 1.08, hjust = 0, gp = gpar(fontsize = 14))
	}
	upViewport()
}

# Taking care of the qqplots
for(i in 1:3) {
	# Generating the data to plot
	if(programs[i] == "kmers") {
		pvalues <- scan(paste0("gwas_results/kmer_data/", trait, "/all_pvalues_sorted.txt"),
				what = numeric(1),
				sep = "\n")
	} else {
		# Formatting the p-values
		pvalues <- readRDS(paste0("gwas_results/", programs[i], "/", trait, "_gwas.rds"))$log10p
		pvalues <- sort(pvalues, decreasing = TRUE)
	}

		# Formatting the expected values
		expected <- -log10(ppoints(length(pvalues)))

		# Setting the viewports
		pushViewport(viewport(layout.pos.row = i, layout.pos.col = 2))
		pushViewport(plotViewport(margins = c(5.1, 4.1, 0.1, 0.1)))
		pushViewport(dataViewport(xData = expected, yData = pvalues))
		pushViewport(dataViewport(xData = expected, yData = pvalues, clip = "on"))

		# Doing the plotting
		grid.rect()
		grid.abline(intercept = 0, slope = 1, gp = gpar(col = "red", lty = 2))
		grid.points(x = expected, y = pvalues, pch = ".", default.units = "native")

		# Going up 1 viewport because we need clipping disabled
		upViewport(1)
		grid.xaxis() ; grid.yaxis()
		grid.text(label = expression(paste(-log[10], "(expected ", italic(p), "-value)")), y = unit(-3, "lines"))
		grid.text(label = expression(paste(-log[10], "(observed ", italic(p), "-value)")), x = unit(-3, "lines"), rot = 90)

		# Going back to the main viewport
		upViewport(3)
}

dev.off()

