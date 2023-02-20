# Code for the figure on pubescence density to include in the main text
# of the manuscript

# Loading the required libraries
suppressMessages(library(grid))
suppressMessages(library(gwask))
suppressMessages(library(GenomicRanges))
suppressMessages(library(Rsamtools))
suppressMessages(library(parallel))
suppressMessages(library(GenomicAlignments))

# Setting the trait and locus for this series of plots
locus <- "pubescence_density_Ps"
trait <- "pubescence_density"

# Loading the grid objects to plot
# DEPENDENCY: figures/grobs/kmers_pubescence_density_manhattan.rds
kmers_gwide_manhattan <- readRDS(paste0("figures/grobs/kmers_", trait, "_manhattan.rds"))

# DEPENDENCY: figures/grobs/kmers_pubescence_density_Ps_gene.rds
kmers_zoomed_manhattan <- readRDS(paste0("figures/grobs/kmers_", locus, "_gene.rds"))

# Sourcing some functions used for plotting
# DEPENDENCY: figures/main_figure_functions.R
source("figures/main_figure_functions.R")

# Computing the average coverage as well as the 25% and 75% quantiles for normal and semisparse pubescence

# Getting the set of samples to analyse
# DEPENDENCY: utilities/srr_id_correspondence.txt
phenodata <- read.table("phenotypic_data/phenotypic_data.csv", sep = ";", header = TRUE)[, c("bayer_id", "pubescence_density")]
phenodata <- phenodata[complete.cases(phenodata), ]
samples <- phenodata$bayer_id
phenotype <- phenodata$pubescence_density

# Getting the names of the bam files to process
# DEPENDENCY: merged bam files for all samples
bam_files <- paste0("illumina_data/merged_bams/", samples, "_merged.bam")

# Getting a vector of average coverage per sample (pre-computed)
# DEPENDENCY: manifest files for all samples
average_cov <- sapply(samples, function(x) {
			      read.table(paste0("sv_genotyping/paragraph/manifest_files/", x, "_manifest.txt"),
					 skip = 1)[[3]]})

# Creating a GRanges object with the coordinates that we are interested in
cnv_range <- GRanges(seqnames = "Gm12", ranges = IRanges(start = kmers_zoomed_manhattan$vp$xscale[1] - 100,
							 end = kmers_zoomed_manhattan$vp$xscale[2] + 100))

cov_vector <- parallel::mclapply(bam_files, function(x, grange) {
					 bam_file <- Rsamtools::BamFile(x)
					 grange_coverage <-
						 GenomicAlignments::coverage(bam_file,
									     param = Rsamtools::ScanBamParam(what = Rsamtools::scanBamWhat(),
													     which = grange))
					 coverage_rle <- grange_coverage[grange]
					 as.numeric(coverage_rle[[1]]) },
					 grange = cnv_range, mc.cores = 30)

normalized_cov <- mapply(function(x, y) x / y, x = cov_vector, y = average_cov)

roll_averaged_cov <- apply(normalized_cov, 2, function(x) zoo::rollmean(x, k = 100, fill = 1))
colnames(roll_averaged_cov) <- samples
roll_averaged_cov <- as.data.frame(roll_averaged_cov)
roll_averaged_cov$pos <- (kmers_zoomed_manhattan$vp$xscale[1] - 100):(kmers_zoomed_manhattan$vp$xscale[2] + 100)

coverage_data <- list(normal = roll_averaged_cov[, samples[phenotype == 1]], semisparse = roll_averaged_cov[, samples[phenotype == 2]])

average <- lapply(coverage_data, function(x) apply(x, 1, median))
lower_bound <- lapply(coverage_data, function(x) apply(x, 1, quantile, probs = 0.25))
upper_bound <- lapply(coverage_data, function(x) apply(x, 1, quantile, probs = 0.75))


# Recovering the positions of the read containing the most significant k-mer for plotting
semisparse_samples <- samples[phenotype == 2]

# For each of these samples, we get the positions of the reads containing the k-mer
# DEPENDENCY: mapped reads extracted from katcher # not yet included in Makefile
kmer_mapping_pos <- lapply(semisparse_samples, function(x) {
				   filename <- paste0("gwas_results/kmer_data/pubescence_density/katcher_results/", x, "/", x, "_pvalues_sorted.bam")
				   records <- system(paste0("samtools view ", filename, " | grep AAAATTAATGATATTTTTTTGTAAAAATTAT | cut -f4-5,9"), intern = TRUE)
				   if(length(records)) {
					   records <- do.call("rbind", strsplit(records, "\t"))
					   records <- as.data.frame(records)
					   names(records) <- c("POS", "MAPQ", "INSERT")
					   records$POS <- as.numeric(records$POS)
					   records$MAPQ <- as.numeric(records$MAPQ)
					   records$INSERT <- as.numeric(records$INSERT)
					   records} else {
						   records <- data.frame(POS = numeric(), MAPQ = numeric(), INSERT = numeric())
					 }
				   return(records)
					 })

kmer_mapping_pos <- do.call("rbind", kmer_mapping_pos)
kmer_mapping_pos <- kmer_mapping_pos[kmer_mapping_pos$MAPQ > 0, ]
kmer_hist_data <- (hist(kmer_mapping_pos$POS, breaks = 35, plot = FALSE))

# Drawing the figure in a PNG device
png(paste0("figures/", locus, "_main_figure.png"), width = 6, height = 5.8, units = "in", res = 200)

# Resetting the plotting page
grid.newpage()

# Creating the layout for the whole figure
pushViewport(viewport(layout = grid.layout(nrow = 8, ncol = 1, heights = unit(c(4, 1.8, 4, 0.2, 4, 0.2, 4, 1.8), "null")),
		      name = "main"))

# Drawing the genome-wide Manhattan plots
pushViewport(viewport(layout.pos.row = 1, name = "gwide_manhattan"))
grid.text("(a)", x = 0.02, y = 0.95)
draw_manhattan(list(kmers_gwide_manhattan), 
	       sigline_regexp = "sigline_3",
	       siglabel_regexp = "siglabel_3",
	       labels = "k-mers",
	       fontsize = 10)


# Drawing the zoomed-in Manhattan plots
seekViewport("main")
pushViewport(viewport(layout.pos.row = 3, name = "zoomed_manhattan"))
grid.text("(b)", x = 0.02, y = 0.95)
draw_zoomed(list(kmers_zoomed_manhattan), labels = "k-mers", fontsize = 10, xaxis = FALSE)

# Push a viewport to plot the coverage data
seekViewport("main")
pushViewport(viewport(layout.pos.row = 5))
grid.text("(c)", x = 0.02, y = 0.95)
pushViewport(plotViewport(c(0.5, 4.5, 0.5, 0.5), xscale = kmers_zoomed_manhattan$vp$xscale, yscale = c(0, 6.9)))
grid.yaxis(gp = gpar(fontsize = 10))
grid.rect()

# Plotting the data
grid.lines(x = unit(roll_averaged_cov$pos, "native"),
	   y = unit(average$normal, "native"),
	   gp = gpar(col = "blue"))
grid.polygon(x = unit(c(roll_averaged_cov$pos, rev(roll_averaged_cov$pos)), "native"),
	     y = unit(c(lower_bound$normal, rev(upper_bound$normal)), "native"),
	     gp = gpar(fill = "blue", alpha = 0.3))

grid.lines(x = unit(roll_averaged_cov$pos, "native"),
	   y = unit(average$semisparse, "native"),
	   gp = gpar(col = "red"))
grid.polygon(x = unit(c(roll_averaged_cov$pos, rev(roll_averaged_cov$pos)), "native"),
	     y = unit(c(lower_bound$semisparse, rev(upper_bound$semisparse)), "native"),
	     gp = gpar(fill = "red", alpha = 0.3))

# Adding the color legend
grid.rect(x = 0.02, y = 0.90, width = 0.03, height = 0.05, just = c(0, 0),
	  gp = gpar(fill = "red"))
grid.text("Semi-sparse pubescence", x = 0.06, y = 0.925, hjust = 0,
	  gp = gpar(fontsize = 8))

grid.rect(x = 0.02, y = 0.81, width = 0.03, height = 0.05, just = c(0, 0),
	  gp = gpar(fill = "blue"))
grid.text("Normal pubescence", x = 0.06, y = 0.835, hjust = 0,
	  gp = gpar(fontsize = 8))

# Adding the y-axis label
grid.text("Copy number", x = unit(-3, "lines"), rot = 90, gp = gpar(fontsize = 10))

# Plotting the positions of the mapped reads
seekViewport("main")
pushViewport(viewport(layout.pos.row = 7))
grid.text("(d)", x = 0.02, y = 0.95)

# Creating a function that plots a histogram using grid functions
grid.histogram <- function(x, xscale, xlabel, ylabel, fontsize) {
	# Setting the plotting region and scale
	pushViewport(plotViewport(c(0.5, 4.5, 0.5, 0.5)))
	pushViewport(dataViewport(yData = x$counts, xscale = xscale, extension = 0.03))

	# Drawing the axes and box
	tick_positions <- axisTicks(usr = xscale, log = FALSE, nint = 5)
	grid.xaxis(at = tick_positions, label = tick_positions / 10^6, gp = gpar(fontsize = fontsize))
	grid.yaxis(gp = gpar(fontsize = fontsize))
	grid.rect()

	# Drawing the histogram bars
	grid.rect(x = x$mids,
		  width = diff(x$breaks),
		  y = 0, height = x$counts,
		  just = c(0.5, 0),
		  default.units = "native",
		  gp = gpar(fill = "gray70"))

	# Adding the x- and y-axis label
	grid.text(xlabel, y = unit(-3, "lines"), gp = gpar(fontsize = fontsize))
	grid.text(ylabel, x = unit(-3, "lines"), rot = 90, gp = gpar(fontsize = fontsize))

	upViewport(2)

}

grid.histogram(kmer_hist_data, xscale = kmers_zoomed_manhattan$vp$xscale,
	       xlabel = "Position along Gm12 (Mb)", ylabel = "Number of reads",
	       fontsize = 10)

dev.off()

