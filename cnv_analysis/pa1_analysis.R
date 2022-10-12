# This script aims to see whether there could be copy number variation
# at the Pa1 locus
# 1- Identify groups of samples with differing copy numbers
#    at that locus based on a normalized value of coverage
#    over the region divided by the average coverage for
#    that sample.
# 2- Add together the coverage values in that region (and
#    extending into non-duplicated regions) for all samples
#    with higher coverage.
# 3- Subject the coverage values for duplicated samples to
#    a segmentation algorithm such that the location of the
#    copy number variation can be identified as precisely
#    as possible

# Loading the required libraries
suppressMessages(library(Rsamtools))
suppressMessages(library(GenomicRanges))
suppressMessages(library(parallel))
suppressMessages(library(robseg))
suppressMessages(library(GenomicAlignments))

# Sourcing some functions used for this analysis
# DEPENDENCY: cnv_analysis/cnv_functions.R
source("cnv_analysis/cnv_functions.R")

# Parameters for the analysis
gene_name <- tolower("Glyma.12g213900")

# Getting the set of samples to analyse
# DEPENDENCY: utilities/srr_id_correspondence.txt
samples <- read.table("utilities/srr_id_correspondence.txt")[[1]]

# Getting the names of the bam files to process
# DEPENDENCY: merged bam files for all samples
bam_files <- paste0("illumina_data/merged_bams/", samples, "_merged.bam")

# Getting a vector of average coverage per sample (pre-computed)
# DEPENDENCY: manifest files for all samples
average_cov <- sapply(samples, function(x) {
			      read.table(paste0("sv_genotyping/paragraph/manifest_files/", x, "_manifest.txt"),
					 skip = 1)[[3]]})

# Loading the set of genes annotated in the reference and extracting the coordinates for the gene of interest
# DEPENDENCY: refgenome/gmax_v4_genes.rds
genes <- readRDS("refgenome/gmax_v4_genes.rds")
genes[gene_name]

# Plotting the coverage over a window larger than the expected to get an idea as to where the coverage shifts
if(interactive()) {
	covplot <- plot_coverage(bam_files, average_coverage = average_cov,
				 yscale = c(0, 10), cores = 30,
				 grange = GRanges(seqnames = "Gm12", ranges = IRanges(start = 38675000, end = 38775000)))
}

# I will use the window plotted above to calculate the average normalized coverage for all samples
rcoverage <- mclapply(bam_files, function(x, grange) get_coverage(x, grange),
		      grange = GRanges(seqnames = "Gm12", ranges = IRanges(start = 38685000, end = 38750000)),
		      mc.cores = 30)

# there seems to be copy number variation going roughly from 36,680,000 to 36,750,000
# Let us quantify copy number and see if there is indeed such variation

# Getting a vector of normalized coverage
rcoverage <- unlist(rcoverage)
names(rcoverage) <- samples
stopifnot(identical(names(rcoverage), names(average_cov)))
normalized_cov <- rcoverage / average_cov

# Visualizing the distribution to identify a threshold separating CNV samples from low copy number samples
if(interactive()) hist(normalized_cov, breaks = 50)
# A threhsold of 1.2 can be clearly used to identify samples with high coverage

# Let us try and see whether there is a link between coverage in that region and the phenotype
phenotypic_data <- read.table("phenotypic_data/phenotypic_data.csv", header = TRUE, sep = ";")
pubform <- phenotypic_data[match(names(normalized_cov), phenotypic_data$bayer_id), "pubescence_form_all"]

boxplot_data <- data.frame(coverage = normalized_cov,
			   pubform = pubform,
			   stringsAsFactors = FALSE)

boxplot_data$pubform <- as.factor(boxplot_data$pubform)

boxplot(coverage ~ pubform, data = boxplot_data)
summary(aov(coverage ~ pubform, data = boxplot_data))

# There is no link between copy number and phenotype in that region
