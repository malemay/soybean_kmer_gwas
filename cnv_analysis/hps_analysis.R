# This script aims to delimitate the extent of the duplicated
# region at the Hps locus. This will done in the following steps:
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
gene_name <- "glyma.15g119800"

# Getting the set of samples to analyse
# DEPENDENCY: utilities/srr_id_correspondence.txt
samples <- read.table("utilities/srr_id_correspondence.txt")[[1]]

# Getting the names of the bam files to process
# DEPENDENCY: bam files for all samples
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
				 yscale = c(0, 15), cores = 30,
				 grange = GRanges(seqnames = "Gm15", ranges = IRanges(start = 9380000, end = 9430000)))
}

# I will use the window plotted above to calculate the average normalized coverage for all samples
rcoverage <- mclapply(bam_files, function(x, grange) get_coverage(x, grange),
		      grange = GRanges(seqnames = "Gm15", ranges = IRanges(start = 9395000, end = 9415000)),
		      mc.cores = 30)

# Getting a vector of normalized coverage
rcoverage <- unlist(rcoverage)
names(rcoverage) <- samples
stopifnot(identical(names(rcoverage), names(average_cov)))
normalized_cov <- rcoverage / average_cov

# Visualizing the distribution to identify a threshold separating CNV samples from low copy number samples
if(interactive()) hist(normalized_cov, breaks = 50)

# Having a look at the link between copy number and phenotype
phenotypic_data <- read.table("phenotypic_data/phenotypic_data.csv", header = TRUE, sep = ";")
luster <- phenotypic_data[match(names(normalized_cov), phenotypic_data$bayer_id), "seed_coat_luster_all"]

boxplot_data <- data.frame(coverage = normalized_cov,
			   luster = luster,
			   stringsAsFactors = FALSE)

boxplot_data$luster <- as.factor(boxplot_data$luster)

boxplot(coverage ~ luster, data = boxplot_data)
summary(aov(coverage ~ luster, data = boxplot_data))
TukeyHSD(aov(coverage ~ luster, data = boxplot_data))

# Also plotting the results with only the dull and shiny phenotypes
boxplot_data$dullshiny <- phenotypic_data[match(names(normalized_cov), phenotypic_data$bayer_id), "seed_coat_luster_dullshiny"]
boxplot(coverage ~ dullshiny, data = boxplot_data)
t.test(coverage ~ dullshiny, data = boxplot_data)

# It is weird that some accessions have shiny seed coats, yet high copy number at the B locus
# Let us make a GWAS contrasting only accessions with high copy number
cnv_luster_data <- boxplot_data[complete.cases(boxplot_data), ]
cnv_luster_data <- cnv_luster_data[cnv_luster_data$coverage > 2,]
cnv_luster_data$id <- rownames(cnv_luster_data)
write.table(cnv_luster_data, file = "cnv_analysis/cnv_luster_data.csv", sep = ",",
	    col.names = TRUE, row.names = FALSE, quote = FALSE)

# Keeping the samples with > 2.5 average coverage over the region
cnv_samples <- names(normalized_cov[normalized_cov > 2.5])
cnv_bampaths <- paste0("illumina_data/merged_bams/", cnv_samples, "_merged.bam")

# Computing the total coverage per position over a larger range for those samples
cnv_coverage <- bam_cov(cnv_bampaths,
			grange = GRanges(seqnames = "Gm15", ranges = IRanges(start = 9350000, end = 9450000)),
			cores = 30)

# Normalizing this coverage based on the total average coverage
normalized_cnv_coverage <- cnv_coverage / sum(average_cov[cnv_samples])

# Making a data.frame out of those positions and normalized coverage
coverage_df <- data.frame(xpos = 9350000:9450000,
			  coverage = normalized_cnv_coverage,
			  stringsAsFactors = FALSE)

segments <- segment(coverage_df$coverage, threshold_param = 4, lambda_param = 2)
coverage_df$segment_value <- segments$smt

if(interactive()) {
	plot(coverage_df$xpos, coverage_df$coverage, type = "l", ylim = c(0, 5))
	lines(x = coverage_df$xpos, y = coverage_df$segment_value, col = "blue")
}

# Identifying the first and last position with segments over 2.5X
coverage_df <- coverage_df[coverage_df$segment_value > 2.5, ]
head(coverage_df)
tail(coverage_df)

cnv_range <- GRanges(seqnames = "Gm15", ranges = IRanges(start = coverage_df[1, "xpos"],
							 end = coverage_df[nrow(coverage_df), "xpos"]))

# Saving that GRanges object to a file
saveRDS(cnv_range, file = "cnv_analysis/hps_cnv_range.rds", compress = FALSE)

