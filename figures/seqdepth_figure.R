#!/home/malem420/.local/bin/Rscript

# This figure shows the proportion of the genome covered by
# only one read as a function of mean coverage for a given sample

# Loading the required libraries
library(parallel)

# Reading the mean coverage per sample
# DEPENDENCY: phenotypic_data/phenotypic_data.csv
samples <- read.table("phenotypic_data/phenotypic_data.csv", sep = ";", header = TRUE)$bayer_id

# DEPENDENCY: sv_genotyping/paragraph/MANIFEST_FILES
average_cov <- sapply(samples, function(x) {
			      read.table(paste0("sv_genotyping/paragraph/manifest_files/", x, "_manifest.txt"),
					 skip = 1)[[3]]})

# Reading the coverage
# DEPENDENCY: illumina_data/SAMTOOLS_STATS
stats_files <- paste0("illumina_data/merged_bams/", samples, "_stats.txt")
stopifnot(all(file.exists(stats_files)))

read_coverage <- function(stats) {
	output <- readLines(stats)
	output <- grep("^COV", output, value = TRUE)
	output <- do.call("rbind", strsplit(output, "\t"))
	output <- as.data.frame(output)[, -1]
	names(output) <- c("range", "mean", "count")
	output$mean <- as.numeric(output$mean)
	output$count <- as.numeric(output$count)
	output
}

one_proportion <- function(stats) {
	coverage <- read_coverage(stats)
	stopifnot(sum(coverage$mean == 1) == 1)
	coverage[coverage$mean == 1, "count"] / sum(coverage$count)
}

# The following vector contains the proportion of positions in the genome that are covered by only one read
kmer_one <- unlist(mclapply(stats_files, one_proportion, mc.cores = 10))

# Saving to file
png("figures/seqdepth.png", width = 6, height = 6, units = "in", res = 300)
plot(average_cov, kmer_one, main = "",
     xlab = "Mean sequencing depth of sample",
     ylab = "Fraction of genome covered by one read")
dev.off()

