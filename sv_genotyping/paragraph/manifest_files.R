# This code prepares the manifest files for use with Paragraph
# The manifest files have the following tab-separated format:
# id	path	depth	read length

# Stats will be obtained from the output of the samtools coverage command for each file

# Setting the working directory
setwd("sv_genotyping/paragraph/")

# Getting a vector of the cultivar names
# DEPENDENCY: utilities/srr_id_correspondence.txt
samples <- read.table("../../utilities/srr_id_correspondence.txt", header = FALSE, stringsAsFactors = FALSE)[[1]]

# Creating and filling a list with the samtools stats output
samtools_stats <- list()

# DEPENDENCY: illumina_data/SAMTOOLS_STATS
for(i in samples) {
	samtools_stats[[i]] <- readLines(paste0("../../illumina_data/merged_bams/", i, "_stats.txt"))
}

# Creating a function that extracts the average read length from a samtools_stats element
parse_length <- function(x) {
	x <- grep("^SN\taverage length:", x, value = TRUE)
	x <- strsplit(x, "\t")[[1]][3]
	as.numeric(x)
}

read_lengths <- sapply(samtools_stats, parse_length)

# DEPENDENCY: illumina_data/SAMTOOLS_COVERAGE
# Reading the output of the samtools coverage command for all samples
samtools_coverage <- list()

for(i in samples) {
	samtools_coverage[[i]] <- read.table(paste0("../../llumina_data/merged_bams/", i, "_coverage.txt"),
					     header = TRUE, stringsAsFactors = FALSE, comment.char = "")
}

# Computing the coverage from the samtools coverage function
compute_mean <- function(x) {
	x <- x[grepl("^Gm[0-9]{2}$", x[[1]]), ]
	weighted.mean(x$meandepth, x$endpos)
}

coverage_means <- sapply(samtools_coverage, compute_mean)

# There will be one manifest file per sample under the manifest_files directory
dir.create("manifest_files", showWarnings = FALSE)

# Printing the manifest file to disk for each sample
for(i in samples) {
	cat("id\tpath\tdepth\tread length\n", file = paste0("manifest_files/", i, "_manifest.txt"))
	cat(paste(i, paste0("../../illumina_data/merged_bams/", i, "_merged.bam"), coverage_means[i], read_lengths[i], "\n", sep = "\t"),
	    file = paste0("manifest_files/", i, "_manifest.txt"),
	    append = TRUE)
}

file.create("MANIFEST_FILES")
