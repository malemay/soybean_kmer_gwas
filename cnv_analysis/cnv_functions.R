# A set of functions used for the analysis of CNV regions
# at specific loci; used to delimitate the region where
# the CNV occurs

# A function used for segmentation
var_diff <- function(x) {
	n <- length(x)

	wei <- c(0.1942, 0.2809, 0.3832, -0.8582)

	mat <- wei %*% t(x)
	mat[2, -n] <- mat[2, -1]
	mat[3, -c(n-1, n)] <- mat[3, -c(1, 2)]
	mat[4, -c(n-2, n-1, n)] <- mat[4, -c(1, 2, 3)]

	sqrt(sum(apply(mat[, -c(n-2, n-1, n)], 2, sum)^2) / (n-3))
}

# A function that extracts the average coverage given a bam file and a GRanges object
get_coverage <- function(bam_path, grange) {
	message("Processing ", bam_path)
	bam_file <- Rsamtools::BamFile(bam_path)
	grange_coverage <- GenomicAlignments::coverage(bam_file,
						       param = Rsamtools::ScanBamParam(what = Rsamtools::scanBamWhat(),
									    which = grange))

	grange_coverage <- grange_coverage[grange]
	average_cov <- sum(grange_coverage) / GenomicRanges::width(grange)

	return(average_cov)
}

# A function that plots the normalized coverage over a given GRanges object for a set of samples
plot_coverage <- function(bam_paths, grange, average_coverage, yscale = c(0, 20), cores = 10) {

	cov_vector <- parallel::mclapply(bam_paths, function(x, grange) {
						 bam_file <- Rsamtools::BamFile(x)
						 grange_coverage <-
							 GenomicAlignments::coverage(bam_file,
										     param = Rsamtools::ScanBamParam(what = Rsamtools::scanBamWhat(),
														     which = grange))
						 coverage_rle <- grange_coverage[grange]
						 as.numeric(coverage_rle[[1]]) },
					  grange = grange, mc.cores = cores)

	normalized_cov <- mapply(function(x, y) x / y, x = cov_vector, y = average_coverage)

	# Initializing a plotting region with appropriate scales
	plot(1, 1, xlim = c(GenomicRanges::start(grange), GenomicRanges::end(grange)), ylim = yscale, type = "n")

	# Adding one coverage line per sample
	for(i in 1:ncol(normalized_cov)) {
		lines(x = start(grange):end(grange), y = normalized_cov[, i], col = sample(colors(), 1))
	}

	normalized_cov
}

# A function that gets the total coverage on a genomic range for a set of samples
bam_cov <- function(bam_paths, grange, cores) {
	coverage_matrix <- parallel::mclapply(cnv_bampaths, function(x, grange) {
					    bam_file <- Rsamtools::BamFile(x)
					    grange_coverage <-
						    GenomicAlignments::coverage(bam_file,
										param = Rsamtools::ScanBamParam(what = Rsamtools::scanBamWhat(),
														which = grange))
					    grange_coverage <- grange_coverage[grange]
					    as.numeric(grange_coverage[[1]])},

					    grange = grange,
					    mc.cores = cores)

	coverage_matrix <- do.call("cbind", coverage_matrix)
	return(rowSums(coverage_matrix))
}

# A function that performs segmentation of a coverage vector using the given parameters
segment <- function(values, threshold_param, lambda_param) {
	# Getting the segmentation parameter
	sd_estimate <- var_diff(values)
	threshold <- threshold_param * sd_estimate

	# Performing the segmentation itself
	segments <- robseg::Rob_seg.std(x = values,
					loss = "Outlier",
					lambda = lambda_param * log(length(values)) * sd_estimate,
					lthreshold = threshold)

	segments
}

