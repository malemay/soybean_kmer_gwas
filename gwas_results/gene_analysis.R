# This code generates four output files for each locus and program of interest
# - A GRanges object that delimits the top 5% markers among the ones that are significantly associated at a given locus
# - A .tsv file with the list of all genes encompassed by the signal and relevant metadata
# - A .tsv file with the list of all genes encompassed by the top 5% region and relevant metadata
# - A .tsv file with containing the gene that is nearest to the most associated variant

# Loading the required libraries
suppressMessages(library(GenomicRanges))

# Getting the name of the trait-locus combination and of the program that called the genotypes from the command line
id <- commandArgs(trailingOnly = TRUE)[1]
program <- commandArgs(trailingOnly = TRUE)[2]

# The distance at which signals overlapping the reference signal will be merged
gapwidth <- 10^6

# Getting the ID of the gene associated with the locus
# DEPENDENCY: utilities/all_signals.rds
target_signal <- readRDS("utilities/all_signals.rds")
target_signal <- target_signal[id]
stopifnot(length(target_signal) == 1)

# Getting the name of the trait associated with the signal, and from it the significance threshold
trait <- target_signal$trait
# DEPENDENCY: significance thresholds (should not need to be included in Makefile, because the dependencs is implicit)
threshold <- -log10(as.numeric(readLines(paste0("gwas_results/", program, "/", trait, "_threshold_5per.txt"))))

# Extracting the gene name and handling the case where several genes are semicolon-separated
gene_name <- target_signal$gene_name_v4

if(is.na(gene_name)) {
	gene_name <- ""
} else {
	gene_name <- strsplit(gene_name, ";")[[1]]
}

### ----------

# Getting the set of reference genes
# DEPENDENCY: refgenome/gmax_v4_genes.rds
genes <- readRDS("refgenome/gmax_v4_genes.rds")

# Reading the reference annotations
# DEPENDENCY: refgenome/soybase_gmax_v4_annotations.rds
gmax_annotation <- readRDS("refgenome/soybase_gmax_v4_annotations.rds")

# Reading the GWAS results and signals from file
# DEPENDENCY: GWAS signals and marker subset
gwas_results <- readRDS(paste0("gwas_results/", program, "/", id, "_gwas_locus.rds"))
gwas_signals <- readRDS(paste0("gwas_results/", program, "/", id, "_signal_locus.rds"))

# Extracting the signal for the locus of interest from those just read
locus_signal <- subsetByOverlaps(gwas_signals, target_signal, ignore.strand = TRUE)
# Merging this signal
locus_signal <- reduce(locus_signal, min.gapwidth = gapwidth, ignore.strand = TRUE)
stopifnot(length(locus_signal) <= 1)

# Subsetting the gwas_results using this signal keeping only significant hits, and sorting by p-value
gwas_markers <- subsetByOverlaps(gwas_results, locus_signal, ignore.strand = TRUE)
gwas_markers <- gwas_markers[gwas_markers$log10p >= threshold]
gwas_markers <- gwas_markers[order(gwas_markers$log10p, decreasing = TRUE)]

# A function that generates a GRanges object of the region encompassing the top associated markers
#' @param markers A Granges object of significantly associated markers. They must have already been subset to a given region
#' @param fraction A numeric between 0 and 1. The fraction of the markers that should be kept (e.g. 0.05 for the top 5% markers)
#'                 the function will minimally return two markers. If there is only one marker or no marker at all, then
#'                 the input markers are returned as is.
top_region <- function(markers, fraction) {
	# Returning an empty GRanges if it is already empty or only has one marker
	if(length(markers) < 2) return(GRanges())
	
	# Getting the number of markers to extract
	top_n <- ceiling(length(markers) * fraction)
	top_n <- max(top_n, 2)

	# Sorting the markers and selecting the top ones
	markers <- markers[order(markers$log10p, decreasing = TRUE)]
	markers <- markers[1:top_n]

	# Reducing the GRanges object to a single region
	output <- reduce(markers, min.gapwidth = 10^7, ignore.strand = TRUE)

	return(output)
}

# A function that returns the top marker in a region
#' @param markers A GRanges object describing a set of markers and their p-values. Should be restricted to a region of interest.
#'                 An empty input will return an empty output.
top_marker <- function(markers) {
	# Returning an empty GRanges object if the input was empty
	if(!length(markers)) return(GRanges())

	# Sort the markers and return the top one
	markers <- markers[order(markers$log10p, decreasing = TRUE)]
	return(markers[1])
}

# A function that subsets the gene annotation data.frame to only include genes within some range
#' @param region:  a GRanges object of length one denoting the region of interest
#' @param markers: a GRanges object containing the set of significant markers at a locus
#' @param genes: a GRanges object of annotated genes
#' @param annotation: a data.frame of gene annotations for the genes in the "genes" object
#' @param causal_genes: a character vector of causal gene(s) name(s) for the locus under investigation
#' @param nearest_only: whether we should extract the gene nearest to the signal, or all overlapping genes
subset_annotation <- function(region, markers, genes, annotation, causal_genes, nearest_only = FALSE) {

	# Creating a template for the metadata to extract from the overlapped genes
	metadata_template <- data.frame(seqnames = character(),
					start = numeric(),
					end = numeric(),
					width = numeric(),
					strand = character(),
					causal_gene = logical(),
					marker_overlap = logical(),
					marker_rank = numeric(),
					locus = character(),
					desc_uniref = character(),
					go_bp_term = character(),
					go_mf_term = character(),
					go_cc_term = character(),
					pfam_id = character(),
					pfam_name = character(),
					alt_tair_id = character(),
					stringsAsFactors = FALSE)

	# If region is empty then we return the template
	if(!length(region)) return(metadata_template)

	# Otherwise the region needs to be a single range
	stopifnot(length(region) == 1)

	# Extract the genes that overlap that region, or the nearest gene only
	if(nearest_only) {
		subset_genes <- genes[nearest(region, genes, ignore.strand = TRUE)]
	} else {
		subset_genes <- subsetByOverlaps(genes, region, ignore.strand = TRUE)
	}

	# If no genes remain then we also return the metadata template
	if(!length(subset_genes)) return(metadata_template)

	# Subset the annotation based on the subset of genes
	subset_annot <- annotation[annotation$locus %in% names(subset_genes), ]

	# Also keeping only the metadata columns of interest
	subset_annot <- subset_annot[, colnames(subset_annot) %in% colnames(metadata_template)]

	# Adding the information regarding the GenomicRanges object to the data.frame
	subset_annot <- cbind(as.data.frame(genes[subset_annot$locus]), subset_annot)

	# Specifying whether each gene corresponds to the causal gene
	subset_annot$causal_gene <- subset_annot$locus %in% causal_genes

	# Specifying whether at least one marker within region overlaps each gene
	subset_annot$marker_overlap <- overlapsAny(genes[subset_annot$locus], markers, ignore.strand = TRUE)

	# Getting the rank of the most associated marker within each gene
	subset_annot$marker_rank <- NA
	markers <- markers[order(markers$log10p, decreasing = TRUE)]

	for(i in 1:nrow(subset_annot)) {
		if(subset_annot[i, "marker_overlap"]) {
			overlapping_markers <- overlapsAny(markers, genes[subset_annot[i, "locus"]], ignore.strand = TRUE)
			subset_annot[i, "marker_rank"] <- min(which(overlapping_markers))
		}
	}

	# Reordering the columns according to the template
	subset_annot <- subset_annot[, colnames(metadata_template)]

	# Returning the data.frame
	return(subset_annot)
}

top_granges <- top_region(gwas_markers, fraction = 0.05)

# Getting a data.frame of all the genes in the region
all_genes_df <- subset_annotation(region = locus_signal, markers = gwas_markers,
				  genes = genes, annotation = gmax_annotation,
				  causal_genes = gene_name, nearest_only = FALSE)

# Getting a data.frame of the genes overlapping the top 5% region
top_genes_df <- subset_annotation(region = top_granges, markers = gwas_markers,
				  genes = genes, annotation = gmax_annotation,
				  causal_genes = gene_name, nearest_only = FALSE)

# Getting a data.frame of the gene nearest to the top marker
nearest_gene_df <- subset_annotation(region = top_marker(gwas_markers), markers = gwas_markers,
				     genes = genes, annotation = gmax_annotation,
				     causal_genes = gene_name, nearest_only = TRUE)

# Writing the results to files
saveRDS(top_granges, file = paste0("gwas_results/", program, "/", id, "_top_markers.rds"), compress = FALSE)

write.table(all_genes_df, quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE,
	    file = paste0("gwas_results/", program, "/", id, "_all_genes.tsv"))

write.table(top_genes_df, quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE,
	    file = paste0("gwas_results/", program, "/", id, "_top_genes.tsv"))

write.table(nearest_gene_df, quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE,
	    file = paste0("gwas_results/", program, "/", id, "_nearest_gene.tsv"))

