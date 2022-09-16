# Loading the required libraries
suppressMessages(library(GenomicRanges))
suppressMessages(library(gwastools))

# Getting the name of the trait and of the program that called the genotypes from the command line
trait <- commandArgs(trailingOnly = TRUE)[1]
program <- commandArgs(trailingOnly = TRUE)[2]

# A parameter that sets the distance between markers for merging signals with extract_signals
signal_distance <- 250000
platypus_extend <- 50000 # A parameter that sets how much platypus signals are extended for computing the p-values of pruned markers

# Reading the GWAS results and threshold
# DEPENDENCY: GWAS association results
# DEPENDENCY: GWAS thresholds
gwas_results <- readRDS(paste0("gwas_results/", program, "/", trait, "_gwas.rds"))
threshold <- -log10(as.numeric(readLines(paste0("gwas_results/", program, "/", trait, "_threshold_5per.txt"))))

# Extracting the signals found for that trait
gwas_signals <- extract_signals(gwas_results,
				threshold = threshold,
				distance = signal_distance)

# Also generating a subset of the GWAS results that only overlaps signals; will be used for the signal and gene plots
gwas_subset <- subsetByOverlaps(gwas_results, gwas_signals, ignore.strand = TRUE)

# If the program is platypus, we need to find the p-values of the markers that were pruned and add them back to the subset for signal and gene plotting
if(program == "platypus" && length(gwas_signals)) {

	# Adding a column in gwas_subset indicating that those marker were not pruned
	gwas_subset$pruned <- FALSE

	# Extending each signal to the right and left
	start(gwas_signals) <- start(gwas_signals) - platypus_extend
	end(gwas_signals) <- end(gwas_signals) + platypus_extend

	# Opening the VCF file that contains the markers before pruning for reading
	vcf_file <- VariantAnnotation::VcfFile("filtered_variants/platypus_full.vcf.gz")

	# Reading the VCF records found in that region; only keeping those that were not already analyzed
	subset_vcf <- VariantAnnotation::readVcf(vcf_file, param = VariantAnnotation::ScanVcfParam(fixed = NA, info = NA, geno = "GT", which = gwas_signals))
	subset_vcf <- subset_vcf[!names(subset_vcf) %in% names(gwas_subset)]

	# It is only worth computing the GWAS if there are any markers left ; also excluding cases with only one marker as it makes GAPIT fail
	if(length(subset_vcf) > 1) {

		# Reading the kinship and pca results
		kinship <- readRDS("filtered_variants/platypus_gapit_kinship.rds")
		pca <- readRDS("filtered_variants/platypus_gapit_pca.rds")

		gapit_results <- gapit_vcf(vcf = subset_vcf,
					   kinship = kinship,
					   pca = pca,
					   phenodata = "phenotypic_data/phenotypic_data.csv",
					   trait = trait,
					   N.sig = min(length(subset_vcf), 20),
					   id_column = "bayer_id",
					   tmproot = "tmpdir")

		stopifnot(all(gapit_results$SNP %in% names(subset_vcf)))
		GenomicRanges::ranges(gapit_results) <- GenomicRanges::ranges(MatrixGenerics::rowRanges(subset_vcf)[gapit_results$SNP])
		GenomeInfoDb::seqlevels(gapit_results) <- GenomeInfoDb::seqlevels(gwas_subset)

		if(length(gapit_results)) gwas_subset <- c(gwas_subset, gapit_results)
	}
}

# Saving the signals to an RDS file for retrieval later on
saveRDS(gwas_signals,
	file = paste0("gwas_results/", program, "/", trait, "_signal.rds"),
	compress = FALSE)

saveRDS(gwas_subset,
	file = paste0("gwas_results/", program, "/", trait, "_gwas_subset.rds"),
	compress = FALSE)

