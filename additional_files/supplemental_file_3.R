# Supplemental file with metadata on the mapped reads for each sample
library(parallel)

# Reading the Bayer IDs and their corresponding SRR IDs
# DEPENDENCY: utilities/srr_id_correspondence.txt
bayer_ids <- read.table("utilities/srr_id_correspondence.txt", header = FALSE, stringsAsFactors = FALSE)
colnames(bayer_ids) <- c("study_ID", "SRA_runs")

# DEPENDENCY: stats obtained from samtools stats for all samples
# Now we read the metadata obtained from the samtools stats command to add the number of reads
# and number of aligned reads
samtools_stats <- mclapply(bayer_ids$study_ID, function(x) {stats <- readLines(paste0("illumina_data/merged_bams/", x, "_stats.txt"))
			   stats <- grep("^SN", stats, value = TRUE)
			   nreads <- as.numeric(strsplit(grep("raw total sequences:", stats, value = TRUE), "\t")[[1]][3])
			   nmapped <- as.numeric(strsplit(grep("reads mapped:", stats, value = TRUE), "\t")[[1]][3])
			   list(nreads = nreads, nmapped = nmapped)},
			   mc.cores = 10)

names(samtools_stats) <- bayer_ids$study_ID

# Adding these data to the data.frame
bayer_ids$nreads <- sapply(samtools_stats, function(x) x$nreads)[bayer_ids$study_ID]
bayer_ids$nmapped <- sapply(samtools_stats, function(x) x$nmapped)[bayer_ids$study_ID]

# Reading the mapping coverage from the manifest files
# DEPENDENCY: manifest files for all samples
average_cov <- sapply(bayer_ids$study_ID, function(x) {
			      read.table(paste0("sv_genotyping/paragraph/manifest_files/", x, "_manifest.txt"),
					 skip = 1)[[3]]})

bayer_ids$average_mapping_depth <- average_cov[bayer_ids$study_ID]

# Adding the results on the concordance rate between the samples and
# DEPENDENCY: illumina_data/soysnp50k_genotyping/gtcheck_results.tsv
concordance_rates <- read.table("illumina_data/soysnp50k_genotyping/gtcheck_results.tsv", sep = "\t",
				header = TRUE, stringsAsFactors = FALSE)

# Adding the concordance rates to the data.frame
bayer_ids$PI_id <- concordance_rates[match(bayer_ids$study_ID, concordance_rates$bayer), "pi_id"]
bayer_ids$soysnp50k_concordance <- concordance_rates[match(bayer_ids$study_ID, concordance_rates$bayer), "conc_rate"]

# Writing the table to a CSV file
write.table(bayer_ids, file = "additional_files/supplemental_file_3.csv", sep = ",",
	    col.names = TRUE, row.names = FALSE, quote = FALSE)
