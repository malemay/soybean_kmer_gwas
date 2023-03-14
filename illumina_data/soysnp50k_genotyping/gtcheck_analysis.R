# This script computes the concordance between the genotypes of
# samples as assessed from the WGS Illumina data and from the
# SoySNP50K chip

# Setting the working directory
setwd("illumina_data/soysnp50k_genotyping")

# Initializing a list to store the results into
# DEPENDENCY: illumina_data/soysnp50k_genotyping/pi_ids.txt
pi_ids <- read.table("pi_ids.txt")

# DEPENDENCY: results from running bcftools gtcheck on the WGS and SoySNP50K data
# Reading the gtcheck results
gtcheck_results <- lapply(1:nrow(pi_ids), function(x, pi_ids) {
				  bayer <- pi_ids[x, 1]
				  id <- pi_ids[x, 2]

				  # Reading the table from file
				  gtcheck <- read.table(paste0("gtcheck/", bayer, "_gtcheck.txt"), stringsAsFactors = FALSE)

				  # Safety check
				  stopifnot(nrow(gtcheck) == 385)

				  # Extracting only the columns we need and naming them
				  gtcheck <- gtcheck[, c(2, 4, 5)]
				  names(gtcheck) <- c("discordance", "n_comparisons", "pi_id")
				  gtcheck$bayer <- bayer

				  # Adding columns for concordance rate
				  gtcheck$conc_rate <- 1 - (gtcheck$discordance / gtcheck$n_comparisons)

				  # Sorting based on concordance
				  gtcheck$rank <- rank(-gtcheck$conc_rate)

				  # Getting the matching sample and returning it as the only row
				  match_index <- which(gtcheck$pi_id == id)
				  stopifnot(length(match_index) == 1)

				  gtcheck[match_index, c("bayer", "pi_id", "rank", "conc_rate", "discordance", "n_comparisons")]},

				  pi_ids = pi_ids)

# Binding all the data.frames together and reformatting the CAD_ID column
gtcheck_results <- do.call("rbind", gtcheck_results)
gtcheck_results <- gtcheck_results[order(gtcheck_results$conc_rate, decreasing = TRUE), ]

# Saving to file for retrieval later on
# OUTPUT: illumina_data/soysnp50k_genotyping/gtcheck_results.tsv
write.table(gtcheck_results, file = "gtcheck_results.tsv", sep = "\t",
	    col.names = TRUE, row.names = FALSE, quote = FALSE)

