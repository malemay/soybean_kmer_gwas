# Reading the trait for which the scaling is to be done from the command line
trait <- commandArgs(trailingOnly = TRUE)[1]

# Reading the threshold from the file
threshold <- as.numeric(readLines(paste0("gwas_results/kmers/", trait, "_threshold_5per")))

# Rescaling the threshold back to its true value (not -log10(p))
threshold <- 10^(-threshold)

# Writing it back to a file with a different name
writeLines(as.character(threshold), con = paste0("gwas_results/kmers/", trait, "_threshold_5per.txt"))

