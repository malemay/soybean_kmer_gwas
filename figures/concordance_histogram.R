# A histogram with the concordance values from the WGS/SoySNP50K comparison

# Reading the bcftools gtcheck results
# DEPENDENCY: illumina_data/soysnp50k_genotyping/gtcheck_results.tsv
gtcheck <- read.table("illumina_data/soysnp50k_genotyping/gtcheck_results.tsv", sep = "\t", header = TRUE)

png("figures/concordance_histogram.png", width = 6, height = 6, units = "in", res = 200)

hist(gtcheck$conc_rate, breaks = 30, main = NULL, xlab = "Concordance rate", ylab = "Number of accessions")
abline(v = 0.9, lty = 2)

dev.off()
