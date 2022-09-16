# This script creates the signal found by de Ronne et al. (2022)
# as a GRanges object for use in downstream analyses
# DOI:10.1002/tpg2.20184

# Loading the GenomicRanges package
library(GenomicRanges)

# Manually entering data for the signal found by de Ronne et al. (2021) with BLINK on Gm15
deronne_signal <- GRanges(seqnames = "Gm15",
			  ranges = IRanges(start = 36000000, end = 40000000))
deronne_signal$signal_id <- "cdwGm15"
deronne_signal$trait <- "corrected_dry_weight"
deronne_signal$locus <- "cdwGm15"
deronne_signal$n_snps <- 31
deronne_signal$log_pvalue <- -log10(7.1 * 10^-13)
deronne_signal$gene_name_v4 <- tolower("Glyma.15G217100")
deronne_signal$common_name <- NA
names(deronne_signal) <- deronne_signal$signal_id

# Saving this signal to file for retrieval in downstream analyses
saveRDS(deronne_signal, "reference_signals/deronne2022_signal.rds")

