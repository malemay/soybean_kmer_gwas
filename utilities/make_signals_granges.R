# This script combines the reference signals from various sources into a GRanges object
# that can be used to extract signals found by GWAS analyses and plot them

library(GenomicRanges)

# DEPENDENCY: reference_signals/bandillo2015_signals.rds
# DEPENDENCY: reference_signals/bandillo2017_signals.rds
bandillo2015_signals <- readRDS("reference_signals/bandillo2015_signals.rds")
bandillo2017_signals <- readRDS("reference_signals/bandillo2017_signals.rds")

# Manually entering data for the signal found by de Ronne et al. (2021) with BLINK on Gm15
deronne_signal <- GRanges(seqnames = "Gm15",
			  ranges = IRanges(start = 36764744, end = 36764744))
deronne_signal$signal_id <- "cdwGm15"
deronne_signal$trait <- "corrected_dry_weight"
deronne_signal$locus <- "cdwGm15"
deronne_signal$n_snps <- 1
deronne_signal$log_pvalue <- -log10(7.1 * 10^-13)
deronne_signal$gene_name_v4 <- NA
deronne_signal$common_name <- NA
names(deronne_signal) <- deronne_signal$signal_id

# We need to duplicate every entry with the various coded traits that were used in analyses
trait_names <- readLines("utilities/trait_names.txt")

all_signals <- lapply(trait_names,
       FUN = function(x, signals) {
	       signals <- signals[substring(x, 1, nchar(signals$trait)) == signals$trait]
	       if(length(signals)) {
		       signals$original_trait <- signals$trait
		       signals$trait <- x
		       signals$signal_id <- paste0(signals$trait, "_", signals$locus)
		       names(signals) <- signals$signal_id
	       }
		       signals
       },
       signals = c(deronne_signal, bandillo2015_signals, bandillo2017_signals))

all_signals <- do.call("c", all_signals)
stopifnot(!anyDuplicated(all_signals$signal_id))

# Saving the list of signal names to file
writeLines(paste(names(all_signals), all_signals$trait, sep = ","), con = "utilities/signal_ids.txt")

# Also saving a list of signals for which the causal gene is known
gene_signals <- all_signals[!is.na(all_signals$gene_name_v4)]
writeLines(names(gene_signals), con = "utilities/gene_signal_ids.txt")

# Saving the GRanges object to an .rds file
saveRDS(all_signals, file = "utilities/all_signals.rds", compress = FALSE)

