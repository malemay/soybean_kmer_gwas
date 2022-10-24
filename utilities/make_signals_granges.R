# This script combines the reference signals from various sources into a GRanges object
# that can be used to extract signals found by GWAS analyses and plot them

suppressMessages(library(GenomicRanges))

# DEPENDENCY: reference_signals/bandillo2015_signals.rds
# DEPENDENCY: reference_signals/bandillo2017_signals.rds
# DEPENDENCY: reference_signals/deronne2022_signal.rds
bandillo2015_signals <- readRDS("reference_signals/bandillo2015_signals.rds")
bandillo2017_signals <- readRDS("reference_signals/bandillo2017_signals.rds")
deronne_signal <- readRDS("reference_signals/deronne2022_signal.rds")
custom_signals <- readRDS("reference_signals/custom_signals.rds")

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
       signals = c(deronne_signal, bandillo2015_signals, bandillo2017_signals, custom_signals))

all_signals <- do.call("c", all_signals)
stopifnot(!anyDuplicated(all_signals$signal_id))

# Saving the list of signal names to file
writeLines(paste(names(all_signals), all_signals$trait, sep = ","), con = "utilities/signal_ids.txt")

# Also saving a list of signals for which the causal gene is known
gene_signals <- all_signals[!is.na(all_signals$gene_name_v4)]
writeLines(names(gene_signals), con = "utilities/gene_signal_ids.txt")

# Saving the GRanges object to an .rds file
saveRDS(all_signals, file = "utilities/all_signals.rds", compress = FALSE)

