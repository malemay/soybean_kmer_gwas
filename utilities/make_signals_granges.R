# This script combines the reference signals from various sources into a GRanges object
# that can be used to extract signals found by GWAS analyses and plot them

library(GenomicRanges)

bandillo_signals <- readRDS("utilities/bandillo_signals.rds")

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
       signals = bandillo_signals)

all_signals <- do.call("c", all_signals)
stopifnot(!anyDuplicated(all_signals$signal_id))

# Saving the list of signal names to file
writeLines(paste(names(all_signals), all_signals$trait, sep = ","), con = "utilities/signal_ids.txt")

# Saving the GRanges object to an .rds file
saveRDS(all_signals, file = "utilities/all_signals.rds", compress = FALSE)

