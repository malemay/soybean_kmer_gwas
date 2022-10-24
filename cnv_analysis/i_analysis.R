# This scripts generates a GRanges object encompassing the
# tandem-inverted I locus based on its alignment with bwa

# Loading the required libraries
suppressMessages(library(Rsamtools))
suppressMessages(library(GenomicRanges))

# Loading the aligned I locus
# DEPENDENCY: cnv_analysis/i_locus.bam
i_locus <- scanBam("cnv_analysis/i_locus.bam")[[1]]

i_ranges <- GRanges(seqnames = as.character(i_locus$rname),
		    ranges = IRanges(start = i_locus$pos, width = i_locus$qwidth))

# We must ensure that this is a single range
stopifnot(length(i_ranges) == 1)

# Saving it to file
saveRDS(i_ranges, file = "cnv_analysis/i_cnv_range.rds", compress = FALSE)

