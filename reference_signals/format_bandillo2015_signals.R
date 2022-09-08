# This script formats the signals found for oil and protein by Bandillo et al. (2015)
# as GRanges objects with coordinates relative to Williams82 genome assembly version 4

# Loading required libraries
library(parallel)
library(GenomicRanges)
library(Rsamtools)
library(Biostrings)

# The padding (number of nucleotides to be extracted from either side of the position) to use as a parameter
padding <- 20

# --- GROUPING ALL EXTERNAL DEPENDENCIES AT THE TOP OF THE FILE

# The path to reference genome version 1
# DEPENDENCY: refgenome/glyma.Wm82.gnm1.FCtY.genome_main.fna
gmax_v1_path <- "refgenome/glyma.Wm82.gnm1.FCtY.genome_main.fna"

# The path to reference genome version 4
# DEPENDENCY: refgenome/Gmax_508_v4.0_mit_chlp.fasta
gmax_v4_path <- "refgenome/Gmax_508_v4.0_mit_chlp.fasta"

# A curated table as obtained from the table in Bandillo et al. (2015) processed with Tabula
# DOI: 10.3835/plantgenome2015.04.0024
# DEPENDENCY: reference_signals/bandillo2015_table1_curated.csv
bandillo2015_snps <- read.table("reference_signals/bandillo2015_table1_curated.csv",
				sep = ",", header = TRUE,
				na.strings = "ns", stringsAsFactors = FALSE)

# Transforming the dataset from a "wide" dataset to a "long" dataset
oil_snps <- protein_snps <- bandillo2015_snps

oil_snps$trait <- "oil"
oil_snps$protein_logp <- oil_snps$protein_allelic_effect <- NULL
oil_snps <- oil_snps[, c("trait", "SNP", "Chromosome", "oil_logp", "oil_allelic_effect")]
colnames(oil_snps) <- c("trait", "SNP", "Chr", "logp", "allelic_effect")

protein_snps$trait <- "protein"
protein_snps$oil_logp <- protein_snps$oil_allelic_effect <- NULL
protein_snps <- protein_snps[, c("trait", "SNP", "Chromosome", "protein_logp", "protein_allelic_effect")]
colnames(protein_snps) <- c("trait", "SNP", "Chr", "logp", "allelic_effect")

bandillo2015_snps <- rbind(oil_snps, protein_snps)

# Keeping only associations of interest on chromsomes 15 and 20, and removing those with missing data
bandillo2015_snps <- bandillo2015_snps[complete.cases(bandillo2015_snps) & bandillo2015_snps$Chr %in% c(15, 20),]

# Formatting a few columns
bandillo2015_snps$Chr <- paste0("glyma.Wm82.gnm1.Gm", bandillo2015_snps$Chr)
bandillo2015_snps$signal_id <- paste0(bandillo2015_snps$trait, "_", bandillo2015_snps$Chr)
bandillo2015_snps$pos <- as.numeric(sub(".*_([0-9]{4,})_.*", "\\1", bandillo2015_snps$SNP))

bandillo2015_signals <- data.frame(signal_id = unique(bandillo2015_snps$signal_id),
				    stringsAsFactors = FALSE)
bandillo2015_signals$trait <- sapply(strsplit(bandillo2015_signals$signal_id, "_"), function(x) x[1])
bandillo2015_signals$Chr <- sapply(strsplit(bandillo2015_signals$signal_id, "_"), function(x) x[2])
bandillo2015_signals$n_snps <- as.numeric(table(bandillo2015_snps$signal_id)[bandillo2015_signals$signal_id])
bandillo2015_signals$first_snp <- tapply(bandillo2015_snps$pos, bandillo2015_snps$signal_id, min)[bandillo2015_signals$signal_id]
bandillo2015_signals$last_snp <- tapply(bandillo2015_snps$pos, bandillo2015_snps$signal_id, max)[bandillo2015_signals$signal_id]
bandillo2015_signals$max_snp <- sapply(split(bandillo2015_snps[, c("pos", "logp")], bandillo2015_snps$signal_id), function(x) x$pos[which.max(x$logp)])[bandillo2015_signals$signal_id]
bandillo2015_signals$log_pvalue <- tapply(bandillo2015_snps$logp, bandillo2015_snps$signal_id, max)[bandillo2015_signals$signal_id]
bandillo2015_signals$locus <- bandillo2015_signals$signal_id

bandillo2015_signals$glyma_name <- NA
bandillo2015_signals$common_name <- NA
bandillo2015_signals$gene_name_v4 <- NA

# Setting the name of the gene for the QTL on chromosome 15
# according to Zhang et al. (2020), DOI:10.1371/journal.pgen.1009114
bandillo2015_signals[bandillo2015_signals$Chr == "glyma.Wm82.gnm1.Gm15", "gene_name_v4"] <- tolower("Glyma.15G049200")

# Reordering the columns
bandillo2015_signals <- bandillo2015_signals[, c("trait", "Chr", "n_snps",
						 "first_snp", "last_snp", "max_snp",
						 "log_pvalue", "locus", "glyma_name",
						 "common_name", "gene_name_v4", "signal_id")]

# Generating the fai index for reference genome version 1 if it does not already exist
if(!file.exists(paste0(gmax_v1_path, ".fai"))) Rsamtools::indexFa(gmax_v1_path)

# Doing the same for reference genome version 4
if(!file.exists(paste0(gmax_v4_path, ".fai"))) Rsamtools::indexFa(gmax_v4_path)

# Creating GRanges objects to store the query positions
first_granges <- GRanges(seqnames = bandillo2015_signals$Chr, 
			 ranges = IRanges(start = bandillo2015_signals$first_snp - padding,
					  end = bandillo2015_signals$first_snp + padding))

last_granges <- GRanges(seqnames = bandillo2015_signals$Chr, 
			ranges = IRanges(start = bandillo2015_signals$last_snp - padding,
					 end = bandillo2015_signals$last_snp + padding))

max_granges <- GRanges(seqnames = bandillo2015_signals$Chr, 
		       ranges = IRanges(start = bandillo2015_signals$max_snp - padding,
					  end = bandillo2015_signals$max_snp + padding))

# Getting the sequences
bandillo2015_signals$first_seq <- as.character(scanFa(gmax_v1_path, first_granges))
bandillo2015_signals$first_rev <- as.character(reverseComplement(scanFa(gmax_v1_path, first_granges)))
bandillo2015_signals$last_seq <- as.character(scanFa(gmax_v1_path, last_granges))
bandillo2015_signals$last_rev <- as.character(reverseComplement(scanFa(gmax_v1_path, last_granges)))
bandillo2015_signals$max_seq <- as.character(scanFa(gmax_v1_path, max_granges))
bandillo2015_signals$max_rev <- as.character(reverseComplement(scanFa(gmax_v1_path, max_granges)))

# Reading the Gmax version 4 sequence
gmax_v4 <- Biostrings::readDNAStringSet(gmax_v4_path)

for(i in c("first", "last", "max")) {
	message("Finding the forward matches for ", i)
	forward_matches <- mclapply(bandillo2015_signals$signal_id,
				    FUN = function(x) vmatchPattern(bandillo2015_signals[bandillo2015_signals$signal_id == x, paste0(i, "_seq")], gmax_v4),
				    mc.cores = 4)

	forward_matches <- lapply(forward_matches, unlist)

	message("Finding the reverse matches for ", i)
	reverse_matches <- mclapply(bandillo2015_signals$signal_id,
				    function(x) vmatchPattern(bandillo2015_signals[bandillo2015_signals$signal_id == x, paste0(i, "_rev")], gmax_v4),
				    mc.cores = 4)

	reverse_matches <- lapply(reverse_matches, unlist)

	# Checking that either forward or reverse matched, but not both; also that all positions had a single match
	if(!all((lengths(forward_matches) + lengths(reverse_matches)) == 1)) {
		warning("At least one query sequence for which no single match was found")
	} else {
		message("Only one match found per query sequence")
	}

	# Getting a data.frame from all the matches
	all_matches <- ifelse(lengths(forward_matches), forward_matches, reverse_matches)
	all_matches <- ifelse(lengths(all_matches) == 1, all_matches, IRanges(start = 1, end = 1, names = "NOMATCH"))
	all_matches <- do.call("rbind", lapply(all_matches, as.data.frame))

	# We only need the chromosome and the position
	all_matches <- all_matches[, c(4, 1)]
	names(all_matches) <- c(paste0(i, "_chrom"), paste0(i, "_pos"))

	# The actual position of the SNP is +padding +1
	all_matches[[2]] <- all_matches[[2]] + padding + 1

	# Adding those columns to the main data
	bandillo2015_signals <- cbind(bandillo2015_signals, all_matches)
}


# Checking that the plot of the widths inferred from both versions still looks the same
#with(bandillo2015_signals, plot(x = abs(first_snp - last_snp), y = abs(first_pos - last_pos)))
#abline(0, 1)

# Let us check that the max_snp is always between first_snp and last_snp
all(bandillo2015_signals$first_pos <= bandillo2015_signals$max_pos)
all(bandillo2015_signals$max_pos <= bandillo2015_signals$last_pos)

# Everything seems fine with the new coordinates on Gmax v4

# Now let us prepare the dataset for saving it as a GenomicRanges object
# We first remove the columns that we no longer need and re-order the columns
bandillo2015_signals <- bandillo2015_signals[, c("signal_id", "trait", "locus", "n_snps",
						 "log_pvalue", "gene_name_v4", "common_name",
						 "first_chrom", "first_pos", "last_pos")]

# Reformatting the signal_id and locus columns
bandillo2015_signals$signal_id <- bandillo2015_signals$locus <- paste0(bandillo2015_signals$trait, bandillo2015_signals$first_chrom)

bandillo2015_signals <- makeGRangesFromDataFrame(bandillo2015_signals,
						 keep.extra.columns = TRUE,
						 ignore.strand = TRUE,
						 seqnames.field = "first_chrom",
						 start.field = "first_pos",
						 end.field = "last_pos")

names(bandillo2015_signals) <- bandillo2015_signals$signal_id

# Saving the output to file
saveRDS(bandillo2015_signals, file = "reference_signals/bandillo2015_signals.rds")

