# This script translates the positions of the signals found by Bandillo et al. (2017)
# from reference genome version 1 to reference genome version 4
#
# This will imply retrieving the positions of the flanking SNPs and maximal SNPs
# in their dataset, extracting the flanking sequences in reference genome version 1
# and look for matches in reference genome version 4
#
# As part of this script we also set the names of the genes that are associated
# with various loci based on literature. The references to the literature can be
# found in the comments.

# Loading the required libraries
library(parallel)
library(Rsamtools)
library(GenomicRanges)
library(Biostrings)

# The padding (number of nucleotides to be extracted from either side of the position) to use as a parameter
padding <- 20

# --- GATHERING ALL EXTERNAL DEPENDENCIES AT THE TOP OF THE FILE

# The path to reference genome version 1
# DEPENDENCY: refgenome/glyma.Wm82.gnm1.FCtY.genome_main.fna
gmax_v1_path <- "refgenome/glyma.Wm82.gnm1.FCtY.genome_main.fna"

# The path to reference genome version 4
# DEPENDENCY: refgenome/Gmax_508_v4.0_mit_chlp.fasta
gmax_v4_path <- "refgenome/Gmax_508_v4.0_mit_chlp.fasta"

# The signals detected by Bandillo et al. (2017), as obtained from the article PDF with Tabula
# DOI: 10.3835/plantgenome2016.06.0054
# DEPENDENCY: reference_signals/bandillo2017_signals_curated.tsv
bandillo_signals <- read.table("reference_signals/bandillo2017_signals_curated.tsv",
			       sep = "\t", header = TRUE,
			       na.strings = "-", stringsAsFactors = FALSE)

# Loading a lookup table of gene names from Williams82 genome assembly version 1 to version 4
# DEPENDENCY: refgenome/lookup_v1_to_v4.rds
gene_name_lookup <- readRDS("refgenome/lookup_v1_to_v4.rds")

# Getting the names of the genes that are Chalcone synthase genes on chromosome Gm08 in version 4
# DEPENDENCY: refgenome/soybase_genome_annotation_v4.0_04-20-2021.txt
chalcone_genes <- system("grep Glyma.08g.*CHS refgenome/soybase_genome_annotation_v4.0_04-20-2021.txt | awk '{print $1}'",
			 intern = TRUE)

# Getting the names of the genes that are HPS genes on chromosome 15
hps_genes <- system("grep Glyma.15g.*Hydrophob_seed refgenome/soybase_genome_annotation_v4.0_04-20-2021.txt | awk '{print $1}'",
		    intern = TRUE)

# Getting a GRanges object of the annotated genes in Williams82 assembly version 4
# DEPENDENCY: refgenome/gmax_v4_genes.rds
gmax_v4_genes <- readRDS("refgenome/gmax_v4_genes.rds")

# Keeping only the columns that are needed for downstream analyses
bandillo_signals <- bandillo_signals[, c("trait", "Chr", "n_snps",
					 "first_snp", "last_snp", "max_snp",
					 "log_pvalue", "locus", "glyma_name",
					 "common_name")]

# Renaming the traits to standard names
trait_lookup <- c("maturity" = "maturity_group",
		  "growth_form" = "stem_termination",
		  "flower_color" = "flower_color",
		  "hair_color" = "pubescence_color",
		  "hair_form" = "pubescence_form",
		  "hair_density" = "pubescence_density",
		  "pod_color" = "pod_color",
		  "luster" = "seed_coat_luster",
		  "coat_color" = "seed_coat_color",
		  "hilum_color" = "hilum_color")

bandillo_signals$trait <- trait_lookup[bandillo_signals$trait]

# We are not interested in signals that do not correspond to named loci
bandillo_signals <- bandillo_signals[!is.na(bandillo_signals$locus), ]

# For processing purposes we rename the loci to only keep the name of the dominant allele
bandillo_signals$locus <- sub("-.*", "", bandillo_signals$locus)
bandillo_signals[bandillo_signals$locus == "B?", "locus"] <- "Bq"

# Renaming the chromosomes so that they match their name in the Gmax_v1 fasta file
bandillo_signals$Chr <- paste0("glyma.Wm82.gnm1.",
			       ifelse(bandillo_signals$Chr < 10, "Gm0", "Gm"),
			       as.character(bandillo_signals$Chr))

# Giving a unique ID to each signal, based on the trait and the locus name
bandillo_signals$signal_id <- paste(bandillo_signals$trait, bandillo_signals$locus, sep = "_")
stopifnot(!any(duplicated(bandillo_signals$signal_id)))

# Adding the name of the gene in the genome version 4
bandillo_signals$glyma_name <- tolower(bandillo_signals$glyma_name)
bandillo_signals$gene_name_v4 <- gene_name_lookup[bandillo_signals$glyma_name]

# Looking at cases where a gene is reported for Gmax_v1 but no match was found in Gmax_v4
bandillo_signals[!is.na(bandillo_signals$glyma_name) & is.na(bandillo_signals$gene_name_v4), ]

# This is the case for the Chalcone synthase genes on chromosome 8
chalcone_genes <- tolower(chalcone_genes)
stopifnot(all(chalcone_genes %in% names(gmax_v4_genes)))
chalcone_genes <- sort(gmax_v4_genes[chalcone_genes], ignore.strand = TRUE)
width(reduce(chalcone_genes, min.gapwidth = 10^5, ignore.strand = TRUE))
# [1] 142062

# Let us add the names of the first and last of those genes in the gene_name_v4 column
chalcone_gene_ids <- paste0(names(chalcone_genes)[c(1, length(chalcone_genes))], collapse = ";")
bandillo_signals[bandillo_signals$locus == "I", "gene_name_v4"] <- chalcone_gene_ids

# The other gene not found in Gmax v1 annotation is glyma19g27460 and its
# Gmax v4 equivalent is Glyma.19G101700 according to Sedivy et al. (2017)
# DOI:10.1111/nph.14418
# bandillo_signals[bandillo_signals$locus == "L1", "gene_name_v4"] <- tolower("Glyma.19G101700")
# this gene, however, does not appear to be the right one. We suggest Glyma.19G120300 instead
bandillo_signals[bandillo_signals$locus == "L1", "gene_name_v4"] <- tolower("Glyma.19G120300")

# Let us also add the genes that have been cloned after Bandillo et al. (2017) and
# therefore not originally listed in their work

# Td locus: Yan et al. (2020), DOI:10.3389/fpls.2019.01809
bandillo_signals[bandillo_signals$locus == "Td", "gene_name_v4"] <- tolower("Glyma.03G258700")

# G locus: Wang et al. (2018), DOI:10.1038/s41588-018-0229-2
bandillo_signals[bandillo_signals$locus == "G", "gene_name_v4"] <- tolower("Glyma.01G198500")

# Pa1 locus: Gilbert (2017), https://hdl.handle.net/11299/193416
bandillo_signals[bandillo_signals$locus == "Pa1", "gene_name_v4"] <- tolower("glyma.12g213900")

# Ps locus: Liu et al. (2020), DOI:10.1016/j.molp.2020.10.004
bandillo_signals[bandillo_signals$locus == "Ps", "gene_name_v4"] <- tolower("Glyma.12G187200")

# B locus: Gijzen et al. (2006), DOI:10.1186/1471-2229-6-6
# This is actually a cluster of HPS genes
# We have obtained the names of the genes corresponding to HPS genes on chromosome 15 above
# However we have to restrict that set to the genes that are located in a cluster around position 9-Mb - 10Mb
hps_genes <- tolower(hps_genes)
stopifnot(all(hps_genes %in% names(gmax_v4_genes)))
hps_genes <- sort(gmax_v4_genes[hps_genes], ignore.strand = TRUE)
hps_genes <- subsetByOverlaps(hps_genes,
			      GRanges(seqnames = "Gm15", ranges = IRanges(start = 9 * 10^6, end = 10 * 10^6)),
			      ignore.strand = TRUE)

# Let us add the names of the first and last of those genes in the gene_name_v4 column
hps_gene_ids <- paste0(names(hps_genes)[c(1, length(hps_genes))], collapse = ";")
bandillo_signals[bandillo_signals$locus == "B", "gene_name_v4"] <- hps_gene_ids

# Generating the fai index for reference genome version 1 if it does not already exist
if(!file.exists(paste0(gmax_v1_path, ".fai"))) Rsamtools::indexFa(gmax_v1_path)

# Doing the same for reference genome version 4
if(!file.exists(paste0(gmax_v4_path, ".fai"))) Rsamtools::indexFa(gmax_v4_path)

# Creating a GRanges object to store the query positions
first_granges <- GRanges(seqnames = bandillo_signals$Chr, 
			 ranges = IRanges(start = bandillo_signals$first_snp - padding,
					  end = bandillo_signals$first_snp + padding))

last_granges <- GRanges(seqnames = bandillo_signals$Chr, 
			ranges = IRanges(start = bandillo_signals$last_snp - padding,
					 end = bandillo_signals$last_snp + padding))

max_granges <- GRanges(seqnames = bandillo_signals$Chr, 
		       ranges = IRanges(start = bandillo_signals$max_snp - padding,
					  end = bandillo_signals$max_snp + padding))

# Getting the sequences
bandillo_signals$first_seq <- as.character(scanFa(gmax_v1_path, first_granges))
bandillo_signals$first_rev <- as.character(reverseComplement(scanFa(gmax_v1_path, first_granges)))
bandillo_signals$last_seq <- as.character(scanFa(gmax_v1_path, last_granges))
bandillo_signals$last_rev <- as.character(reverseComplement(scanFa(gmax_v1_path, last_granges)))
bandillo_signals$max_seq <- as.character(scanFa(gmax_v1_path, max_granges))
bandillo_signals$max_rev <- as.character(reverseComplement(scanFa(gmax_v1_path, max_granges)))

# Reading the Gmax version 4 sequence
gmax_v4 <- Biostrings::readDNAStringSet(gmax_v4_path)

# Finding the matches of sequences in Gmax_v4

for(i in c("first", "last", "max")) {
	message("Finding the forward matches for ", i)
	forward_matches <- mclapply(bandillo_signals$signal_id,
				    FUN = function(x) vmatchPattern(bandillo_signals[bandillo_signals$signal_id == x, paste0(i, "_seq")], gmax_v4),
				    mc.cores = 33)

	forward_matches <- lapply(forward_matches, unlist)

	message("Finding the reverse matches for ", i)
	reverse_matches <- mclapply(bandillo_signals$signal_id,
				    function(x) vmatchPattern(bandillo_signals[bandillo_signals$signal_id == x, paste0(i, "_rev")], gmax_v4),
				    mc.cores = 33)

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
	bandillo_signals <- cbind(bandillo_signals, all_matches)
}

# Checking a few things about the output
# First, have the same chromosomes been identified for the three positions
bandillo_signals[bandillo_signals$first_chrom != bandillo_signals$last_chrom, ]

# In these two cases, one of the matches appears to have failed
# Let us correct these cases
bandillo_signals[bandillo_signals$signal_id == "seed_coat_luster_B", "last_chrom"] <- NA
bandillo_signals[bandillo_signals$signal_id == "seed_coat_luster_B", "last_pos"] <- NA

bandillo_signals[bandillo_signals$signal_id == "seed_coat_color_T", "first_chrom"] <- NA
bandillo_signals[bandillo_signals$signal_id == "seed_coat_color_T", "first_pos"] <- NA

# Checking whether everytyhing is okay now
stopifnot(sum(bandillo_signals$first_chrom != bandillo_signals$last_chrom, na.rm = TRUE) == 0)
stopifnot(sum(bandillo_signals$first_chrom != bandillo_signals$max_chrom, na.rm = TRUE) == 0)
stopifnot(sum(bandillo_signals$last_chrom != bandillo_signals$max_chrom, na.rm = TRUE) == 0)

# The max_chrom column does not have any NA so we can use it to identify rows where the
# Gmax_v1 chromosome does not match the one on Gmax_v4
stopifnot(!any(is.na(bandillo_signals$max_chrom)))
bandillo_signals[sub(".*(Gm[0-9]{2}$)", "\\1", bandillo_signals$Chr) != bandillo_signals$max_chrom, ]

# Now we should take care of filling in the values for which there were NAs introduced
bandillo_signals[!complete.cases(bandillo_signals[, c("first_chrom", "first_pos", 
						      "last_chrom", "last_pos",
						      "max_chrom", "max_pos")]), ]

# For the first signal, we can use the width of the original signal to estimate the width and therefore position of the last SNP
# This signal is luster_glyma.Wm82.gnm1.Gm15_10416352
index1 <- which(bandillo_signals$signal_id == "seed_coat_luster_B")
direction <- sign(bandillo_signals[index1, "max_pos"] - bandillo_signals[index1, "first_pos"])
bandillo_signals[index1, "last_pos"] <- 
	bandillo_signals[index1, "first_pos"] +
	direction * (bandillo_signals[index1, "last_snp"] - bandillo_signals[index1, "first_snp"])
bandillo_signals[index1, "last_chrom"] <- bandillo_signals[index1, "first_chrom"]

# For the second signal with NAs, we can similarly estimate the position of first_pos
# This signal is coat_color_glyma.Wm82.gnm1.Gm06_18766611
index2 <- which(bandillo_signals$signal_id == "seed_coat_color_T")
direction <- sign(bandillo_signals[index2, "max_pos"] - bandillo_signals[index2, "last_pos"])
bandillo_signals[index2, "first_pos"] <- 
	bandillo_signals[index2, "last_pos"] +
	direction * (bandillo_signals[index2, "last_snp"] - bandillo_signals[index2, "first_snp"])
bandillo_signals[index2, "first_chrom"] <- bandillo_signals[index2, "last_chrom"]

# Now there shouldn't remain any data points with incomplete data
stopifnot(nrow(bandillo_signals[!complete.cases(bandillo_signals[, c("first_chrom", "first_pos", 
								      "last_chrom", "last_pos",
								      "max_chrom", "max_pos")]), ]) == 0)

# Checking whether all intervals are still in the same orientation
table((bandillo_signals$last_pos - bandillo_signals$first_pos) >= 0)
# 
# FALSE  TRUE 
#     2    31 

# Let us look at the ones that are inverted relative to the original values
bandillo_signals[(bandillo_signals$last_pos - bandillo_signals$first_pos) < 0, ]

# They are all unsurprisingly located on Gm13 as parts of this chromosome were inverted between genome assemblies
# For these we will invert the start and end positions
invert_indices <- which((bandillo_signals$last_pos - bandillo_signals$first_pos) < 0)

for(i in invert_indices) {
	tmp <- bandillo_signals[i, "first_pos"]
	bandillo_signals[i, "first_pos"] <- bandillo_signals[i, "last_pos"]
	bandillo_signals[i, "last_pos"] <- tmp
}

# Re-checking the orientations
table((bandillo_signals$last_pos - bandillo_signals$first_pos) >= 0)
# 
# TRUE 
#   33 

# Now let us prepare the dataset for saving it as a GenomicRanges object
# We first remove the columns that we no longer need and re-order the columns
bandillo_signals <- bandillo_signals[, c("signal_id", "trait", "locus", "n_snps",
					 "log_pvalue", "gene_name_v4", "common_name",
					 "first_chrom", "first_pos", "last_pos")]

bandillo_signals <- makeGRangesFromDataFrame(bandillo_signals,
					     keep.extra.columns = TRUE,
					     ignore.strand = TRUE,
					     seqnames.field = "first_chrom",
					     start.field = "first_pos",
					     end.field = "last_pos")

names(bandillo_signals) <- bandillo_signals$signal_id

# Saving the output to file
saveRDS(bandillo_signals, file = "reference_signals/bandillo2017_signals.rds")

