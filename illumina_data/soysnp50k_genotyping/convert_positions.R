# This script converts the genomic coordinates of the SoySNP50K chip from
# genome assembly version 2 to genome assembly version 4

# Loading the required libraries
library(parallel)
library(Rsamtools)
library(GenomicRanges)
library(Biostrings)

# Setting the working directory
setwd("illumina_data/soysnp50k_genotyping/")

# The padding (number of nucleotides to be extracted from either side of the position) to use as a parameter
padding <- 20

# The path to reference genome version 2
# DEPENDENCY: refgenome/Gmax_nuclv2_mit_chlp.fasta
gmax_v2_path <- "../../refgenome/Gmax_nuclv2_mit_chlp.fasta"

# The path to reference genome version 4
# DEPENDENCY: refgenome/Gmax_508_v4.0_mit_chlp.fasta
gmax_v4_path <- "../../refgenome/Gmax_508_v4.0_mit_chlp.fasta"

# Reading the positions of the SoySNP50K file
# DEPENDENCY: illumina_data/soysnp50k_genotyping/soysnp50K_maf10.vcf
soysnp <- read.table("soysnp50K_maf10.vcf", sep = "\t", header = FALSE)
names(soysnp)[1:3] <- c("chrom", "pos", "id")
rownames(soysnp) <- soysnp$id

# Getting the header and outputting it to file separately
soysnp_header <- grep("^#", readLines("soysnp50K_maf10.vcf"), value = TRUE)
writeLines(soysnp_header, con = "soysnp50K_gmax_v4.vcf")


# Creating a GRanges object to store the query positions
snp_ranges <- GRanges(seqnames = soysnp$chr, 
		      ranges = IRanges(start = soysnp$pos - padding,
				       end = soysnp$pos + padding))

# Getting the sequences
# we neglect the reverse sequence to speed up the query
# as we do not need to translate all possible positions anyway
soysnp$seq <- as.character(scanFa(gmax_v2_path, snp_ranges))

# Reading the Gmax version 4 sequence
gmax_v4 <- Biostrings::readDNAStringSet(gmax_v4_path)

# Finding the matches of sequences in Gmax_v4
matches <- mclapply(1:nrow(soysnp),
		    FUN = function(x, snp_data) {
			    chr <- snp_data[x, "chrom"]
			    output <- vmatchPattern(snp_data[x, "seq"], gmax_v4[sub("Chr", "Gm", chr)])
		    },
		    snp_data = soysnp,
		    mc.cores = 30)

# Setting the names to the id so they can be matched later on
names(matches) <- soysnp$id

# Keeping only SNPs for which a single match in genome version 4 was found
matches <- matches[lengths(matches) == 1]

# Gathering all the matches in a single data.frame
matches <- do.call("rbind", lapply(matches, function(x) as.data.frame(unlist(x))))

# Formatting the data
matches$id <- rownames(matches)
matches$pos <- matches$start + padding
matches <- matches[, c("names", "pos", "id")]
colnames(matches)[1] <- "chrom"

# Setting the positions to the original VCF file to those on Gmax version 4
soysnp$chrom <- NA
soysnp$pos <- NA

soysnp[, "chrom"] <- matches[rownames(soysnp), "chrom"]
soysnp[, "pos"] <- matches[rownames(soysnp), "pos"]

soysnp <- soysnp[complete.cases(soysnp[, c("chrom", "pos")]), ]

# Saving this new VCF file
# OUTPUT: illumina_data/soysnp50k_genotyping/soysnp50K_gmax_v4.vcf
write.table(soysnp, file = "soysnp50K_gmax_v4.vcf", quote = FALSE, sep = "\t",
	    row.names = FALSE, col.names = FALSE, append = TRUE)

