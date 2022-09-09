# Loading the required libraries
library(grid)
library(Biostrings)
library(IRanges)
library(stringr)

max_kmers <- 10000
kmer_length <- 31

# A function that maps a set of numeric values onto a color palette
# Values at the minimum are set to black
map_color <- function(values, max, min = 0, pal = "YlOrRd", n.colors = 9) {
	# Breaking the numeric values into a factor
	index <- cut(values, breaks = seq(min, max, length.out = n.colors))
	palette <- RColorBrewer::brewer.pal(n = n.colors, name = pal)
	output <- palette[as.integer(index)]
	# Setting the values at the minimum to black
	output[values <= min] <- "black"
	output
}

# This script generates a k-mer plot from the consensus sequences of a group of samples
locus <- commandArgs(trailingOnly = TRUE)[1]
locus_data <- read.table("utilities/kmer_plot_ranges.txt")
trait <- locus_data[locus_data[[1]] == locus, 2]

# DEPENDENCY: phenotypic_data/phenotypic_data.csv
phenotypes <- read.table("phenotypic_data/phenotypic_data.csv",
			 sep = ";", header = TRUE,)

# DEPENDENCY: phenotypic_data/lookup_tables.rds
lookup_tables <- readRDS("phenotypic_data/lookup_tables.rds")
lookup_table  <- lookup_tables[[paste0(trait, "_lookup")]]

# DEPENDENCY: k-mer p-values
kmer_pvalues <- read.table(paste0("gwas_results/kmer_data/", trait, "/katcher_results/pass_threshold_5per_sorted.txt"),
			   colClasses = c("NULL", "character", rep("NULL", 6), "numeric"),
			   nrows = max_kmers)
names(kmer_pvalues) <- c("kmer_id", "pvalue")

kmer_pvalues$kmer <- sub("[^ATGC]+", "", kmer_pvalues$kmer_id)
kmer_pvalues$kmer_reverse <- as.character(reverseComplement(DNAStringSet(kmer_pvalues$kmer)))

# DEPENDENCY: consensus sequence for a given locus
fasta_file <- readLines(paste0("gwas_results/kmer_consensus/", locus, "_sequences.fa"))

# Using conensus sequences from bwa in case the assemblies did not yield great results
if(length(fasta_file) < 200) {
	warning("Using consensus sequences from read alignment by bwa for locus ", locus)
	fasta_file <- readLines(paste0("gwas_results/kmer_consensus/", locus, "_bwa_sequences.fa"))
}

stopifnot(length(fasta_file) > 200)

# Extracting the names of the samples from the fasta file
sample_names <- grep("^>", fasta_file, value = TRUE)
sample_names <- sub("^>", "", sample_names)

# Extracting the sequences themselves from the fasta file
sequences <- grep("^>", fasta_file, value = TRUE, invert = TRUE)
names(sequences) <- sample_names

# Getting the number of times that a given sequence occurs across the whole dataset
# I will call these haplotypes
haplotypes <- sort(table(sequences))

# Keeping only haplotypes that occur at least five times
haplotypes <- haplotypes[haplotypes >= 5]

# Creating a data.frame linking haplotypes to phenotypes
haplotype_data <- data.frame(sample = names(sequences),
			     haplotype = unname(sequences))

# Assigning a haplotype to each sample
haplotype_data$id <- NA
for(i in 1:nrow(haplotype_data)) {
	if(haplotype_data[i, "haplotype"] %in% names(haplotypes)) {
		haplotype_data[i, "id"] <- which(names(haplotypes) == haplotype_data[i, "haplotype"])
	}
}

# Adding the phenotype for each sample
stopifnot(all(haplotype_data$sample[haplotype_data$sample != "Williams82"] %in% phenotypes$bayer_id))
haplotype_data$phenotype_numeric <- phenotypes[match(haplotype_data$sample, phenotypes$bayer_id), trait]
haplotype_data$phenotype_character <- NA
for(i in 1:nrow(haplotype_data)) {
	i_pheno <- haplotype_data[i, "phenotype_numeric"]
	if(!is.na(i_pheno)) {
	haplotype_data[i, "phenotype_character"] <- names(lookup_table)[which(i_pheno == lookup_table)]
	}
}

# Finding the matches of the significant k-mers in the haplotypes
kmer_overlaps <- lapply(names(haplotypes), function(x, kmers) {
				kmers$fmatch <- sapply(kmer_pvalues$kmer, function(kmer) regexpr(kmer, x, fixed = TRUE)) 
				kmers$rmatch <- sapply(kmer_pvalues$kmer_reverse, function(kmer) regexpr(kmer, x, fixed = TRUE)) 
				kmers$matchpos <- pmax(kmers$fmatch, kmers$rmatch)
				kmers
			     },
			     kmers = kmer_pvalues)
names(kmer_overlaps) <- names(haplotypes)

# Keeping only the k-mers for which there is at least one overlapping k-mer with a p-value
# Also removing irrelevant columns
kmer_overlaps <- lapply(kmer_overlaps, function(x) x[x$matchpos != -1, c("pvalue", "matchpos")])

# Coercing to an IRanges object
kmer_overlaps <- lapply(kmer_overlaps, function(x) {
				IRanges(start = x$matchpos, width = kmer_length, log10p = -log10(x$pvalue))
			     })

# Generating a data.frame with separate nucleotides and their associated p-value for each haplotype
plotting_data <-
	lapply(names(haplotypes), function(x, kmer_overlaps) {
		       output_df <- data.frame(pos = 1:nchar(x),
					       nuc = strsplit(x, "")[[1]],
					       log10p = 0)

		       # Getting the maximum -log10(p-value) for that nucleotide
		       for(i in 1:nrow(output_df)) {
			       i_range <- IRanges(start = output_df[i, "pos"], width = 1)
			       overlapped_pos <- subsetByOverlaps(kmer_overlaps[[x]], i_range)
			       if(length(overlapped_pos)) {
				       output_df[i, "log10p"] <- max(mcols(overlapped_pos)$log10p)
			       }
		       }
		       
		       output_df
			     },
		       kmer_overlaps = kmer_overlaps)

# Determining the plotting color from a common palette for all haplotypes
max_pvalue <- max(do.call("rbind", plotting_data)$log10p)

plotting_data <- lapply(plotting_data, function(x) {
		       x$color <- map_color(values = x$log10p, min = 0, max = max_pvalue, n.colors = 9, pal = "YlOrRd")
		       x
		       })

# Performing multiple alignment of the sequences to find the gaps
output_fasta <- file(paste0("gwas_results/kmer_consensus/", locus, "_haplotypes.fa"), open = "w+")
on.exit(close(output_fasta))

for(i in 1:length(haplotypes)) {
	cat(paste0(">hap", i, "\n"), file = output_fasta)
	cat(names(haplotypes)[i], "\n", file = output_fasta)
}

alignment <- system(paste0("mafft --auto gwas_results/kmer_consensus/", locus, "_haplotypes.fa"), intern = TRUE)
alignment <- strsplit(paste0(alignment, collapse = ""), split = ">hap[0-9]+")[[1]]
alignment <- alignment[nchar(alignment) > 0]

# Adjusting the positions for plotting depending on the gaps in the alignment
gaps <- str_locate_all(alignment, "-")
gaps <- lapply(gaps, function(x) IRanges(start = x[, 1], end = x[, 2]))
gaps <- lapply(gaps, function(x) reduce(x))

# A function that takes a data.frame of plotting data and
# a set of gaps, and adjusts the plotting positions accordingly
adjust_gaps <- function(hapdata, gaps) {
	if(!length(gaps)) return(hapdata)

	for(i in 1:length(gaps)) {
		hapdata[hapdata$pos >= start(gaps[i]), "pos"] <- hapdata[hapdata$pos >= start(gaps[i]), "pos"] + width(gaps[i])
	}

	return(hapdata)
}

plotting_data <- mapply(FUN = adjust_gaps, hapdata = plotting_data, gaps = gaps, SIMPLIFY = FALSE)

# A function that fills the deleted positions with dashes
fill_gaps <- function(x) {
	# We need to get the set of positions covered from 1 to the maximum
	maxpos <- max(do.call("rbind", x)$pos)

	# Then we loop over all the data.frames and fill the missing positions with dashes
	for(i in 1:length(x)) {
		indices <- which(! 1:maxpos %in% x[[i]]$pos)
		if(!length(indices)) next
		new_rows <- data.frame(pos = indices,
				       nuc = "-",
				       log10p = 0,
				       color = "black",
				       stringsAsFactors = FALSE)
		x[[i]] <- rbind(x[[i]], new_rows)
		x[[i]] <- x[[i]][order(x[[i]]$pos), ]
	}

	x
}

plotting_data <- fill_gaps(plotting_data)

# A function that finds the positions where the nucleotides differ between two consecutive sequences
nucdiff <- function(hapdata) {
	if(length(hapdata) <= 1) stop("nucdiff needs at least two haplotypes to compare")

	# Creating a list that will contain the positions that differ
	output <- list()

	for(i in 1:(length(hapdata) - 1)) {
		if(!all(hapdata[[i]]$pos == hapdata[[i + 1]]$pos)) stop("All positions must be shared between all haplotypes")

		output[[i]] <- which(hapdata[[i]]$nuc != hapdata[[i + 1]]$nuc & hapdata[[i]]$nuc != "-" & hapdata[[i + 1]]$nuc != "-")
	}

	output
}

difflist <- nucdiff(plotting_data)

# A function that sets up viewports for plotting haplotypes and draws them
grid.haplotypes <- function(plotting_data, difflist, fontsize = 8) {

	# Determining the number of columns as the maximum position to plot
	ncolumns <- max(do.call("rbind", plotting_data)$pos)

	# Dividing the viewport into rows
	# Haplotype rows are interleaved with rows used to show the differences
	#  between consecutive haplotypes, hence the number of rows being 2 * (number of haplotypes) - 1
	grid::pushViewport(grid::viewport(layout = grid::grid.layout(nrow = length(plotting_data) * 2 - 1,
								     ncol = ncolumns)))

	# Plotting each haplotype in its row
	for(i in 1:length(plotting_data)) {
		hapdata <- plotting_data[[i]]

		for(j in 1:nrow(hapdata)) {
			grid.text(hapdata[j, "nuc"],
				  gp = grid::gpar(fontsize = fontsize,
						  fontfamily = "mono",
						  fontface = ifelse(hapdata[j, "log10p"] != 0, "bold", "plain"),
						  col = hapdata[j, "color"]),
				  vp = grid::viewport(layout.pos.row = i * 2 - 1,
						      layout.pos.col = hapdata[j, "pos"]))
		}
	}

	# Also plotting the positions where the nucleotides differ
	for(i in 1:length(difflist)) {
		if(!length(difflist[[i]])) next

		for(j in difflist[[i]]) {
			grid::grid.lines(x = 0.5, y = c(0, 1), vp = grid::viewport(layout.pos.row = i * 2, layout.pos.col = j))
		}
	}

	# Moving back into the top viewport
	grid::upViewport()

}

grid.phenotable <- function(haplotype_data) {
	# Computing a table of the phenotypes observed per character
	phenotable <- table(haplotype_data$phenotype_character, haplotype_data$id)

	# Creating a viewport with the required cells to print the data
	grid::pushViewport(grid::viewport(layout = grid::grid.layout(nrow = nrow(phenotable) + 1,
								     ncol = ncol(phenotable) + 1)))

	# Printing the IDs of the phenotypes
	for(i in 1:ncol(phenotable)) grid.text(colnames(phenotable)[i], vp = grid::viewport(layout.pos.row = 1, layout.pos.col = i + 1))

	# Printing the phenotype values themselves
	for(i in 1:nrow(phenotable)) grid.text(rownames(phenotable)[i], vp = grid::viewport(layout.pos.row = i + 1, layout.pos.col = 1))

	# Printing the counts
	for(i in 1:nrow(phenotable)) {
		for(j in 1:ncol(phenotable)) {
			grid.text(as.character(phenotable[i, j]), vp = grid::viewport(layout.pos.row = i + 1, layout.pos.col = j + 1))
		}
	}

	# Returning to the top viewport
	upViewport()

	return(invisible(NULL))
}

# Outputting to a png file
png(paste0("figures/", locus, "_kmers.png"), width = 8, height = 2, units = "in", res = 100)

# Initializing the device
grid.newpage()

# Drawing a box around the plotting region
grid.rect()

# Dividing the viewport into three viewports:
# - The first viewport will be used for showing the haplotype sequences
# - The second viewport will be used for the colour scale of the p-values
# - The third viewport will be used for showing a table of the phenotypes observed per haplotype
pushViewport(viewport(layout = grid.layout(nrow = 3, heights = unit(c(0.3, 0.1, 0.6), "npc"))))


# Moving into the viewport associated with the sequences and drawing them
pushViewport(viewport(layout.pos.row = 1))
grid.haplotypes(plotting_data, difflist, fontsize = 8)
upViewport()

# Moving into the viewport for the color scale
pushViewport(viewport(layout.pos.row = 2))
grid.text("Color scale")
upViewport()

# Moving into the viewport for the table
pushViewport(viewport(layout.pos.row = 3))
grid.phenotable(haplotype_data)
upViewport()

dev.off()

