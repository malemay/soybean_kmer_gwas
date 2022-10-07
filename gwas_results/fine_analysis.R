# This script performs fine analyses for the results to mention
# in the manuscript on a trait-specific basis.

# Loading the required libraries
library(GenomicRanges)
library(VariantAnnotation)
library(Biostrings)

# First we will create a few functions that will be used throughout

# A function that reads the gwas results for all programs and
# sorts them be decreasing p-value
read_gwas <- function(trait, programs = c("platypus", "vg", "paragraph", "kmers")) {
	output <- lapply(programs, function(x) readRDS(paste0("gwas_results/", x, "/", trait, "_gwas_subset.rds")))
	output <- lapply(output, function(x) x[order(x$log10p, decreasing = TRUE)])
	output <- lapply(output, function(x) {x$rank <- rank(-x$log10p, ties.method = "min"); x})
	names(output) <- programs
	output
}

# A function that extracts the markers overlapping the gene associated with a locus
gene_markers <- function(markers, signals, genes, locus) {
	stopifnot(locus %in% names(signals))
	signal <- signals[locus]

	stopifnot(signal$gene_name_v4 %in% names(genes))
	gene <- genes[signal$gene_name_v4]

	lapply(markers, function(x) subsetByOverlaps(x, gene, ignore.strand = TRUE))
}

# A function that extracts a marker from the original vcf file
vcf_extract <- function(vcf_path, variant) {
	vcf_file <- VcfFile(vcf_path)
	vcf_data <- readVcf(vcf_file, param = variant)
	vcf_data[names(variant)]
}

# Loading the dataset of all signals
all_signals <- readRDS("utilities/all_signals.rds")

# Loading the set of genes
genes <- readRDS("refgenome/gmax_v4_genes.rds")

# The paths to the filtered variant files
platypus_file <- "filtered_variants/platypus/filtered_variants.vcf.gz"
paragraph_file <- "filtered_variants/paragraph/filtered_variants.vcf.gz"
vg_file <- "filtered_variants/vg/filtered_variants.vcf.gz"




### ---------- FLOWER COLOR
# Loading the gwas results for this trait
flower_color <- read_gwas("flower_color")
gene_markers(flower_color, all_signals, genes, "flower_color_W1")

# Having a look at the top variants for vg and Paragraph
vcfFixed(vg_w1 <- vcf_extract(vg_file, flower_color$vg[1]))
vcfInfo(vg_w1)
vcfFixed(paragraph_w1 <- vcf_extract(paragraph_file, flower_color$paragraph[1]))
vcfInfo(paragraph_w1)

### ---------- PUBESCENCE COLOR ALL
pubescence_color_all <- read_gwas("pubescence_color_all")
pc_all_top_markers <- gene_markers(pubescence_color_all, all_signals, genes, "pubescence_color_all_T")


# The indel reported by Zabala and Vodkin (2003) comes 188th in Platypus but is
# the most significantly associated within the gene.

# The most associated k-mers in the gene are in 24th rank overall, but 4th and 5th
# rank within the gene. What I would like to find out is what variation underlies
# the top three k-mers in the gene (1st, 2nd and 4th overall)

# Let us look at the sequence of those k-mers
pc_all_top_markers$kmers$kmer_canon[1:3]
# [1] "AATGCGATGACAAAGATTGTAATTTTAAAGC" "AGAAAATAAAAATAAAAAATGAAAACAGAAG" "AAAAAATAAAAATAAAAAATGAAAACAGAAG"
# It could also match their reverse complement
as.character(reverseComplement(DNAStringSet(pc_all_top_markers$kmers$kmer_canon[1:3])))
# [1] "GCTTTAAAATTACAATCTTTGTCATCGCATT" "CTTCTGTTTTCATTTTTTATTTTTATTTTCT" "CTTCTGTTTTCATTTTTTATTTTTATTTTTT"

# Let us see what the reads that they match to look like
system("samtools view gwas_results/kmer_data/pubescence_color_all/katcher_results/USB-018/USB-018_pvalues_sorted.bam | grep AATGCGATGACAAAGATTGTAATTTTAAAGC | less -S")
system("samtools view gwas_results/kmer_data/pubescence_color_all/katcher_results/USB-018/USB-018_pvalues_sorted.bam | grep GCTTTAAAATTACAATCTTTGTCATCGCATT | less -S")
# this first k-mer does not seem to show anything special


system("samtools view gwas_results/kmer_data/pubescence_color_all/katcher_results/USB-029/USB-029_pvalues_sorted.bam | grep AGAAAATAAAAATAAAAAATGAAAACAGAAG | less -S")
system("samtools view gwas_results/kmer_data/pubescence_color_all/katcher_results/USB-029/USB-029_pvalues_sorted.bam | grep CTTCTGTTTTCATTTTTTATTTTTATTTTCT | less -S")
# Nor does this one


system("samtools view gwas_results/kmer_data/pubescence_color_all/katcher_results/USB-512/USB-512_pvalues_sorted.bam | grep AAAAAATAAAAATAAAAAATGAAAACAGAAG | less -S")
system("samtools view gwas_results/kmer_data/pubescence_color_all/katcher_results/USB-512/USB-512_pvalues_sorted.bam | grep CTTCTGTTTTCATTTTTTATTTTTATTTTTT | less -S")
# There seems to be an insertion associated with this last k-mer
