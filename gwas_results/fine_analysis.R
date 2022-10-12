# This script performs fine analyses for the results to mention
# in the manuscript on a trait-specific basis.

# Loading the required libraries
library(GenomicRanges)
library(VariantAnnotation)
library(Biostrings)
library(GenomicFeatures)
library(Rsamtools)

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
platypus_file <- "filtered_variants/platypus_full.vcf.gz"
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

# I haven't been able to link those k-mers to any variant

# Let us have a look at the two significant markers found using Paragraph
vcfFixed(paragraph_t <- vcf_extract(paragraph_file, pc_all_top_markers$paragraph[1:2]))
vcfInfo(paragraph_t)
width(vcfFixed(paragraph_t)$REF)
nchar(vcfFixed(paragraph_t)$ALT)




### ---------- PUBESCENCE COLOR NO GRAY
pubescence_color_nogray <- read_gwas("pubescence_color_nogray")
pc_nogray_top_markers <- gene_markers(pubescence_color_nogray, all_signals, genes, "pubescence_color_nogray_Td")
vcfFixed(platypus_td <- vcf_extract(platypus_file, pc_nogray_top_markers$platypus[1]))


### ---------- SEED COAT COLOR GREEN YELLOW
seed_coat_color_greenyellow <- read_gwas("seed_coat_color_greenyellow")
scc_gy <- gene_markers(seed_coat_color_greenyellow, all_signals, genes, "seed_coat_color_greenyellow_G")
vcfFixed(platypus_g <- vcf_extract(platypus_file, scc_gy$platypus[1]))

# We need to check that this position indeed corresponds to the one reported by Wang et al. (2018)
# their position corresponded to Chr01:53229579 in version 2 of the genome
system("samtools faidx ~/refgenome/Gmax_nuclv2_mit_chlp.fasta Chr01:53229560-53229600")


### ---------- STEM TERMINATION ALL
stem_termination_all <- read_gwas("stem_termination_all")
stall_gene <- gene_markers(stem_termination_all, all_signals, genes, "stem_termination_all_Dt1")
vcfFixed(platypus_dt1 <- vcf_extract(platypus_file, stall_gene$platypus[1]))

### ---------- STEM TERMINATION semi-determinate/indeterminate
stem_termination_sn <- read_gwas("stem_termination_sn")

# For this trait we found signals on Gm19 (Dt1) plus three unknown signals on Gm11, Gm16 and Gm18
# Given that Dt2 is found on Gm18, let us see if that might be the signal we found
all_signals[all_signals$locus == "Dt2"]
genes[all_signals[all_signals$locus == "Dt2" & all_signals$trait == "stem_termination_sn"]$gene_name_v4]
stem_termination_sn$kmers[seqnames(stem_termination_sn$kmers) == "Gm18"]
stem_termination_sn_signals <- readRDS("gwas_results/kmers/stem_termination_sn_signal.rds")
stem_termination_sn_signals[seqnames(stem_termination_sn_signals) == "Gm18"]
# This does not seem to correspond to the location of neither the Dt2 signal nor the Dt2 putative gene
# We will create a new signal for it and plot the results

# Now let us have a look at the signal on chromosome 16
stem_termination_sn_signals[seqnames(stem_termination_sn_signals) == "Gm16"]
# It spans roughly 100 kb ; we will add it to the list of signals as we did for the one on Gm18

# the signal on chromsome Gm11 will also be added to the list of signals
stem_termination_sn_signals[seqnames(stem_termination_sn_signals) == "Gm11"]


### ---------- HILUM COLOR
hilum_color_all <- read_gwas("hilum_color_all")
hilum_color_gene_t <- gene_markers(hilum_color_all, all_signals, genes, "hilum_color_all_T")

# I would like to see if the k-mers near the boundary of the tandem duplication/inversion at the I locus
# may be associated with the causal variant
subsetByOverlaps(hilum_color_all$kmers,
		 GRanges(seqnames = "Gm08", ranges = IRanges(start = 8525000, end = 8530000)),
		 ignore.strand = TRUE)[1]$kmer_canon
# [1] "CTGTCTGCCCATAAGCACCTGAATTCAGCCA"
as.character(reverseComplement(DNAString("CTGTCTGCCCATAAGCACCTGAATTCAGCCA")))

# Let us see that these k-mers match in the reads
system("samtools view gwas_results/kmer_data/hilum_color_all/katcher_results/USB-022/USB-022_pvalues_sorted.bam | grep CTGTCTGCCCATAAGCACCTGAATTCAGCCA | less -S")
system("samtools view gwas_results/kmer_data/hilum_color_all/katcher_results/USB-022/USB-022_pvalues_sorted.bam | grep TGGCTGAATTCAGGTGCTTATGGGCAGACAG | less -S")

# There is no obvious variant associated with this k-mer


### ---------- HILUM COLOR (black/brown)
hilum_color_blbr <- read_gwas("hilum_color_blackbrown")
hilum_color_gene_r <- gene_markers(hilum_color_blbr, all_signals, genes, "hilum_color_blackbrown_R")

# Trying to see whether there is a variant associated with the top k-mer
hilum_color_blbr$kmers[1:5]
# It looks like the 1st and 4th k-mer are associated with the peak
hilum_color_blbr$kmers[c(1, 4)]$kmer_canon
# [1] "CCAGCTATGGCAGTGAAAACGATGGTAATGG" "ACCACTGCCACCATTACCATCGTTTTCACTA"
# Extracting their reverse complement
as.character(reverseComplement(DNAStringSet(hilum_color_blbr$kmers[c(1, 4)]$kmer_canon)))
# [1] "CCATTACCATCGTTTTCACTGCCATAGCTGG" "TAGTGAAAACGATGGTAATGGTGGCAGTGGT"
system("samtools view gwas_results/kmer_data/hilum_color_blackbrown/katcher_results/USB-038/USB-038_pvalues_sorted.bam | grep CCAGCTATGGCAGTGAAAACGATGGTAATGG | less -S")


system("samtools view gwas_results/kmer_data/hilum_color_blackbrown/katcher_results/USB-029/USB-029_pvalues_sorted.bam | grep ACCACTGCCACCATTACCATCGTTTTCACTA | less -S")
system("samtools view gwas_results/kmer_data/hilum_color_blackbrown/katcher_results/USB-029/USB-029_pvalues_sorted.bam | grep TAGTGAAAACGATGGTAATGGTGGCAGTGGT | less -S")

# The underlying variation does not appear clear from the results. Maybe a k-mer plot will given clearer results
# The k-mer plot showed that there were three haplotypes at that site, with one being clearly associated with the trait
# The presence of three different haplotypes was due to two SNPs occurring right next to one another

# Let us try and see if we can find any of the variations identified by Gillman et al. (2011) in our data
hilum_color_gene_r$platypus
# It looks like the 5th most significant within the gene body might correspond to one of those variants
hilum_color_gene_r$platypus[5]
vcfFixed(platypus_r <- vcf_extract(platypus_file, hilum_color_gene_r$platypus[5]))


### ---------- POD COLOR (all phenotypes)
# There are no results sufficiently clear regarding pod color to analyse in detail


### ---------- POD COLOR (black/brown)
pod_color_blbr <- read_gwas("pod_color_blbr")
pod_color_blbr_signals <- readRDS("gwas_results/kmers/pod_color_blbr_signal.rds")
pod_color_blbr_signals[seqnames(pod_color_blbr_signals) == "Gm15"]
# Using the coordinates of this signal to define the signal



# PUBESCENCE FORM (all)
pub_form <- read_gwas("pubescence_form_all")
pub_form_pa1 <- gene_markers(pub_form, all_signals, genes, "pubescence_form_all_Pa1")

# Let us have a look at the top paragraph variant
vcfFixed(paragraph_pa1 <- vcf_extract(paragraph_file, pub_form$paragraph[1]))


# I want to check whether the variants overlapping this gene could have a potential impact on gene structure
platypus_variants_pa1 <- vcf_extract(platypus_file, pub_form_pa1$platypus)[1:2]
annotations <- makeTxDbFromGFF("~/refgenome/Gmax_v4/Wm82.a4.v1/annotation/Gmax_508_Wm82.a4.v1.gene_exons.gff3")
predictCoding(query = platypus_variants_pa1,
	      subject = annotations,
	      seqSource = FaFile("~/refgenome/Gmax_v4/Gmax_508_v4.0_mit_chlp.fasta"))
# Both variants are non-synonymous
# Let us have a look at the phenotypes depending on the genotype at BOTH loci
pa1_genotypes <- data.frame(id = colnames(geno(platypus_variants_pa1)$GT),
			    haplotype = paste0(geno(platypus_variants_pa1)$GT[1, ], ":", geno(platypus_variants_pa1)$GT[2, ]))

# Reading the phenotypic data
phenodata <- read.table("phenotypic_data/phenotypic_data.csv", sep = ";", header = TRUE, stringsAsFactors = FALSE)
pa1_genotypes$pubform <- phenodata[match(pa1_genotypes$id, phenodata$bayer_id), "pubescence_form_all"]

# Now let us look at a contingency table
table(pa1_genotypes$pubform, pa1_genotypes$haplotype)



# PUBESCENCE DENSITY
pubescence_density <- read_gwas("pubescence_density")
pubescence_density$kmers[1]$kmer_canon
# [1] "AAAATTAATGATATTTTTTTGTAAAAATTAT"

# We want to look at the reads associated with the most significant k-mer
system("samtools view gwas_results/kmer_data/pubescence_density/katcher_results/USB-672/USB-672_pvalues_sorted.bam | grep AAAATTAATGATATTTTTTTGTAAAAATTAT | less -S")
# This k-mer is clearly associated with the CNV at the Ps locus

