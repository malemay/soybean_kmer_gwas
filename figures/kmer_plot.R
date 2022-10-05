# This script generates a k-mer plot from the consensus sequences of a group of samples at a given locus

# Loading the required libraries
suppressMessages(library(grid))
suppressMessages(library(gwastools))
suppressMessages(library(GenomicRanges))

# Setting some analysis parameters
max_kmers <- -1 # the maximum number of k-mers to map onto the haplotype sequences; -1 means all k-mers
kmer_length <- 31 # the length of the k-mers used in the analysis
min_sequences <- 80 # min number of samples that must have a consensus sequence from assembly; otherwise obtained from bwa
min_frequency <- 5 # min number of times that a haplotype must occur to be used for plotting
xextend <- 0.1 # The expansion factor on either side of the gene for the transcript plot x-scale

# Reading the locus from the command line
locus <- commandArgs(trailingOnly = TRUE)[1]

# Getting the trait that corresponds to that locus
# DEPENDENCY: utilities/kmer_plot_ranges.txt
locus_data <- read.table("utilities/kmer_plot_ranges.txt")
trait <- locus_data[locus_data[[1]] == locus, 2]

# Also extracting the plotting range
plotting_range <- locus_data[locus_data[[1]] == locus, 3]
chrom <- sub(":.*", "", plotting_range)
grange <- as.numeric(strsplit(sub(".*:", "", plotting_range), "-")[[1]])

# Getting the lookup table that links the trait name used in the USDA database to the one used for our analyses
# DEPENDENCY: phenotypic_data/trait_names.rds
trait_names <- readRDS("phenotypic_data/trait_names.rds")
usda_trait <- names(trait_names)[sapply(trait_names, function(x) grepl(paste0("^", x), trait))]

# DEPENDENCY: phenotypic_data/pheno_names_lookup.rds
pheno_names_lookup <- readRDS("phenotypic_data/pheno_names_lookup.rds")[[usda_trait]]

# Loading the reference gene and transcript data for plotting the transcript
genes <- readRDS("refgenome/gmax_v4_genes.rds")
transcripts <- readRDS("refgenome/gmax_v4_transcripts.rds")
cds <- readRDS("refgenome/gmax_v4_cds.rds")
exons <- readRDS("refgenome/gmax_v4_exons.rds")

# Loading the data corresponding to that signal so we can extract the causal gene
target_signal <- readRDS("utilities/all_signals.rds")[locus]
stopifnot(length(target_signal) == 1)

gene_name <- target_signal$gene_name_v4

# Handling the special case where the gene_name value contains more than one gene or no gene at all
if(!is.na(gene_name) && grepl(";", gene_name)) {
	gene_names <- strsplit(gene_name, ";")[[1]]
	causal_gene <- genes[gene_names]
	causal_gene <- reduce(causal_gene, min.gapwidth = 10^6, ignore.strand = TRUE)
} else if (!is.na(gene_name)) {
	causal_gene <- genes[gene_name]
} else {
	gene <- NULL
}

if(is.null(causal_gene)) stop("No gene found where there should be one")

# Extending the window for the causal gene
window_size <- GenomicRanges::width(causal_gene)
start(causal_gene) <- start(causal_gene) - xextend * window_size
end(causal_gene) <- end(causal_gene) + xextend * window_size

# Getting the phenotypic data
# DEPENDENCY: phenotypic_data/phenotypic_data.csv
phenotypes <- read.table("phenotypic_data/phenotypic_data.csv", sep = ";", header = TRUE,)

# Setting the values of the phenotype to that in pheno_names_lookup
phenotypes[[usda_trait]] <- pheno_names_lookup[phenotypes[[usda_trait]]]

# Setting to NA the values that were NA in the analysis being considered
phenotypes[[usda_trait]][is.na(phenotypes[[trait]])] <- NA

# Reading the k-mers associated with that analysis and their p-values
# DEPENDENCY: k-mer p-values
kmer_pvalues <- read_kmer_pvalues(kmer_file = paste0("gwas_results/kmer_data/", trait,
						     "/katcher_results/pass_threshold_5per_sorted.txt"),
				  max_kmers = max_kmers,
				  kmer_length = 31)

# Reading the consensus sequence from the assemblies
# DEPENDENCY: consensus sequence computed from the assembly of significant reads
sequences <- read_fasta(input_fasta = paste0("gwas_results/kmer_consensus/", locus, "_sequences.fa"))

# DEPENDENCY: consensus sequence computed from the alignment using BWA
# Using consensus sequences from bwa in case the assemblies did not yield great results
if(length(sequences) < min_sequences) {
	warning("Only ", length(sequences), "have a consensus sequence from assembly.",
	       	"Using consensus sequences from read alignment by bwa for locus ", locus)
	sequences <- read_fasta(paste0("gwas_results/kmer_consensus/", locus, "_bwa_sequences.fa"))
}

# Get a character vector of the haplotypes found with min_frequency in the dataset
haplotypes <- get_haplotypes(sequences = sequences, min_frequency = min_frequency)

# Creating a data.frame linking haplotypes to phenotypes
haplotype_data <- link_phenotypes(sequences = sequences,
				  haplotypes = haplotypes,
				  phenotypes = phenotypes,
				  id_column = "bayer_id",
				  phenotype_column = usda_trait)

# Writing this data.frame to file for retrieval later
saveRDS(haplotype_data, file = paste0("gwas_results/kmers/", locus, "_phenodata.rds"))

# Get the positions in the haplotypes that overlap any of the significant k-mers
kmer_overlaps <- match_kmers(sequences = haplotypes,
			     kmers = kmer_pvalues,
			     kmer_length = kmer_length,
			     data_columns = "log10p")

# Generating a data.frame with separate nucleotides and their associated p-value for each haplotype
plotting_data <- format_haplotypes(haplotypes = haplotypes,
				   overlaps = kmer_overlaps)

# Performing multiple alignment of the sequences to find the gaps
alignment <- mafft_align(fasta_path = paste0("gwas_results/kmer_consensus/", locus, "_haplotypes.fa"),
			 sequences = haplotypes,
			 mafft_path = "mafft",
			 mafft_options = "--auto --quiet")

# Adjusting the positions for plotting depending on the gaps in the alignment
plotting_data <- adjust_gaps(hapdata = plotting_data,
			     alignment = alignment)

# Finding the positions where the aligned nucleotides differ between haplotypes
difflist <- nucdiff(plotting_data)

# Outputting to a png file
png(paste0("figures/", locus, "_kmers.png"), width = 8, height = 8, units = "in", res = 100)

# Initializing the device
grid.newpage()

# Drawing a box around the plotting region
grid::grid.rect()

# Dividing the viewport into five row viewports:
# - The first viewport is used to plot the transcript for the gene and the highlighted region
# - The second viewport acts as a buffer between the first and third
# - The third viewport is used for showing the haplotype sequences
# - The fourth viewport acts as a buffer between the third and fifth
# - The fifth viewport is used for showing a table of the phenotypes observed per haplotype
grid::pushViewport(grid::viewport(layout = grid::grid.layout(nrow = 5, heights = grid::unit(c(0.1, 0.15, 0.3, 0.15, 0.3), "npc"))))

# Plotting the transcript in the first viewport
grid::pushViewport(grid::viewport(width = 0.84, layout.pos.row = 1))
grid.draw(transcriptsGrob(genes = genes,
			  transcripts = transcripts,
			  exons = exons,
			  cds = cds,
			  xscale = causal_gene,
			  highlight = GenomicRanges::GRanges(seqnames = chrom,
							     ranges = IRanges::IRanges(start = grange[1], end = grange[2])),
			  draw_arrows = TRUE,
			  first_tx_only = TRUE))
grid::upViewport()

# Moving into the viewport associated with the sequences and drawing them
grid::pushViewport(grid::viewport(layout.pos.row = 3))
grid.haplotypes(hapdata = plotting_data, difflist = difflist, position = plotting_range)
grid::upViewport()

# Moving into the viewport for the table
grid::pushViewport(grid::viewport(layout.pos.row = 5))
grid.phenotable(phenodata = haplotype_data)
grid::upViewport()

dev.off()

