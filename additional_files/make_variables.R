# Creating a list of key-value pairs that will be made available to
# additional file 1 so that they can be read in LaTeX and output
# in the file

suppressMessages(library(GenomicRanges))

# First making a list with the number of k-mers used for LD plots
# DEPENDENCY: LD analysis
ld_files <- dir("gwas_results/kmers", pattern = "clustered_ld.txt$", full.names = TRUE)

ld_length <- sapply(ld_files, function(x) length(readLines(x)))
names(ld_length) <- sub("_clustered_ld.txt", "", basename(names(ld_length)))
names(ld_length) <- paste0(names(ld_length), "_ldkmers")

# Then making a list of the number of samples used for GWAS analyses
# DEPENDENCY: phenotypic_data/phenotypic_data.csv
phenodata <- read.table("phenotypic_data/phenotypic_data.csv", sep = ";", header = TRUE, stringsAsFactors = FALSE)
# DEPENDENCY: utilities/trait_names.txt
trait_names <- readLines("utilities/trait_names.txt")
trait_length <- sapply(trait_names, function(x, dataset) sum(!is.na(dataset[[x]])), dataset = phenodata)
names(trait_length) <- paste0(names(trait_length), "_nsamples")

# Now making a list of the genes associated with each trait/locus combination
# DEPENDENCY: utilities/all_signals.rds
all_signals <- readRDS("utilities/all_signals.rds")
all_signals <- all_signals[!is.na(all_signals$gene_name_v4)]
genes <- all_signals$gene_name_v4
genes <- sub("^g", "G", genes)
names(genes) <- paste0(names(all_signals), "_gene")

# Grouping all the key-value pairs together
key_values <- c(ld_length, trait_length, genes)

# Reformatting them using the " = " separator
key_values <- paste0(names(key_values), " = ", key_values)

writeLines(key_values, con = "additional_files/variables.txt")
