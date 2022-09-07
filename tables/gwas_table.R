# Creating a .csv file of phenotypic observations suitable for inclusion as a supplementary table

# Getting the name of the trait from the command line
trait <- commandArgs(trailingOnly = TRUE)[1]

# Loading the lookup table of correspondence between GRIN trait names and the ones we use
# DEPENDENCY: phenotypic_data/trait_names.rds
trait_names <- readRDS("phenotypic_data/trait_names.rds")
trait_prefix <- trait_names[trait]

# DEPENDENCY: phenotypic_data/phenotypic_data.csv
phenotypes <- read.table("phenotypic_data/phenotypic_data.csv", sep = ";", header = TRUE)

# DEPENDENCY: phenotypic_data/lookup_tables.rds
lookup_tables <- readRDS("phenotypic_data/lookup_tables.rds")
lookup_tables <- lookup_tables[grepl(paste0("^", trait_prefix), names(lookup_tables))]

# A function that takes a data.frame of phenotypic observations and a set of lookup tables
# and outputs a data.frame indicating the frequency of each phenotype and the code used in
# various GWAS analyses
format_table <- function(phenodata, phenotype, lookup_tables) {
	x <- as.data.frame(table(phenodata[[phenotype]], useNA = "ifany"), stringsAsFactors = FALSE)
	colnames(x) <- c("Value", "Frequency")

	for(i in 1:length(lookup_tables)) {
		x[[paste0("GWAS", letters[i])]] <- lookup_tables[[i]][x$Value]
	}

	return(x)
}

output_table <- format_table(phenotypes, trait, lookup_tables)

write.table(output_table, sep = ";", col.names = TRUE, quote = FALSE, row.names = FALSE,
	    file = paste0("tables/", trait, "_gwas_table.csv"))

