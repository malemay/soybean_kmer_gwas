# Setting the working directory
setwd("gwas_results/kmers/")

# Looping over all traits
# DEPENDENCY: utilities/trait_names.txt
for(phenotype in readLines("../../utilities/trait_names.txt")) {

	# Formatting the phenotypic data for input to the k-mers GWAS program
	# DEPENDENCY: phenotypic_data/phenotypic_data.csv
	phenotypic_data <- read.table("../../phenotypic_data/phenotypic_data.csv", sep = ";", header = TRUE,
				      stringsAsFactors = FALSE)

	# Keeping only the ID column and the phenotype of interest
	phenotypic_data <- phenotypic_data[, c("bayer_id", phenotype)]
	colnames(phenotypic_data) <- c("accession_id", "phenotype_value")
	phenotypic_data <- phenotypic_data[complete.cases(phenotypic_data), ]

	# Writing the table to file
	write.table(phenotypic_data, file = paste0(phenotype, ".pheno"),
		    sep = "\t", col.names = TRUE, row.names = FALSE,
		    quote = FALSE)

	# A function that launches a kmers_gwas analysis
	kmers_gwas <- function(pheno_file, kmers_table, kmer_length, n_kmers, maf,
			       threads, output_dir, kmers_gwas_path, gemma_path) {

		command <- paste0("python2.7 ", kmers_gwas_path,
				  " --pheno ", pheno_file,
				  " --kmers_table ", kmers_table,
				  " -l ", kmer_length,
				  " -k ", n_kmers,
				  " --maf ", maf,
				  " -p ", threads,
				  " --outdir ", output_dir,
				  " --gemma_path ", gemma_path)

		system(command)
	}

	# DEPENDENCY: kmers_table/kmers_table.names
	# DEPENDENCY: kmers_table/kmers_table.table
	# DEPENDENCY: kmers_table/kmers_table.kinship
	kmers_gwas(pheno_file = paste0(phenotype, ".pheno"),
		   kmers_table = "../../kmers_table/kmers_table",
		   kmer_length = 31,
		   n_kmers = 1000001,
		   maf = 0.02,
		   threads = 16,
		   output_dir = phenotype,
		   kmers_gwas_path = "kmers_gwas.py",
		   gemma_path = "gemma_0_96")
}

file.create("KMER_GWAS")

