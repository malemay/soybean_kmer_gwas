# Loading the GAPIT3 package
library(GAPIT3)

# Setting the working directory
setwd("gwas_results/paragraph/")

# Using maturity as a phenotype
phenotype <- commandArgs(trailingOnly = TRUE)[1]

# Loading phenotype data (from the same file as TASSEL)
# DEPENDENCY: phenotypic_data/phenotypic_data.csv
phenotypes <- read.table("../../phenotypic_data/phenotypic_data.csv",
			 header = TRUE, sep = ";",
			 stringsAsFactors = FALSE)[, c("bayer_id", phenotype)]

phenotypes <- phenotypes[complete.cases(phenotypes), ]

colnames(phenotypes) <- c("Taxa", phenotype)

# Loading the genotype data
# DEPENDENCY: sv_genotyping/paragraph/paragraph_formatted.hmp.txt
genotypes <- read.table("../../sv_genotyping/paragraph/paragraph_formatted.hmp.txt",
			comment.char = "",
			header = FALSE)

# Running the GAPIT analysis using 9 PCs and VanRaden kinship
dir.create(phenotype, recursive = TRUE)
setwd(phenotype)

gapit_model <- GAPIT(Y = phenotypes,
		     G = genotypes,
		     kinship.algorithm = "VanRaden",
		     PCA.total = 9,
		     Model.selection = FALSE, 
		     model = "MLM",
		     file.output = TRUE)

save(gapit_model, file = paste0("gapit_MLM_", phenotype, ".RData"))

# Copying the output CSV file to the curent directory
# OUTPUT: gwas_results/paragraph/%_gwas.csv
file.copy(paste0("GAPIT.MLM.", phenotype, ".GWAS.Results.csv"), paste0("../", phenotype, "_gwas.csv"))

