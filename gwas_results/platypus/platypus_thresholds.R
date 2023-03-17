# Loading the GAPIT3 package
library(GAPIT3)
library(parallel)

# Setting the working directory
setwd("gwas_results/platypus/")

# Getting the phenotype from the command line
phenotype <- commandArgs(trailingOnly = TRUE)[1]

# Loading phenotype data (from the same file as TASSEL)
# DEPENDENCY: phenotypic_data/phenotypic_data.csv
phenotypes <- read.table("../../phenotypic_data/phenotypic_data.csv",
			 header = TRUE, sep = ";",
			 stringsAsFactors = FALSE)[, c("bayer_id", phenotype)]

phenotypes <- phenotypes[complete.cases(phenotypes), ]

colnames(phenotypes) <- c("Taxa", phenotype)

# Loading the genotype data
# DEPENDENCY: variant_calling/platypus/platypus_formatted.hmp.txt
genotypes <- read.table("../../variant_calling/platypus/platypus_formatted.hmp.txt",
			comment.char = "",
			header = FALSE)

# Creating a directory to store the permutations
dir.create(paste0(phenotype, "_permutations"), recursive = TRUE)
setwd(paste0(phenotype, "_permutations"))

# Creating the phenotype permutations outside of the mclapply function call
# because it does not deal well with random number generation in parallel
shuffled_phenotypes <- list()

for(i in 1:100) {
	shuffled_phenotypes[[i]] <- phenotypes
	shuffled_phenotypes[[i]][[phenotype]] <- sample(shuffled_phenotypes[[i]][[phenotype]])
}

# Using the mclapply function from the parallel package to run the permutations on 20 different cores
mclapply(1:100,
	 FUN = function(x, phenotype_data, phenotype, genotypes) {
		 # Creating the output directory and moving to it
		 dir.create(paste0("perm_", x))
		 setwd(paste0("perm_", x))
		 # Launching the GAPIT analysis
		 gapit_model <- GAPIT(Y = phenotype_data[[x]],
				      G = genotypes,
				      kinship.algorithm = "VanRaden",
				      PCA.total = 9,
				      Model.selection = FALSE, 
				      model = "MLM",
				      file.output = TRUE)

		 save(gapit_model, file = paste0("gapit_MLM_", phenotype, "_perm", x, ".RData"))
		 setwd("..")
	 },
	 phenotype_data = shuffled_phenotypes,
	 phenotype = phenotype,
	 genotypes = genotypes,
	 mc.cores = 5)

# Getting the threshold for every trait by picking the 5th lowest p-value
# from the top p-values of all 100 permutations

# Going back to the main directory
setwd("..")

results_files <- dir(paste0(phenotype, "_permutations"),
		     pattern = "GAPIT\\.MLM\\..*\\.GWAS\\.Results\\.csv",
		     recursive = TRUE,
		     full.names = TRUE)

stopifnot(length(results_files) == 100)

# Reading all the results
top_pvalues <- mclapply(results_files,
			FUN = function(x) {
				pvalues <- read.table(x, header = TRUE, sep = ",",
						      colClasses = c(rep("NULL", 3), "numeric", rep("NULL", 6)))
				pvalues <- pvalues[[1]]
				min(pvalues)
			},
			mc.cores = 10)

top_pvalues <- unlist(top_pvalues)

# Making sure that randomization worked properly
stopifnot(length(unique(top_pvalues)) > 90)

top_pvalues <- sort(top_pvalues)
threshold <- top_pvalues[5]

# Writing the threshold to a file in the appropriate directory
# OUTPUT: gwas_results/platypus/%_threshold_5per.txt
cat(threshold, "\n", file = paste0(phenotype, "_threshold_5per.txt"))

