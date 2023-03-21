# Setting the working directory
setwd("gwas_results/platypus/")

# DEPENDENCY: Platypus GWAS analyses

# Looping over the output of all GAPIT analyses, checking if the PCA and kinship results are the same for all, and outputting to file
rdata_files <- dir(".", pattern = ".*\\.RData$", recursive = TRUE)
stopifnot(length(rdata_files) == 22)
kinship <- list()
pca <- list()

for(i in 1:length(rdata_files)) {
	message("Reading ", rdata_files[i])
	load(rdata_files[i])
	kinship[[i]] <- gapit_model$KI
	pca[[i]] <- gapit_model$PCA
}

for(i in 2:length(kinship)) stopifnot(identical(kinship[[1]], kinship[[i]]))
for(i in 2:length(pca)) stopifnot(identical(pca[[1]], pca[[i]]))

# Saving the kinship and pca to an .rds file
setwd("../../filtered_variants/")
saveRDS(kinship[[1]], file = "platypus_gapit_kinship.rds")
saveRDS(pca[[1]], file = "platypus_gapit_pca.rds")

