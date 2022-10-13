# Loading the GAPIT3 package
library(GAPIT3)

# Loading phenotype data
phenotypes <- read.table("cnv_analysis/cnv_luster_data.csv",
			 header = TRUE, sep = ",",
			 stringsAsFactors = FALSE)[, c("id", "dullshiny")]

colnames(phenotypes) <- c("Taxa", "dullshiny")

# Loading the genotype data
genotypes <- read.table("~/sv_gwas/usda_lines/gwas/platypus/platypus_formatted.hmp.txt",
			comment.char = "",
			header = FALSE)

# Running the GAPIT analysis using 9 PCs and VanRaden kinship
dir.create("cnv_analysis/cnv_luster_gwas", recursive = TRUE)
setwd("cnv_analysis/cnv_luster_gwas")

gapit_model <- GAPIT(Y = phenotypes,
		     G = genotypes,
		     kinship.algorithm = "VanRaden",
		     PCA.total = 9,
		     Model.selection = FALSE, 
		     model = "MLM",
		     file.output = TRUE)

save(gapit_model, file = paste0("gapit_MLM_", phenotype, ".RData"))

