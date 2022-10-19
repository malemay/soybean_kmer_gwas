# Preparing a supplemental table with the phenotypic data used bu our GWAS analyses

# Reading the main phenotypic data table
phenodata <- read.table("phenotypic_data/phenotypic_data.csv", sep = ";",
			header = TRUE, stringsAsFactors = FALSE)

# Changing the names of some columns
colnames(phenodata)[colnames(phenodata) == "bayer_id"] <- "study_ID"

# Removing some other columns
phenodata <- phenodata[, !colnames(phenodata) %in% c("TAXONOMY", "N_CDW")]

# Writing the resulting table to file
write.table(phenodata, file = "additional_files/supplemental_file_3.csv", sep = ";",
	    col.names = TRUE, row.names = FALSE, quote = FALSE)
