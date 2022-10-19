# Generating the supplemental data on the SRA accessions metadata

# Reading the metadata from file
# DEPENDENCY: utilities/correct_sra_metadata.csv
sra <- read.table("utilities/correct_sra_metadata.csv", sep = ",", header = TRUE, stringsAsFactors = FALSE)

# Keeping only the necessary columns for the output file
sra <- sra[, c("Run", "spots", "bases", "avgLength", "BioProject", "Model", "BioSample", "SampleName")]

# Adding a column with the ID used as part of this work
sra$study_ID <- NA

# Loading the correspondence between SRA IDs and sample IDs
# DEPENDENCY: utilities/srr_id_correspondence.txt
bayer_ids <- read.table("utilities/srr_id_correspondence.txt", header = FALSE, stringsAsFactors = FALSE)

# Filling the study_ID column
for(i in 1:nrow(sra)) {
	matching_row <- which(grepl(sra[i, "Run"], bayer_ids[[2]]))
	stopifnot(length(matching_row) == 1)
	sra[i, "study_ID"] <- bayer_ids[matching_row, 1]
}

# Checking that all SRA IDs matched a sample and that all 389 samples are represented
stopifnot(all(!is.na(sra$study_ID)))
stopifnot(length(unique(sra$study_ID)) == 389)

# Writing the resulting table to file in CSV format
write.table(sra, file = "additional_files/supplemental_file_1.csv", sep = ",",
	    col.names = TRUE, row.names = FALSE, quote = FALSE)

