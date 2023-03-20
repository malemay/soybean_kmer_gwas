# Using a custom R function to re-select the SVs clustered by SVmerge by preferentially
# choosing realigned variants

# A function that takes a vcf file merged by SVmerge and uses the information
# in it to identify variant clusters and re-select variants from those. The
# criteria to select variants will be as follows:

# - If there is at least one REALIGNED variant in the cluster, select randomly from those
# - Otherwise, select randomly from all variants in the cluster

# svmerge_file is a vcf file merged by SVmerge from which the names of the variants in the cluster can be extracted
# input_vcfs is a character vector of vcf file names from which the IDs are to be queried
# output_vcf is the name of the file to which the output will be written
merge_aligned_variants <- function(svmerge_file, input_vcfs, output_vcf) {

	# Reading the svmerge file
	svmerge <- scan(svmerge_file, what = character(), sep = "\n", quiet = TRUE, 
			comment.char = "#")

	# Extracting the ClusterIDs INFO field for each record
	svmerge <- sub(".*ClusterIDs=([^;]*);.*", "\\1", svmerge)

	# Reading all the vcf files in memory and putting them in a single data.frame
	vcf_data <- list()

	for(i in input_vcfs) {
		vcf_data[[i]] <- read.table(i, header = FALSE, sep = "\t", comment.char = "#", stringsAsFactors = FALSE)
	}

	vcf_data <- do.call("rbind", vcf_data)
	colnames(vcf_data) <- c("#CHROM", "POS", "ID", "REF", "ALT", 
				"QUAL", "FILTER", "INFO", "FORMAT", "GT")

	# Opening the output file for writing, because lines will be written one by one
	output <- file(output_vcf,, open = "w")
	on.exit(close(output), add = FALSE)

	# We read the header from the first of input_vcfs so we can write it to output
	header <- scan(input_vcfs[1], what = character(), sep = "\n", quiet = TRUE)
	header <- grep("^##", header, value = TRUE)

	cat(header, sep = "\n", file = output)
	cat(names(vcf_data)[-c(9, 10)], "\n", sep = "\t", file = output)

	# Looping over the lines in the svmerge vector
	for(i in svmerge) {
		ids <- strsplit(i, ":")[[1]]

		# Extracting the lines with a matching ID
		i_lines <- vcf_data[vcf_data$ID %in% ids, ]

		# Sanity check
		stopifnot(length(ids) == nrow(i_lines))

		# Checking if any of the lines have been realigned
		i_realigned <- grep("REALIGNED", i_lines$INFO)

		# We sample from those if any has been realigned
		if(length(i_realigned)) {
			i_lines <- i_lines[i_realigned, ]
		}

		# Sampling one line from those available
		i_line <- i_lines[sample(nrow(i_lines), 1), ]	

		# Writing it to the output file
		cat(as.character(i_line[, -c(9, 10)]), "\n", sep = "\t", file = output)

	}

	return(invisible(svmerge))
}

# Generating a character vector of input vcfs
input_vcfs <- scan("files.txt", what = character(), sep = "\n", quiet = TRUE)

# Lauching the command
merge_aligned_variants("svmerged_preliminary.clustered.vcf",
		       input_vcfs,
		       "svmerged.clustered.vcf")

