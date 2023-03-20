# The objective of this code is to read the vcf file output by the default behaviour
# of SVmerge and output two different vcf files from it:
# - One for which Illumina variants are systematically favoured
# - Another for which Oxford Nanopore variants are systematically favoured

# Only the one which systematically favours illumina will be used
# The one that favours Nanopore is legacy from a different project
select_svs <- function(svmerge_file, input_illumina, input_nanopore, output_illumina, output_nanopore) {

	# Reading the Illumina and Nanopore vcf files in memory as data.frames
	vcfs <- list()
	vcfs[["illumina"]] <- read.table(input_illumina, header = FALSE, sep = "\t", comment.char = "#", stringsAsFactors = FALSE,
					 col.names = c("CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO"))
	vcfs[["nanopore"]] <- read.table(input_nanopore, header = FALSE, sep = "\t", comment.char = "#", stringsAsFactors = FALSE,
					 col.names = c("CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO"))
	vcfs <- do.call("rbind", vcfs)

	# Opening the output files for writing, because lines will be written one by one
	output_illumina <- file(output_illumina, open = "w")
	on.exit(close(output_illumina), add = TRUE)

	output_nanopore <- file(output_nanopore, open = "w")
	on.exit(close(output_nanopore), add = TRUE)

	# We read the input from the svmerge file line by line because it is a large file
	svmerge <- file(svmerge_file, open = "rt")
	on.exit(close(svmerge), add = TRUE)

	while(length(current_line <- scan(svmerge, what = character(), sep = "\n", n = 1, quiet = TRUE))) {

		# If it is a header line, we write it to both output files and jump to the next line
		if(grepl("^#", current_line)) {
			cat(current_line, "\n", sep = "", file = output_illumina)
			cat(current_line, "\n", sep = "", file = output_nanopore)
			next
		}

		# Extracting the cluster_ids from that line
		cluster_ids <- sub(".*ClusterIDs=([^;]*);.*", "\\1", current_line)
		cluster_ids_split <- strsplit(cluster_ids, ":")[[1]]

		# The first case (most common one) is the one for which there is only a single variant
		if(length(cluster_ids_split) == 1) {
			# Then we must get the row corresponding to that ID
			output_row <- vcfs[vcfs$ID == cluster_ids_split, ]
			output_row$INFO <- paste0(output_row$INFO, ";ClusterIDs=", cluster_ids)
			# Outputting to both files and going to the next line
			cat(as.character(output_row), "\n", sep = "\t", file = output_illumina)
			cat(as.character(output_row), "\n", sep = "\t", file = output_nanopore)
			next

		# Then we take care of the situation where there is only one type of variant
		} else if(all(grepl("illumina", cluster_ids_split)) || all(grepl("nanopore", cluster_ids_split))) {
			# Then we randomly select a variant which we output to both files
			output_row <- vcfs[vcfs$ID == sample(cluster_ids_split, 1), ]
			output_row$INFO <- paste0(output_row$INFO, ";ClusterIDs=", cluster_ids)
			# Outputting to both files and going to the next line
			cat(as.character(output_row), "\n", sep = "\t", file = output_illumina)
			cat(as.character(output_row), "\n", sep = "\t", file = output_nanopore)
			next

		# Finally we take care of the situation where we have a mixture of both types
		# Then we sample an illumina variant for the illumina file and a Nanopore variant
		# for the nanopore file
		} else {
			illumina_ids <- grep("illumina", cluster_ids_split, value = TRUE)
			nanopore_ids <- grep("nanopore", cluster_ids_split, value = TRUE)

			# First we take care of the illumina output
			output_row <- vcfs[vcfs$ID == sample(illumina_ids, 1), ]
			output_row$INFO <- paste0(output_row$INFO, ";ClusterIDs=", cluster_ids)
			cat(as.character(output_row), "\n", sep = "\t", file = output_illumina)

			# And now of the nanopore file
			output_row <- vcfs[vcfs$ID == sample(nanopore_ids, 1), ]
			output_row$INFO <- paste0(output_row$INFO, ";ClusterIDs=", cluster_ids)
			cat(as.character(output_row), "\n", sep = "\t", file = output_nanopore)
		}
	}

	return(invisible(NULL))
}

# Lauching the command
select_svs(svmerge_file = "nanopore_svmerged.clustered.vcf", input_illumina = "illumina_svs.vcf", 
	   input_nanopore = "nanopore_svs.vcf", output_illumina = "illumina_merged.vcf", output_nanopore = "nanopore_merged.vcf")

