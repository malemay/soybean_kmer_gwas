# A function that takes a SVABA vcf as input and outputs a vcf with SVs that are recoded as
#  DEL, DUP or INV depending of the characteristics of the breakpoints. All SVABA SVs are
#  initially encoded as BND, so it is necessary to undergo these steps to convert from
#  BND to a format that genotypers can understand

svaba_classifier <- function(vcf_file, output_file) {
  
  # Opening a connection for the output file
  output <- file(output_file, open = "w")
  on.exit(close(output, type = "w"), add = TRUE)
  
  # We first read the whole file in memory (since it is not too big)
  records <- scan(vcf_file, what = character(), sep = "\n", quiet = TRUE)
  
  # Now we iterate over all the lines in the file
  for(i in 1:length(records)) {
    
    # Assigning the current line to the variable "line" 
    line <- records[i]
    
    # Header lines are just printed as is
    if(grepl("^#", line)) {

	    # We add a metainfo line for MATEALT when we reach the INFO lines
	    if(grepl("^##INFO=<ID=SVTYPE,", line)) {
		    cat('##INFO=<ID=MATE_ALT,Number=1,Type=String,Description="ALT allele of the mate (for debugging purposes)">\n', file = output)
	    }

      cat(line, "\n", sep = "", file = output)
      next
    }
    
    # Splitting the line into its tab-separated fields
    vcf_fields <- strsplit(line, "\t")[[1]]
    
    # Fields are named for easier data extraction
    stopifnot(length(vcf_fields) == 10) # Only supports one sample
    names(vcf_fields) <- c("CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", "SAMPLE")
    
    # We only take care of the first part of the breakpoint pair
    if(!grepl(":1", vcf_fields["ID"])) next
      
    # It is a deletion if it matches a pattern of the kind N[GmXY:123456[
    if (grepl("^[ATGCN]+\\[", vcf_fields["ALT"])) {
      svtype <- "DEL"
      vcf_fields["INFO"] <-  sub("SVTYPE=[A-Z]*([;]?)", "SVTYPE=DEL\\1", vcf_fields["INFO"])
      
      # It is an insertion/duplication if it matches a pattern of the kind ]GmXY:123456]N
    } else if (grepl("\\][ATGCN]+$", vcf_fields["ALT"])) {
      svtype <- "DUP"
      vcf_fields["INFO"] <-  sub("SVTYPE=[A-Z]*([;]?)", "SVTYPE=DUP\\1", vcf_fields["INFO"])
      
      # It is an inversion if it matches either [GmXY:123456[N or N]GmXY:123456]
    } else if (grepl("\\[[ATGCN]+$", vcf_fields["ALT"]) || grepl("^[ATGCN]+\\]", vcf_fields["ALT"])) {
      svtype <- "INV"
      vcf_fields["INFO"] <-  sub("SVTYPE=[A-Z]*([;]?)", "SVTYPE=INV\\1", vcf_fields["INFO"])
      
      # Otherwise it is an unknown SVTYPE
    } else {
      svtype <- "UNKNOWN"
      vcf_fields["INFO"] <-  sub("SVTYPE=[A-Z]*([;]?)", "SVTYPE=UNKNOWN\\1", vcf_fields["INFO"])
    }
    
    # We need this breakpoint's mate to add its ALT allele in the INFO field
    mate_id <- sub(".*MATEID=([0-9]+:[0-9]);.*", "\\1", vcf_fields["INFO"])
    
    # Then we loop from line i+1 until we find the mate
    found_mate <- FALSE
    j <- i + 1
    
    while(!found_mate) {
      # We extract the line and its corresponding fields
      line2 <- records[j]
      line2_fields <- strsplit(line2, "\t")[[1]]
      names(line2_fields) <- c("CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", "SAMPLE")
      
      # We check if the MATEID corresponds
      # If it does we add it to the info field of the i line we are looping on
      if(line2_fields["ID"] == mate_id) {
        found_mate <- TRUE
        vcf_fields["INFO"] <- paste0(paste0("MATE_ALT=", line2_fields["ALT"], ";"), vcf_fields["INFO"])
      } else {
        # Otherwise we go to the next j line
        j <- j + 1
        if(j > length(records)) {
		message("Mate not found for ", vcf_fields["ID"])
		break
	}
      }
    }
    
    # The modified i line is output to the output file
    cat(paste(vcf_fields, collapse = "\t"), "\n", sep = "", file = output)
  }
  
  invisible(NULL)
}




# A function to get the reverse complement of a sequence 
# (useful for determining the sequence of inversions)
rev_comp <- function(seq) {
  # A lookup table for getting the complement
  comp_table <- c("A" = "T",
                  "T" = "A",
                  "G" = "C",
                  "C" = "G",
                  "N" = "N")
  
  paste0(rev(comp_table[strsplit(seq, "")[[1]]]), collapse = "")
}




# Now writing a function that will take the files output by the function svaba_classifier
#  as input and output a vcf where the full sequence of the alternative allele is specified
svaba_converter <- function(vcf_file, output_file, refgenome, contig_pattern, max_span = 500000) {
  
  # Opening a connection for the output file
  output <- file(output_file, open = "w")
  on.exit(close(output), add = TRUE)
  
  # Creating a scalar value for the number of variants skipped because of a too large span
  n_skipped <- 0
  
  # We first read the whole file in memory (since it is not too big)
  records <- scan(vcf_file, what = character(), sep = "\n", quiet = TRUE)
  
  # Now we iterate over all the lines in the file
  for(i in 1:length(records)) {
    
    # Assigning the current line to the variable "line" 
    line <- records[i]
    
    # Header lines are just printed as is
    if(grepl("^#", line)) {

	    # We add a metainfo line for ORALT when we reach the INFO lines
	    if(grepl("^##INFO=<ID=SVTYPE,", line)) {
		    cat('##INFO=<ID=ORALT,Number=1,Type=String,Description="ALT allele before modification by svaba_converter">\n', file = output)
	    }

      cat(line, "\n", sep = "", file = output)
      next
    }
    
    # Splitting the line into its tab-separated fields
    vcf_fields <- strsplit(line, "\t")[[1]]
    
    # Fields are named for easier data extraction
    stopifnot(length(vcf_fields) == 10) # Only supports one sample
    names(vcf_fields) <- c("CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", "SAMPLE")
    
    # Extracting the SPAN numeric value
    span <- as.numeric(sub(".*SPAN=(-?[0-9]+);.*", "\\1", vcf_fields["INFO"]))
    # This vcf line is skipped if the variant span is >= max_span or if span is < 0 (interchromosomal)
    if(span >= max_span || span < 0) {
      n_skipped <- n_skipped + 1
      message("Line ", i, " skipped, span = ", span)
      next
    }

    # We do not process variants whose position is 0
    if(vcf_fields["POS"] <= 0) {
      n_skipped <- n_skipped + 1
      message("Line ", i, " skipped, POS = ", vcf_fields["POS"])
      next
    }
    
    # We also remove variants that don't match contig_pattern
    if(!(grepl(contig_pattern, vcf_fields["CHROM"]))) {
      n_skipped <- n_skipped + 1
      message("Line ", i, " skipped, contig = ", vcf_fields["CHROM"])
      next
    }
    
    # Adding the ALT allele to the INFO field before modifying it (for debugging purposes)
    vcf_fields["INFO"] <- paste0("ORALT=", vcf_fields["ALT"], ";", vcf_fields["INFO"])
    
    # We will also need the sequence that is inserted at the breakpoint (if present)
    if(grepl("INSERTION", vcf_fields["INFO"])) {
      insertion <- sub(".*INSERTION=([ATGCN]+);.*", "\\1", vcf_fields["INFO"])
    } else {
      # If there is no insertion, we will use a placeholder for convenience
      insertion <- ""
    }
    
    # The case for deletions; deletions are of the type N[GmXY:123456[
    if (grepl("SVTYPE=DEL", vcf_fields["INFO"])) {
      # We first extract the initial nucleotide
      nuc <- sub("^([ATGCN]+)\\[.*", "\\1",  vcf_fields["ALT"])
      # This nucleotide should always be equal to the REF
      stopifnot(nuc == vcf_fields["REF"])
      
      ## Then we extract the chromosome on which the sequence is located
      chr <- sub("^[ATGCN]+\\[([a-zA-Z]+[0-9]+):[0-9]+\\[", "\\1", vcf_fields["ALT"])
      
      ## And the position of the breakpoint on that chromosome
      bp  <- as.numeric(sub("^[ATGCN]+\\[[a-zA-Z]+[0-9]+:([0-9]+)\\[", "\\1", vcf_fields["ALT"]))
        
      # For deletions, we actually need to transform the REF allele to encode the deletion
      # We need to paste REF + the sequence that has been deleted from pos + 1 until bp - 1
      # We get the sequence corresponding to this region from the reference genome
      #ref_sequence <- substr(get(chr, envir = as.environment("package:Gmaxv4")),
       #                      as.numeric(vcf_fields["POS"]) + 1, 
        #                     bp - 1)
      ref_sequence <- Rsamtools::scanFa(refgenome,
					GRanges(seqnames = chr, 
						IRanges(start = as.numeric(vcf_fields["POS"]) + 1, 
							end = bp - 1)))
      ref_sequence <- as.character(ref_sequence)

      # We modify the reference allele accordingly :
      vcf_fields["REF"] <- paste0(nuc, ref_sequence)
      
      # We also modify the alternate allele by adding the inserted sequence
      vcf_fields["ALT"] <- paste0(nuc, insertion)
      

     # The case for duplications; follow the pattern ]GmXY:123456]N
    } else if (grepl("SVTYPE=DUP", vcf_fields["INFO"])) {

      # We first extract the initial nucleotide
      nuc <- sub(".*\\]([ATGCN]+)$", "\\1",  vcf_fields["ALT"])
      # This nucleotide should always be equal to the REF
      stopifnot(nuc == vcf_fields["REF"])
      
      ## Then we extract the chromosome on which the sequence is located
      chr <- sub("^\\]([a-zA-Z]+[0-9]+):[0-9]+\\][ATGCN]+", "\\1", vcf_fields["ALT"])
      
      ## And the position of the breakpoint on that chromosome
      bp  <- as.numeric(sub("^\\][a-zA-Z]+[0-9]+:([0-9]+)\\][ATGCN]+", "\\1", vcf_fields["ALT"]))
      
      # We need the reference sequence from POS to bp as all these nucleotides are duplicated
      #ref_sequence <- substr(get(chr, envir = as.environment("package:Gmaxv4")),
       #                      as.numeric(vcf_fields["POS"]), 
        #                     bp)
      
      ref_sequence <- Rsamtools::scanFa(refgenome,
					GRanges(seqnames = chr, 
						IRanges(start = as.numeric(vcf_fields["POS"]), 
							end = bp)))
      ref_sequence <- as.character(ref_sequence)


      # The reference sequence does not need to be updated
      
      # We modify the alternate allele by adding the duplication and inserted sequence
      vcf_fields["ALT"] <- paste0(ref_sequence, insertion, nuc)
      
      
      # The case for inversions ; inversions match either [GmXY:123456[N or N]GmXY:123456]
      # Note : I am not entirely certain about the treatment of inversions ;
      # I was largely inspired from the section 5.4.7 of the VCf spec 4.2 to figure this out
    } else if (grepl("SVTYPE=INV", vcf_fields["INFO"])) {

      if(grepl("^[ATGCN]+\\]", vcf_fields["ALT"])) {
        # Let's start by extracting the initial nucleotide, chromosome and breakpoint position
        ## We first extract the initial nucleotide
        nuc <- sub("^([ATGCN]+)\\][a-zA-Z]+[0-9]+:[0-9]+\\]", "\\1", vcf_fields["ALT"])
        # This nucleotide should always be equal to the REF
        stopifnot(nuc == vcf_fields["REF"])
        ## Then we extract the chromosome on which the sequence is located
        chr <- sub("^[ATGCN]+\\]([a-zA-Z]+[0-9]+):[0-9]+\\]", "\\1", vcf_fields["ALT"])
        
        ## And the position of the breakpoint on that chromosome
        bp  <- as.numeric(sub("^[ATGCN]+\\][a-zA-Z]+[0-9]+:([0-9]+)\\]", "\\1", vcf_fields["ALT"]))
        
        # We need to identify to start and end positions of the inverted sequence
        # In this case it is POS + 1 and bp
        #ref_sequence <- substr(get(chr, envir = as.environment("package:Gmaxv4")),
         #                      as.numeric(vcf_fields["POS"]) + 1, 
          #                     bp)
	ref_sequence <- Rsamtools::scanFa(refgenome,
					  GRanges(seqnames = chr, 
						  IRanges(start = as.numeric(vcf_fields["POS"]) + 1, 
							  end = bp)))
	ref_sequence <- as.character(ref_sequence)
        
        # We update the reference allele with the sequence over the full inverted section
        vcf_fields["REF"] <- paste0(nuc, ref_sequence)
        
        # We update the alternate allele with the full reverse-complemented sequence of the inverted section, and add the insertion
        vcf_fields["ALT"] <- paste0(nuc, insertion, rev_comp(ref_sequence))
        
        
        
        
      } else if(grepl("\\[[ATGCN]+$", vcf_fields["ALT"])) {
        # In this case we need the chromsome and bp position but not the initial nucleotide
        ## Then we extract the chromosome on which the sequence is located
        chr <- sub("^\\[([a-zA-Z]+[0-9]+):[0-9]+\\[[ATGCN]+", "\\1", vcf_fields["ALT"])
        
        ## And the position of the breakpoint on that chromosome
        bp  <- as.numeric(sub("^\\[[a-zA-Z]+[0-9]+:([0-9]+)\\[[ATGCN]+", "\\1", vcf_fields["ALT"]))
        
        # The reference position should be POS - 1
        # The position in the vcf will be updated accordingly
        vcf_fields["POS"] <- as.character(as.numeric(vcf_fields["POS"]) - 1)
        # Thus nuc will be extracted from the reference sequence
      #  nuc <- substr(get(chr, envir = as.environment("package:Gmaxv4")),
       #               as.numeric(vcf_fields["POS"]), 
        #              as.numeric(vcf_fields["POS"]))
	nuc <- Rsamtools::scanFa(refgenome,
				 GRanges(seqnames = chr, 
					 IRanges(start = as.numeric(vcf_fields["POS"]), 
						 end =   as.numeric(vcf_fields["POS"]))))
	nuc <- as.character(nuc)
        
        # We get the reference sequence over the area affected by the inversion
        #ref_sequence <- substr(get(chr, envir = as.environment("package:Gmaxv4")),
         #                      as.numeric(vcf_fields["POS"]) + 1, 
        #                       bp - 1)
	ref_sequence <- Rsamtools::scanFa(refgenome,
					  GRanges(seqnames = chr, 
						  IRanges(start = as.numeric(vcf_fields["POS"]) + 1, 
							  end =   bp - 1)))
	ref_sequence <- as.character(ref_sequence)
        
        # The reference allele is updated
        vcf_fields["REF"] <- paste0(nuc, ref_sequence)
        
        # We update the alternate allele with the full reverse-complemented sequence of the inverted section, and add the insertion
        vcf_fields["ALT"] <- paste0(nuc, rev_comp(ref_sequence), insertion)
        
      } else {
        stop("Unknown inversion type")
      }
      
      
      # Otherwise it is an unknown SVTYPE
    } else {
      stop("SVTYPE not recognized")
    }
    
    # The modified i line is output to the output file
    cat(paste(vcf_fields, collapse = "\t"), "\n", sep = "", file = output)
  }
  
  invisible(NULL)
}

