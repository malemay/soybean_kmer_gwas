# Creating a table that shows which previously discovered signals were found by our analyses

# Loading required libraries
library(GenomicRanges)

# Loading the GRanges object of known signals
# DEPENDENCY: utilities/all_signals.rds
all_signals <- readRDS("utilities/all_signals.rds")

# A vector of program names to query for overlaps with signals
programs <- c("platypus", "vg", "paragraph", "kmers")

# List dependencies for this analysis
# DEPENDENCY: signals found for each program and trait

# Creating a data.frame with one row per locus
signals_table <- unique(as.data.frame(all_signals)[, c("original_trait", "locus", "seqnames",
						       "start", "end", "log_pvalue", "gene_name_v4")])
rownames(signals_table) <- NULL

# Creating columns for the log10 p-values observed for each program
signals_table$platypus <- 0
signals_table$vg <- 0
signals_table$paragraph <- 0
signals_table$kmers <- 0

# For each of the loci, we want to get the most significant p-value detected by each approach
for(i in 1:nrow(signals_table)) {
	message("Processing row ", i)
	i_locus <- signals_table[i, "locus"]
	locus_signals <- all_signals[all_signals$locus == i_locus & all_signals$original_trait == signals_table[i, "original_trait"]]

	for(j in 1:length(locus_signals)) {
		j_signal <- locus_signals[j]

		for(program in programs) {
			psignals <- readRDS(paste0("gwas_results/", program, "/", j_signal$trait, "_signal.rds"))
			psignals <- subsetByOverlaps(psignals, j_signal)
			if(length(psignals)) {
				signals_table[i, program] <- max(signals_table[i, program], psignals$log10_p)
			}
		}
	}
}

# Formatting the p-values columns for output
for(program in programs) {
	signals_table[, program] <- round(signals_table[, program], 1)
}

signals_table[signals_table$platypus == 0, "platypus"] <- "-"
signals_table[signals_table$vg == 0, "vg"] <- "-"
signals_table[signals_table$paragraph == 0, "paragraph"] <- "-"
signals_table[signals_table$kmers == 0, "kmers"] <- "-"

# Formatting the log_pvalue column
signals_table$log_pvalue <- round(signals_table$log_pvalue, 1)

# Formatting the gene name column; also changing the separator
signals_table$gene_name_v4 <- gsub("g", "G", signals_table$gene_name_v4)
signals_table$gene_name_v4 <- gsub(";", ",", signals_table$gene_name_v4)

# The I and B loci are special cases
signals_table[signals_table$locus == "I", "gene_name_v4"] <- "CHS gene cluster"
signals_table[signals_table$locus == "B", "gene_name_v4"] <- "Duplicated HPS genes"

# Formatting the coordinates columns
signals_table$start <- prettyNum(signals_table$start, big.mark = ",")
signals_table$end <- prettyNum(signals_table$end, big.mark = ",")

# Adding a column indicating which study reported the signal
signals_table$study <- NA_character_
signals_table[signals_table$original_trait %in% c("stem_termination", "pubescence_color", "pubescence_form",
						  "pubescence_density", "seed_coat_luster", "pod_color",
						  "flower_color", "seed_coat_color",
						  "hilum_color", "maturity_group"), "study"] <- "Bandillo et al. (2017)"

signals_table[signals_table$original_trait %in% c("oil", "protein"), "study"] <- "Bandillo et al. (2015)"
signals_table[signals_table$original_trait == "corrected_dry_weight", "study"] <- "de Ronne et al. (2022)"
signals_table[signals_table$locus %in% c("stGm11", "pcGm16", "pfGm04", "pfGm15",
					 "sclGm09", "sclGm20", "pdcGm15"), "study"] <- "This study"

# Formatting trait names
signals_table$original_trait <- gsub("_", " ", signals_table$original_trait)
signals_table$original_trait <- sapply(strsplit(signals_table$original_trait, ""), function(x) {x[1] <- toupper(x[1]); paste0(x, collapse = "")})

# Formatting loci names and removing names for the ones for which I added custom names
signals_table[signals_table$locus == "Bq", "locus"] <- "B?"
signals_table$locus <- paste0("\\emph{", signals_table$locus, "}")
signals_table[grepl("Gm[0-9]{2}", signals_table$locus), "locus"] <- "-"

# Creating a single column with the position on the reference
signals_table$position <- paste0(signals_table$seqnames, ":", signals_table$start, "-", signals_table$end)

# Changing the column names
colnames(signals_table) <- c("Trait", "Locus", "Chromosome", "Start", "End",
			     "Pvalues", "Gene", "Platypus", "Vg",
			     "Paragraph", "Kmers", "Study", "Position")

# Sorting the columns for display
signals_table <- signals_table[order(signals_table$Trait, signals_table$Pvalues, decreasing = c(FALSE, TRUE), method = "radix"), ]

# Writing this table to a csv file for use in the supplemental data file
write.table(signals_table, file = "tables/signals_table.csv",
	    quote = FALSE, sep = ";", row.names = FALSE, col.names = TRUE)

