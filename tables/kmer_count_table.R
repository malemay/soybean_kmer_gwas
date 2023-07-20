# Generating a table with the number of significant k-mers found
# for each GWAS analysis
library(GenomicRanges)

# DEPENDENCY: utilities/all_signals.rds
all_signals <- readRDS("utilities/all_signals.rds")

# Initializing the output table
output_table <- unique(as.data.frame(all_signals)[, c("trait", "original_trait")])
rownames(output_table) <- NULL

# Creating a lookup table from id to the name to display
id_lookup <- c("stem_termination_all" = "Stem termination type (GWAS \\#1)",
	       "stem_termination_sn" = "Stem termination type (GWAS \\#2)",
	       "pubescence_color_all" = "Pubescence color (GWAS \\#1)",
	       "pubescence_color_nogray" = "Pubescence color (GWAS \\#2)",
	       "protein" = "Protein",
	       "pubescence_form_all" = "Pubescence form (GWAS \\#1)",
	       "pubescence_form_noerect" = "Pubescence form (GWAS \\#2)",
	       "oil" = "Oil",
	       "pubescence_density" = "Pubescence density",
	       "seed_coat_luster_all" = "Seed coat luster (GWAS \\#1)",
	       "seed_coat_luster_nointermediate" = "Seed coat luster (GWAS \\#2)",
	       "seed_coat_luster_dullshiny"      = "Seed coat luster (GWAS \\#3)",
	       "pod_color_all"   = "Pod color (GWAS  \\#1)",
	       "pod_color_blbr"  = "Pod color (GWAS  \\#2)",
	       "flower_color" = "Flower color",
	       "seed_coat_color_all" = "Seed coat color (GWAS \\#1)",
	       "seed_coat_color_greenyellow" = "Seed coat color (GWAS \\#2)",
	       "hilum_color_all"         = "Hilum color (GWAS \\#1)",
	       "hilum_color_blackbrown"  = "Hilum color (GWAS \\#2)",
	       "hilum_color_rbr"         = "Hilum color (GWAS \\#3)",
	       "maturity_group"          = "Maturity group",
	       "corrected_dry_weight" = "Resistance to \\textit{P. sojae}")

# Adding that column
output_table$gwas <- id_lookup[output_table$trait]

output_table$original_trait <- factor(output_table$original_trait,
				       levels = c("flower_color", "pubescence_color", "seed_coat_color", "stem_termination",
						  "hilum_color", "pod_color", "pubescence_form", "pubescence_density",
						  "seed_coat_luster", "maturity_group", "corrected_dry_weight",
						  "oil", "protein"))

output_table <- output_table[order(output_table$original_trait, output_table$gwas), ]
output_table$original_trait <- as.character(output_table$original_trait)

# Adding a colum with the number of raw significantr k-mers found
output_table$raw <- 0
output_table$filtered <- 0

# DEPENDENCY: k-mer GWAS analysis
# DEPENDENCY: k-mer signals
for(i in 1:nrow(output_table)) {
output_table[i, "raw"] <- length(readLines(paste0("gwas_results/kmer_data/", output_table[i, "trait"], "/kmers/pass_threshold_5per"))) - 1
output_table[i, "filtered"] <- length(readRDS(paste0("gwas_results/kmers/", output_table[i, "trait"], "_kmer_positions.rds")))
}

output_table$nsignals <- 0

# Finally we need to count the number of k-mers that overlap known signals
for(i in 1:nrow(output_table)) {
	i_signals <- readRDS(paste0("gwas_results/kmers/", output_table[i, "trait"], "_signal.rds"))
	i_subset <- subsetByOverlaps(i_signals, all_signals[all_signals$trait == output_table[i, "trait"]])
	output_table[i, "nsignals"] <- sum(i_subset$n_markers)
}

output_table$percent <- round(output_table$nsignals / output_table$filtered * 100, 1)
output_table$percent <- ifelse(is.finite(output_table$percent), as.character(output_table$percent), "NA")

# Formatting the numbers prior to outût
output_table$raw <- prettyNum(output_table$raw, big.mark = ",")
output_table$filtered <- prettyNum(output_table$filtered, big.mark = ",")
output_table$nsignals <- prettyNum(output_table$nsignals, big.mark = ",")

# Write output
write.table(output_table, file = "tables/kmer_count_table.csv", sep = ";", quote = FALSE, col.names = TRUE, row.names = FALSE)

