# Creating a table that shows the performance of various approaches in
# identifying causal genes/variants at loci for which the gene is known

# Loading required libraries
suppressMessages(library(GenomicRanges))

# Loading the GRanges object of known signals
# DEPENDENCY: utilities/all_signals.rds
all_signals <- readRDS("utilities/all_signals.rds")

# Keeping only the information that we need
loci_table <- unique(as.data.frame(all_signals)[, c("original_trait", "locus", "gene_name_v4")])

# Ordering the rows of the table in the order in which they are presented in the manuscript
loci_table$original_trait <- factor(loci_table$original_trait,
				    levels = c("flower_color", "pubescence_color", "seed_coat_color", "stem_termination",
					       "hilum_color", "pod_color", "pubescence_form", "pubescence_density",
					       "seed_coat_luster", "maturity_group", "corrected_dry_weight",
					       "oil", "protein"))

loci_table <- loci_table[order(loci_table$original_trait), ]
loci_table$original_trait <- as.character(loci_table$original_trait)

# Adding an ID for each of the trait/loci combinations
loci_table$id <- paste0(loci_table$original_trait, "_", loci_table$locus)

# Keeping only loci for which there is a known cloned gene
loci_table <- loci_table[!is.na(loci_table$gene_name_v4), ]

# Removing a few cases for which we only have candidates
loci_table <- loci_table[!loci_table$locus %in% c("Pa1", "cdwGm15"), ]

# Setting the GWAS to use for various loci using a lookup table
gwas_lookup <- c("stem_termination_Dt2" = "stem_termination_all_Dt2",
		 "stem_termination_Dt1" = "stem_termination_all_Dt1",
		 "stem_termination_E3" = "stem_termination_all_E3",
		 "pubescence_color_Td" = "pubescence_color_nogray_Td",
		 "pubescence_color_T" = "pubescence_color_all_T",
		 "protein_proteinGm20" = "protein_proteinGm20",
		 "protein_proteinGm15" = "protein_proteinGm20",
		 "oil_oilGm20" = "oil_oilGm20",
		 "oil_oilGm15" = "oil_oilGm15",
		 "pubescence_density_Pd1" = "pubescence_density_Pd1",
		 "pubescence_density_Ps" = "pubescence_density_Ps",
		 "pubescence_density_P1" = "pubescence_density_P1",
		 "seed_coat_luster_I" = "seed_coat_luster_dullshiny_I",
		 "seed_coat_luster_B" = "seed_coat_luster_dullshiny_B",
		 "flower_color_W1" = "flower_color_W1",
		 "seed_coat_color_G" = "seed_coat_color_greenyellow_G",
		 "seed_coat_color_T" = "seed_coat_color_all_T",
		 "seed_coat_color_I" = "seed_coat_color_all_I",
		 "seed_coat_color_R" = "seed_coat_color_all_R",
		 "hilum_color_T" = "hilum_color_all_T",
		 "hilum_color_I" = "hilum_color_all_I",
		 "hilum_color_R" = "hilum_color_blackbrown_R",
		 "hilum_color_W1" = "hilum_color_all_W1",
		 "maturity_group_E1" = "maturity_group_E1",
		 "maturity_group_E2" = "maturity_group_E2",
		 "maturity_group_E3" = "maturity_group_E3",
		 "maturity_group_E4" = "maturity_group_E4")


# Creating a list of GRanges objects comprising the genes or causal variants
# DEPENDENCY: refgenome/gmax_v4_genes.rds
genes <- readRDS("refgenome/gmax_v4_genes.rds")

gene_list <- list()

for(i in 1:nrow(loci_table)) {
	i_locus <- loci_table[i, "id"]
	i_gene <- loci_table[i, "gene_name_v4"]
	gene_list[[i_locus]] <- if(!grepl(";", i_gene)) genes[i_gene] else NA
}

# Reading the GRanges objects of the inversion at locus I and CNV at loci B and Ps
# DEPENDENCY : cnv_analysis/i_cnv_range.rds
# DEPENDENCY : cnv_analysis/hps_cnv_range.rds
# DEPENDENCY : cnv_analysis/ps_cnv_range.rds
i_range <- readRDS("cnv_analysis/i_cnv_range.rds")
b_range <- readRDS("cnv_analysis/hps_cnv_range.rds")
ps_range <- readRDS("cnv_analysis/ps_cnv_range.rds")

# Adding those ranges to the gene list
gene_list[["seed_coat_luster_I"]] <- i_range
gene_list[["seed_coat_luster_B"]] <- b_range
gene_list[["seed_coat_color_I"]] <- i_range
gene_list[["hilum_color_I"]] <- i_range
gene_list[["pubescence_density_Ps"]] <- ps_range

# For each range, we find out whether it overlaps various features
# from each approach. We return a value that depends on how
# well a given approach succeded: 
# - 3 if the gene overlaps the top region
# - 2 if the gene overlaps a signal
# - 1 otherwise
check_overlap <- function(locus, approach, gene_list, trait_lookup) {
	# First we extract the GRanges corresponding to the gene
	gene <- gene_list[[locus]]
	trait <- trait_lookup[locus]

	# First we check whether the gene overlaps the top region
	top_markers <- readRDS(paste0("gwas_results/", approach, "/", trait, "_top_markers.rds"))
	if(overlapsAny(gene, top_markers)) return(3)

	# Next we check if there is any signal overlapping the gene/variant
	signals <- readRDS(paste0("gwas_results/", approach, "/", trait, "_signal_locus.rds"))
	if(overlapsAny(gene, signals)) return(2)

	return(1)
}

# Filling columns for each approach depending on their performance
loci_table$platypus <- sapply(loci_table$id, check_overlap, approach = "platypus", gene_list = gene_list, trait_lookup = gwas_lookup)
loci_table$paragraph <- sapply(loci_table$id, check_overlap, approach = "paragraph", gene_list = gene_list, trait_lookup = gwas_lookup)
loci_table$kmers <- sapply(loci_table$id, check_overlap, approach = "kmers", gene_list = gene_list, trait_lookup = gwas_lookup)

# Manually adjusting the cases where the most associated variant was the causal variant
loci_table[loci_table$id == "flower_color_W1", c("paragraph", "kmers")] <- 4
loci_table[loci_table$id == "pubescence_color_Td", "kmers"] <- 4
loci_table[loci_table$id == "seed_coat_color_G", "kmers"] <- 4
loci_table[loci_table$id == "pubescence_density_Ps", "kmers"] <- 4
# And also a single case where a gene overlap is falsely reported
loci_table[loci_table$id == "seed_coat_color_R", "kmers"] <- 1


# Converting the numbers to asteriks for display in a table
significance_lookup <- c("-", "*", "**", "***")
for(i in c("platypus", "paragraph", "kmers")) {
	loci_table[[i]] <- significance_lookup[loci_table[[i]]]
}

# Formatting the table for output
rownames(loci_table) <- NULL

loci_table$original_trait <- gsub("_", " ", loci_table$original_trait)
loci_table$original_trait <- sapply(strsplit(loci_table$original_trait, ""), function(x) {x[1] <- toupper(x[1]); paste0(x, collapse = "")})
loci_table[loci_table$original_trait == "Stem termination", "original_trait"] <- "Stem termination type"

# The I and B loci are special cases
loci_table[loci_table$locus == "I", "gene_name_v4"] <- "CHS gene cluster"
loci_table[loci_table$locus == "B", "gene_name_v4"] <- "Duplicated HPS genes"

# Renaming the protein loci
loci_table$locus <- sub("protein", "pro", loci_table$locus)

# Setting the loci to italics
loci_table$locus <- paste0("\\textit{", loci_table$locus, "}")

# Formatting the gene names with an initial capital
loci_table$gene_name_v4 <- sub("^g", "G", loci_table$gene_name_v4)

# Adding the type of causal variant using a lookup table
causal_lookup <- c("flower_color_W1" = "65-bp deletion \\citep{zabala2007}",
		   "pubescence_color_T" = "1-bp deletion \\citep{zabala2003}",
		   "pubescence_color_Td" = "SNP \\citep{yan2020}",
		   "seed_coat_color_I" = "10.9-kb inversion/duplication \\citep{tuteja2008}",
		   "seed_coat_color_G" = "SNP \\citep{wang2018}",
		   "seed_coat_color_T" = "1-bp deletion \\citep{zabala2003}",
		   "seed_coat_color_R" = "SNPs and indels \\citep{gillman2011}",
		   "stem_termination_Dt1" = "SNPs \\citep{liu2010}",
		   "stem_termination_Dt2" = "SNPs \\citep{ping2014}",
		   "stem_termination_E3" = "Various types \\citep{tardivel2019}",
		   "hilum_color_T" = "1-bp deletion \\citep{zabala2003}",
		   "hilum_color_I" = "10.9-kb inversion/duplication \\citep{tuteja2008}",
		   "hilum_color_R" = "SNPs and indels \\citep{gillman2011}",
		   "hilum_color_W1" = "65-bp deletion \\citep{zabala2007}",
		   "pubescence_density_Ps" = "25.6-kb CNV \\citep{liu2020ps}",
		   "pubescence_density_Pd1" = "SNP \\citep{liu2020ps}",
		   "pubescence_density_P1" = "SNP \\citep{liu2020ps}",
		   "seed_coat_luster_B" = "CNV (several kb) \\citep{gijzen2006}",
		   "seed_coat_luster_I" = "10.9-kb inversion/duplication \\citep{tuteja2008}",
		   "maturity_group_E1" = "Various types \\citep{tardivel2019}",
		   "maturity_group_E2" = "Various types \\citep{tardivel2019}",
		   "maturity_group_E3" = "Various types \\citep{tardivel2019}",
		   "maturity_group_E4" = "Various types \\citep{tardivel2019}",
		   "oil_oilGm20" = "321-bp deletion \\citep{fliege2022}",
		   "oil_oilGm15" = "SNPs and indels \\citep{zhang2020}",
		   "protein_proteinGm20" = "321-bp deletion \\citep{fliege2022}",
		   "protein_proteinGm15" = "SNPs and indels \\citep{zhang2020}")

# Add the column for the causal variant
loci_table$causal <- causal_lookup[loci_table$id]

# Adding a note for loci that are not discussed in the main text
nomain <- c("seed_coat_color_T", "seed_coat_color_R", "stem_termination_Dt2",
	    "stem_termination_E3", "hilum_color_W1", "pubescence_density_Pd1",
	    "pubescence_density_P1", "seed_coat_luster_I", "maturity_group_E1",
	    "maturity_group_E2", "maturity_group_E3", "maturity_group_E4",
	    "oil_oilGm20", "oil_oilGm15", "protein_proteinGm20", "protein_proteinGm15",
	    "hilum_color_T", "hilum_color_I")

loci_table$original_trait <- paste0(loci_table$original_trait, ifelse(loci_table$id %in% nomain, "\\tnote{c}", ""))

# Reordering the table based on the causal_lookup vector
loci_table$id <- factor(loci_table$id, levels = names(causal_lookup))
loci_table <- loci_table[order(loci_table$id), ]

loci_table$id <- NULL

# Formatting the columns
colnames(loci_table) <- c("Trait", "Locus", "Gene", "platypus", "paragraph", "kmers", "Causal")

# Writing the table to file
write.table(loci_table, file = "tables/loci_table.csv", sep = ",", row.names = FALSE,
	    col.names = TRUE, quote = FALSE)

