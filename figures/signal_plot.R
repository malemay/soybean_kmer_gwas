# Loading the required packages
suppressMessages(library(grid))
suppressMessages(library(gwask))
suppressMessages(library(GenomicRanges))

# Reading the ID being analyzed from the command line
id <- commandArgs(trailingOnly = TRUE)[1]

# Reading the information on genes and transcripts
genes <- readRDS("refgenome/gmax_v4_genes.rds")
transcripts <- readRDS("refgenome/gmax_v4_transcripts.rds")
cds <- readRDS("refgenome/gmax_v4_cds.rds")
exons <- readRDS("refgenome/gmax_v4_exons.rds")

# Getting the ID of the gene associated with the trait and the name of the trait from the all_signals.rds file
# DEPENDENCY: utilities/all_signals.rds
target_signal <- readRDS("utilities/all_signals.rds")
target_signal <- target_signal[id]

gene_name <- target_signal$gene_name_v4

# Handling the special case where the gene_name value contains more than one gene or no gene at all
if(!is.na(gene_name) && grepl(";", gene_name)) {
	gene_names <- strsplit(gene_name, ";")[[1]]
	gene <- genes[gene_names]
	gene <- reduce(gene, min.gapwidth = 10^6, ignore.strand = TRUE)
} else if (!is.na(gene_name)) {
	gene <- genes[gene_name]
} else {
	gene <- NULL
}

# DEPENDENCY: Paragraph signal subplot
# DEPENDENCY: Platypus signal subplot
# DEPENDENCY: kmers signal subplot
pvalue_grobs <- list(readRDS(paste0("figures/grobs/platypus_", id, "_signal.rds")),
		     readRDS(paste0("figures/grobs/paragraph_", id, "_signal.rds")),
		     readRDS(paste0("figures/grobs/kmers_", id, "_signal.rds")))

output_grob <- pvalue_tx_grob(pvalue_grobs = pvalue_grobs,
			      xrange = NULL,
			      xchrom = seqnames(target_signal),
			      genes = genes,
			      transcripts = transcripts,
			      exons = exons,
			      cds = cds,
			      strand_colors = c("dodgerblue", "dodgerblue"),
			      first_tx_only = TRUE,
			      margins = c(3.5, 3.6, 2.5, 1.1))

# Outputting to a PNG file
png(paste0("figures/", id, "_signal.png"), width = 9, height = 10, units = "in", res = 100)

# Resetting the viewport
grid.newpage()

grid.draw(output_grob)

# Add labels for the panels
grid.text(c("(a)", "(b) SNPs/indels", "(c) SVs"),
	  x = unit(0.09, "npc"),
	  y = unit(c(0.993, 0.8725, 0.57), "npc"),
	  hjust = 0,
	  gp = gpar(fontsize = 14))
grid.text(expression(paste("(d) ", italic(k)-mers)),
	  x = unit(0.09, "npc"),
	  y = unit(0.2675, "npc"),
	  hjust = 0,
	  gp = gpar(fontsize = 14))

dev.off()

