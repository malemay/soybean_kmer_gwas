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

### ----------
# Getting the chromosome of the gene associated with the trait from the all_signals.rds file
# DEPENDENCY: utilities/all_signals.rds
target_signal <- readRDS("utilities/all_signals.rds")
xchrom <- seqnames(target_signal[id])

# DEPENDENCY: Paragraph gene subplot
# DEPENDENCY: vg gene subplot
# DEPENDENCY: Platypus gene subplot
# DEPENDENCY: kmers gene subplot
pvalue_grobs <- list(readRDS(paste0("figures/grobs/platypus_", id, "_gene.rds")),
		     readRDS(paste0("figures/grobs/vg_", id, "_gene.rds")),
		     readRDS(paste0("figures/grobs/paragraph_", id, "_gene.rds")),
		     readRDS(paste0("figures/grobs/kmers_", id, "_gene.rds")))

output_grob <- pvalue_tx_grob(pvalue_grobs = pvalue_grobs,
			      xchrom = xchrom,
			      genes = genes,
			      transcripts = transcripts,
			      exons = exons,
			      cds = cds,
			      draw_arrows = TRUE,
			      first_tx_only = TRUE,
			      strand_colors = c("dodgerblue", "dodgerblue"),
			      margins = c(3.5, 3.6, 2.5, 1.1))

# Outputting to a PNG file
png(paste0("figures/", id, "_gene.png"), width = 9, height = 10, units = "in", res = 100)

# Resetting the viewport
grid.newpage()

grid.draw(output_grob)

# Add labels for the panels
grid.text(c("(a)", "(b) Platypus", "(c) vg", "(d) Paragraph"),
	  x = unit(0.09, "npc"),
	  y = unit(c(0.993, 0.8725, 0.645, 0.4175), "npc"),
	  hjust = 0,
	  gp = gpar(fontsize = 14))
grid.text(expression(paste("(e) ", italic(k)-mers)),
	  x = unit(0.09, "npc"),
	  y = unit(0.19, "npc"),
	  hjust = 0,
	  gp = gpar(fontsize = 14))

dev.off()

