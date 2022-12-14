# Loading the required packages
suppressMessages(library(grid))
suppressMessages(library(gwask))
suppressMessages(library(GenomicRanges))


# Reading in the command line arguments; the first one is the trait, the second is the genotyping program
trait <- commandArgs(trailingOnly = TRUE)[1]
program <- commandArgs(trailingOnly = TRUE)[2]

# Reading the reference signals for this trait
# DEPENDENCY: utilities/all_signals.rds
signals <- readRDS("utilities/all_signals.rds")
# converting the name of the Bq locus back to B? for plotting purposes
signals[signals$locus == "Bq"]$locus <- "B?"
signals <- signals[signals$trait == trait]
if(!length(signals)) signals <- NULL

# Reading the signals that were found for this trait and program
# DEPENDENCY: trait-wise signals
found_signals <- readRDS(paste0("gwas_results/", program, "/", trait, "_signal.rds"))
if(!length(signals)) found_signals <- NULL

# It is a special case if we are analyzing k-mers
threshold <- -log10(as.numeric(readLines(paste0("gwas_results/", program, "/", trait, "_threshold_5per.txt"))))

# Loading the results from the rds file
gwas_results <- readRDS(paste0("gwas_results/", program, "/", trait, "_gwas.rds"))

# If the program is the k-mers, we filter out variants that are located on unachored scaffolds
if(program == "kmers") {
	GenomeInfoDb::seqlevels(gwas_results, pruning.mode = "coarse") <- grep("^Gm[0-9]{2}$", GenomeInfoDb::seqlevels(gwas_results), value = TRUE)
}
	
# A vector of traits for which the signal labels should be plotted on the left of the signal itself
left_labels <- c("seed_coat_luster_all", "seed_coat_luster_nointermediate", "seed_coat_luster_dullshiny", "oil", "protein")
# Plotting the results using the gwask::manhattan_plot function
gwas_plot <- manhattanGrob(gwas_results,
			   threshold = threshold,
			   ref_signals = signals,
			   new_signals = found_signals,
			   label_offset = if(trait %in% left_labels) -10^7 else 10^7,
			   label_hjust  = if(trait %in% left_labels) 1 else 0,
			   point_colors = c("blue", "skyblue3"),
			   numeric_chrom = TRUE,
			   margins = c(5.1, 4.1, 0.1, 0.1))

# Outputting to an rds file; actual plotting will be done in a different script
saveRDS(gwas_plot, file = paste0("figures/grobs/", program, "_", trait, "_manhattan.rds"), compress = FALSE)

