# This script creates signals that have not been documented before
# but were found through our analyses as a GRanges object for use in
# downstream analyses

# Loading the GenomicRanges package
library(GenomicRanges)

# Manually entering data for the signal found on chromosome 16 for pubescence color
# This signal goes from 1429763-1436787 for kmers, 1430025-1434520 for Paragraph, and 1380791-1485581 for Platypus
pubescence_color_Gm16_signal <- GRanges(seqnames = "Gm16", ranges = IRanges(start = 1380791, end = 1485581))
pubescence_color_Gm16_signal$signal_id <- "pubescence_color_all_Gm16"
pubescence_color_Gm16_signal$trait <- "pubescence_color"
pubescence_color_Gm16_signal$locus <- "pcGm16"
pubescence_color_Gm16_signal$n_snps <- NA
pubescence_color_Gm16_signal$log_pvalue <- 10 # Choosing a random p-value just so it shows up on Manhattan plots
pubescence_color_Gm16_signal$gene_name_v4 <- NA
pubescence_color_Gm16_signal$common_name <- NA
names(pubescence_color_Gm16_signal) <- pubescence_color_Gm16_signal$signal_id

# Manually entering data for the signal found on chromosome 11 for stem_termination_sn
# This signal goes from 15127026-15402077 for kmers and 15159161-15259163 for Platypus
stem_termination_Gm11_signal <- GRanges(seqnames = "Gm11", ranges = IRanges(start = 15127026, end = 15402077))
stem_termination_Gm11_signal$signal_id <- "stem_termination_sn_Gm11"
stem_termination_Gm11_signal$trait <- "stem_termination"
stem_termination_Gm11_signal$locus <- "stGm11"
stem_termination_Gm11_signal$n_snps <- NA
stem_termination_Gm11_signal$log_pvalue <- 10 # Choosing a random p-value just so it shows up on Manhattan plots
stem_termination_Gm11_signal$gene_name_v4 <- NA
stem_termination_Gm11_signal$common_name <- NA
names(stem_termination_Gm11_signal) <- stem_termination_Gm11_signal$signal_id

# Manually entering data for the signal found on chromosome 15 for pod_color_blbr
# This signal goes from 3908733-4166999 for kmers and 3878655-4084300 for Platypus
pod_color_Gm15_signal <- GRanges(seqnames = "Gm15", ranges = IRanges(start = 3878655, end = 4166999))
pod_color_Gm15_signal$signal_id <- "pod_color_blbr_Gm15"
pod_color_Gm15_signal$trait <- "pod_color"
pod_color_Gm15_signal$locus <- "pdcGm15"
pod_color_Gm15_signal$n_snps <- NA
pod_color_Gm15_signal$log_pvalue <- 10 # Choosing a random p-value just so it shows up on Manhattan plots
pod_color_Gm15_signal$gene_name_v4 <- NA
pod_color_Gm15_signal$common_name <- NA
names(pod_color_Gm15_signal) <- pod_color_Gm15_signal$signal_id

# Manually entering data for the signal found on chromosome 4 for pubescence_form_all
# This signal goes from 39281548-39287059 for kmers and 39232038-39336408 for Platypus
pubescence_form_Gm04_signal <- GRanges(seqnames = "Gm04", ranges = IRanges(start = 39232038, end = 39336408))
pubescence_form_Gm04_signal$signal_id <- "pubescence_form_all_Gm04"
pubescence_form_Gm04_signal$trait <- "pubescence_form"
pubescence_form_Gm04_signal$locus <- "pfGm04"
pubescence_form_Gm04_signal$n_snps <- NA
pubescence_form_Gm04_signal$log_pvalue <- 10 # Choosing a random p-value just so it shows up on Manhattan plots
pubescence_form_Gm04_signal$gene_name_v4 <- NA
pubescence_form_Gm04_signal$common_name <- NA
names(pubescence_form_Gm04_signal) <- pubescence_form_Gm04_signal$signal_id

# Manually entering data for the signal found on chromosome 15 for pubescence_form_all
# This signal goes from 23304848-23307410 for kmers and 23254708-23355585 for Platypus
pubescence_form_Gm15_signal <- GRanges(seqnames = "Gm15", ranges = IRanges(start = 23254708, end = 23355585))
pubescence_form_Gm15_signal$signal_id <- "pubescence_form_all_Gm15"
pubescence_form_Gm15_signal$trait <- "pubescence_form"
pubescence_form_Gm15_signal$locus <- "pfGm15"
pubescence_form_Gm15_signal$n_snps <- NA
pubescence_form_Gm15_signal$log_pvalue <- 10 # Choosing a random p-value just so it shows up on Manhattan plots
pubescence_form_Gm15_signal$gene_name_v4 <- NA
pubescence_form_Gm15_signal$common_name <- NA
names(pubescence_form_Gm15_signal) <- pubescence_form_Gm15_signal$signal_id

# Manually entering data for the signal found on chromosome 09 for seed coat luster ; may correspond to B? locus
# This signal goes from 1084563-1453277 for kmers
seed_coat_luster_Gm09_signal <- GRanges(seqnames = "Gm09", ranges = IRanges(start = 1084563, end = 1453277))
seed_coat_luster_Gm09_signal$signal_id <- "seed_coat_luster_all_Gm09"
seed_coat_luster_Gm09_signal$trait <- "seed_coat_luster"
seed_coat_luster_Gm09_signal$locus <- "sclGm09"
seed_coat_luster_Gm09_signal$n_snps <- NA
seed_coat_luster_Gm09_signal$log_pvalue <- 10 # Choosing a random p-value just so it shows up on Manhattan plots
seed_coat_luster_Gm09_signal$gene_name_v4 <- NA
seed_coat_luster_Gm09_signal$common_name <- NA
names(seed_coat_luster_Gm09_signal) <- seed_coat_luster_Gm09_signal$signal_id

# Manually entering data for the signal found on chromosome 20 for seed coat luster (dull/shiny)
# This signal goes from 32023980-32030600 for kmers and 32024214-32025823
seed_coat_luster_Gm20_signal <- GRanges(seqnames = "Gm20", ranges = IRanges(start = 32023980, end = 32030600))
seed_coat_luster_Gm20_signal$signal_id <- "seed_coat_luster_dullshiny_Gm20"
seed_coat_luster_Gm20_signal$trait <- "seed_coat_luster"
seed_coat_luster_Gm20_signal$locus <- "sclGm20"
seed_coat_luster_Gm20_signal$n_snps <- NA
seed_coat_luster_Gm20_signal$log_pvalue <- 10 # Choosing a random p-value just so it shows up on Manhattan plots
seed_coat_luster_Gm20_signal$gene_name_v4 <- NA
seed_coat_luster_Gm20_signal$common_name <- NA
names(seed_coat_luster_Gm20_signal) <- seed_coat_luster_Gm20_signal$signal_id

# Bundling all those signals together and saving them to file
custom_signals <- c(pubescence_color_Gm16_signal,
		    stem_termination_Gm11_signal,
		    pod_color_Gm15_signal,
		    pubescence_form_Gm04_signal,
		    pubescence_form_Gm15_signal,
		    seed_coat_luster_Gm09_signal,
		    seed_coat_luster_Gm20_signal)

# Saving this signal to file for retrieval in downstream analyses
saveRDS(custom_signals, file = "reference_signals/custom_signals.rds")

