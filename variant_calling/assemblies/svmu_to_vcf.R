# This script converts all svmu output files to VCF using the svmu_to_vcf
# function of the svmutools library

# Loading the svmutools library
library(svmutools)
library(GenomicRanges)

samples <- read.table("assembly_samples.txt", header = FALSE, stringsAsFactors = FALSE)

for(i in 1:nrow(samples)) {
	sample_path <- samples[i, 1]
	sample_name <- samples[i, 2]

	svmu_to_vcf(svmu_file = paste0(sample_name, "/sv.txt"),
		    sample_name = sample_name,
		    ref_fasta = "../../refgenome/Gmax_508_v4.0_mit_chlp.fasta",
		    query_fasta = paste0("../../external_data/genome_assemblies/", sample_path),
		    output_file = paste0(sample_name, ".vcf"),
		    mc.cores = 5, logging = TRUE,
		    remove_duplicates = TRUE, insdup = TRUE, min_insdup = 50,
		    process_inversions = TRUE, min_inv_distance = 5000, inv_maxgap = 10,
		    apply_filters = TRUE, min_size = 50, sv_types = c("INS", "INSDUP", "DEL"),
		    maxcov = 1, max_size = 500000)
}

