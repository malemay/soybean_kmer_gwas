RSCRIPT := ~/.local/bin/Rscript --quiet
PDFLATEX := ~/.local/texlive/2022/bin/x86_64-linux/pdflatex
BIBTEX := ~/.local/texlive/2022/bin/x86_64-linux/bibtex

# htslib/1.10.2
# ASV_VariantDetector
# add_pvalues
# bamaddrg/1.0
# bayestypertools v. 1.5
# bcftools/1.8
# bcftools/1.10
# bedtools/2.26.0
# bwa/0.7.17
# delta-filter
# edlib-aligner
# filter_kmers
# gemma_0_96
# katcher
# kmers_gwas.py
# LAST
# list_kmers
# manta
# mummer/3.23
# nucmer
# samtools/1.8
# samtools/1.12
# samtools/1.13
# smoove
# SOAPdenovo/2.04
# SvABA
# SVmerge
# svmu
# svmutools
# bbduk
# vcftools/0.1.16
# plink/1.90b5.3
# run_pipeline.pl (TASSEL)
# python/2.7
# htslib/1.8
# htslib/1.10.2
# platypus/0.8.1.1
# python/3.7
# multigrmpy.py

SDIR := additional_files

# Creating a few variables for improving the readability of rules
txdb := refgenome/gmax_v4_genes.rds refgenome/gmax_v4_transcripts.rds refgenome/gmax_v4_exons.rds refgenome/gmax_v4_cds.rds
gmax_db := refgenome/Gmax4_db.bck refgenome/Gmax4_db.des refgenome/Gmax4_db.prj refgenome/Gmax4_db.sds refgenome/Gmax4_db.ssp refgenome/Gmax4_db.suf refgenome/Gmax4_db.tis
refgen := refgenome/Gmax_508_v4.0_mit_chlp.fasta
signals_gr := utilities/all_signals.rds
grobdir := figures/grobs
addfile := additional_files/additional_file_1.tex

maintables := tables/loci_table.csv

suptables := tables/FLOWER.COLOR_gwas_table.csv \
	tables/PUBESCENCE.COLOR_gwas_table.csv \
	tables/SEED.COAT.COLOR_gwas_table.csv \
	tables/STEM.TERMINATION.TYPE_gwas_table.csv \
	tables/HILUM.COLOR_gwas_table.csv \
	tables/POD.COLOR_gwas_table.csv \
	tables/PUBESCENCE.FORM_gwas_table.csv \
	tables/PUBESCENCE.DENSITY_gwas_table.csv \
	tables/SEED.COAT.LUSTER_gwas_table.csv \
	tables/MATURITY.GROUP_gwas_table.csv \
	tables/signals_table.csv

mainfigures := figures/flower_color_W1_main_figure.png \
	figures/pubescence_color_nogray_Td_main_figure.png \
	figures/seed_coat_color_greenyellow_G_main_figure.png \
	figures/pubescence_density_Ps_main_figure.png \
	figures/pubescence_form_all_Pa1_main_figure.png \
	figures/stem_termination_sn_main_figure.png

# Getting the LD png figures to generate dynamically from additional_files/additional_file_1.tex
manhattanplots := $(shell grep '^[%]*\\manhattanplot' $(addfile) | sed -E 's/[%]*\\manhattanplot\{([a-zA-Z_]*)\}.*$$/\1/' | xargs -I {} echo figures/{}_manhattan.png)
ldplots := $(shell grep '^[%]*\\ldplot' $(addfile) | sed -E 's/[%]*\\ldplot\{([a-zA-Z_]*)\}.*$$/\1/' | xargs -I {} echo figures/{}_ld.png)
geneplots := $(shell grep '^[%]*\\geneplot' $(addfile) | sed -E 's/[%]*\\geneplot\{([a-zA-Z0-9_]*)\}.*$$/\1/' | xargs -I {} echo figures/{}_gene.png)
kmerplots := $(shell grep '^[%]*\\kmerplot' $(addfile) | sed -E 's/[%]*\\kmerplot\{([a-zA-Z0-9_]*)\}.*$$/\1/' | xargs -I {} echo figures/{}_kmers.png)
signalplots := $(shell grep '^[%]*\\signalplot' $(addfile) | sed -E 's/[%]*\\signalplot\{([a-zA-Z0-9_]*)\}.*$$/\1/' | xargs -I {} echo figures/{}_signal.png)

# Grouping all the supplementary figures together
supfigures := $(manhattanplots) \
	$(signalplots) \
	$(geneplots) \
	$(kmerplots) \
	$(ldplots) \
	figures/concordance_histogram.png

# Creating variables for objects that are used in the interpretation of the results
# We only create these tables for loci that have a corresponding signal plot
topgranges := $(foreach prog,platypus paragraph kmers,$(shell grep '^[%]*\\signalplot' $(addfile) | sed -E 's/[%]*\\signalplot\{([a-zA-Z0-9_]*)\}.*$$/\1/' | xargs -I {} echo gwas_results/$(prog)/{}_top_markers.rds))

allgenes := $(foreach prog,platypus paragraph kmers,$(shell grep '^[%]*\\signalplot' $(addfile) | sed -E 's/[%]*\\signalplot\{([a-zA-Z0-9_]*)\}.*$$/\1/' | xargs -I {} echo gwas_results/$(prog)/{}_all_genes.tsv))

topgenes := $(foreach prog,platypus paragraph kmers,$(shell grep '^[%]*\\signalplot' $(addfile) | sed -E 's/[%]*\\signalplot\{([a-zA-Z0-9_]*)\}.*$$/\1/' | xargs -I {} echo gwas_results/$(prog)/{}_top_genes.tsv))

nearestgene := $(foreach prog,platypus paragraph kmers,$(shell grep '^[%]*\\signalplot' $(addfile) | sed -E 's/[%]*\\signalplot\{([a-zA-Z0-9_]*)\}.*$$/\1/' | xargs -I {} echo gwas_results/$(prog)/{}_nearest_gene.tsv))

genetables := $(allgenes) $(topgenes) $(nearestgene)

phenodata := $(shell cut -f1 utilities/kmer_plot_ranges.txt | xargs -I {} echo gwas_results/kmers/{}_phenodata.rds)


# The first target of the Makefile
# It generates the manuscript as well as the supplementary files, and data for interpretation
all: $(SDIR)/manuscript.pdf \
	$(genetables) $(phenodata) \
	$(SDIR)/supplemental_file_2.csv \
	$(SDIR)/supplemental_file_3.csv \
	$(SDIR)/supplemental_file_4.csv

# Making sure that intermediate files do not get deleted
.PRECIOUS: $(grobdir)/platypus_%_manhattan.rds \
	$(grobdir)/paragraph_%_manhattan.rds \
	$(grobdir)/kmers_%_manhattan.rds \
	$(grobdir)/platypus_%_signal.rds \
	$(grobdir)/paragraph_%_signal.rds \
	$(grobdir)/kmers_%_signal.rds \
	$(grobdir)/platypus_%_gene.rds \
	$(grobdir)/paragraph_%_gene.rds \
	$(grobdir)/kmers_%_gene.rds \
	gwas_results/kmers/%_threshold_5per.txt \
	gwas_results/kmer_consensus/%_sequences.fa \
	gwas_results/kmers/%_clustered_ld.txt \
	cnv_analysis/hps_cnv_range.rds \
	cnv_analysis/ps_cnv_range.rds

# Adding some more intermediate files to .PRECIOUS
$(foreach signal,$(shell cut -d "," -f1 utilities/signal_ids.txt | xargs -I {} echo gwas_results/%/{}_gwas_locus.rds),$(eval .PRECIOUS: $(signal)))
$(foreach signal,$(shell cut -d "," -f1 utilities/signal_ids.txt | xargs -I {} echo gwas_results/%/{}_signal_locus.rds),$(eval .PRECIOUS: $(signal)))
$(foreach prog,platypus paragraph kmers,$(eval .PRECIOUS : gwas_results/$(prog)/%_signal.rds))
$(foreach prog,platypus paragraph kmers,$(eval .PRECIOUS : gwas_results/$(prog)/%_gwas.rds))

.PRECIOUS : $(topgranges)

# Compiling the manuscript and additional file 1 (in 1 pdf)
$(SDIR)/manuscript.pdf: $(SDIR)/manuscript.tex \
	$(SDIR)/main_text.tex \
	$(SDIR)/additional_file_1.tex \
	$(supfigures) \
	$(suptables) \
	$(mainfigures) \
	$(maintables) \
	additional_files/references.bib \
	$(SDIR)/variables.txt
	cd $(SDIR) ; $(PDFLATEX) manuscript.tex ; $(BIBTEX) main_text ; $(BIBTEX) additional_file_1 ; $(PDFLATEX) manuscript.tex ; $(PDFLATEX) manuscript.tex

# Creating the list of variables stored in additional_files/variables.txt, for retrival in additional file 1
$(SDIR)/variables.txt: $(SDIR)/make_variables.R \
	$(ldplots) \
	phenotypic_data/phenotypic_data.csv \
	utilities/trait_names.txt \
	$(signals_gr)
	$(RSCRIPT) $<

# Creating the supplemental CSV files
$(SDIR)/%.csv: $(SDIR)/%.R
	$(RSCRIPT) $<

# Additional dependencies for supplemental files
$(SDIR)/supplemental_file_2.csv: utilities/correct_sra_metadata.csv utilities/srr_id_correspondence.txt
# DEPENDENCIES on samtools stats and manifest files are missing for this rule
$(SDIR)/supplemental_file_3.csv: utilities/srr_id_correspondence.txt illumina_data/soysnp50k_genotyping/gtcheck_results.tsv
$(SDIR)/supplemental_file_4.csv: phenotypic_data/phenotypic_data.csv

# MAIN FIGURES --------------------------------------------------
figures/flower_color_W1_main_figure.png: figures/flower_color_main_figure.R \
	utilities/kmer_plot_ranges.txt \
	gwas_results/kmer_consensus/flower_color_W1_plotting_data.rds \
	gwas_results/kmer_consensus/flower_color_W1_difflist.rds \
	gwas_results/kmer_consensus/flower_color_W1_causal_gene.rds \
	gwas_results/kmers/flower_color_W1_phenodata.rds \
	$(txdb) \
	figures/main_figure_functions.R \
	figures/grobs/paragraph_flower_color_manhattan.rds \
	figures/grobs/kmers_flower_color_manhattan.rds \
	figures/grobs/paragraph_flower_color_W1_signal.rds \
	figures/grobs/kmers_flower_color_W1_signal.rds
	$(RSCRIPT) $< flower_color_W1

figures/pubescence_color_nogray_Td_main_figure.png: figures/pubescence_color_nogray_main_figure.R \
	utilities/kmer_plot_ranges.txt \
	gwas_results/kmer_consensus/pubescence_color_nogray_Td_plotting_data.rds \
	gwas_results/kmer_consensus/pubescence_color_nogray_Td_difflist.rds \
	gwas_results/kmer_consensus/pubescence_color_nogray_Td_causal_gene.rds \
	gwas_results/kmers/pubescence_color_nogray_Td_phenodata.rds \
	$(txdb) \
	figures/main_figure_functions.R \
	figures/grobs/platypus_pubescence_color_nogray_manhattan.rds \
	figures/grobs/kmers_pubescence_color_nogray_manhattan.rds \
	figures/grobs/platypus_pubescence_color_nogray_Td_signal.rds \
	figures/grobs/kmers_pubescence_color_nogray_Td_signal.rds
	$(RSCRIPT) $< pubescence_color_nogray_Td

figures/seed_coat_color_greenyellow_G_main_figure.png: figures/seed_coat_color_greenyellow_main_figure.R \
	utilities/kmer_plot_ranges.txt \
	gwas_results/kmer_consensus/seed_coat_color_greenyellow_G_plotting_data.rds \
	gwas_results/kmer_consensus/seed_coat_color_greenyellow_G_difflist.rds \
	gwas_results/kmer_consensus/seed_coat_color_greenyellow_G_causal_gene.rds \
	gwas_results/kmers/seed_coat_color_greenyellow_G_phenodata.rds	 \
	$(txdb) \
	figures/main_figure_functions.R \
	figures/grobs/platypus_seed_coat_color_greenyellow_manhattan.rds \
	figures/grobs/kmers_seed_coat_color_greenyellow_manhattan.rds \
	figures/grobs/platypus_seed_coat_color_greenyellow_G_signal.rds \
	figures/grobs/kmers_seed_coat_color_greenyellow_G_signal.rds
	$(RSCRIPT) $< seed_coat_color_greenyellow_G

# Main figure for pubescence density
# The mapped reads processed through katcher have yet to be included in this rule
figures/pubescence_density_Ps_main_figure.png: figures/pubescence_density_main_figure.R \
	figures/grobs/kmers_pubescence_density_manhattan.rds \
	figures/grobs/kmers_pubescence_density_Ps_gene.rds \
	figures/main_figure_functions.R \
	phenotypic_data/phenotypic_data.csv \
	$(shell cut -d ' ' -f1 utilities/srr_id_correspondence.txt | xargs -I {} echo illumina_data/merged_bams/{}_merged.bam) \
	$(shell cut -d ' ' -f1 utilities/srr_id_correspondence.txt | xargs -I {} echo illumina_data/merged_bams/{}_merged.bam.bai) \
	$(shell cut -d ' ' -f1 utilities/srr_id_correspondence.txt | xargs -I {} echo sv_genotyping/paragraph/manifest_files/{}_manifest.txt)
	$(RSCRIPT) $<

figures/pubescence_form_all_Pa1_main_figure.png: figures/pubescence_form_all_main_figure.R \
	utilities/kmer_plot_ranges.txt \
	$(txdb) \
	figures/main_figure_functions.R \
	figures/grobs/platypus_pubescence_form_all_manhattan.rds \
	figures/grobs/kmers_pubescence_form_all_manhattan.rds \
	figures/grobs/platypus_pubescence_form_all_Pa1_signal.rds \
	figures/grobs/kmers_pubescence_form_all_Pa1_signal.rds \
	figures/grobs/platypus_pubescence_form_all_Pa1_gene.rds \
	figures/grobs/kmers_pubescence_form_all_Pa1_gene.rds
	$(RSCRIPT) $< pubescence_form_all_Pa1

figures/stem_termination_sn_main_figure.png: figures/stem_termination_sn_main_figure.R \
	figures/main_figure_functions.R \
	figures/grobs/kmers_stem_termination_sn_manhattan.rds \
	gwas_results/kmers/stem_termination_sn_clustered_ld.txt \
	gwas_results/kmers/stem_termination_sn_gwas.rds
	$(RSCRIPT) $<

# SIGNALS --------------------------------------------------

# This script prepares the reference signals GRanges object
$(signals_gr) utilities/signal_ids.txt utilities/gene_signal_ids.txt: utilities/make_signals_granges.R \
	reference_signals/bandillo2017_signals.rds \
	reference_signals/bandillo2015_signals.rds \
	reference_signals/deronne2022_signal.rds \
	reference_signals/custom_signals.rds \
	utilities/trait_names.txt
	$(RSCRIPT) utilities/make_signals_granges.R

# Generating the reference signals obtained from Bandillo et al. (2017) for qualitative traits
reference_signals/bandillo2017_signals.rds: reference_signals/format_bandillo2017_signals.R \
	refgenome/glyma.Wm82.gnm1.FCtY.genome_main.fna \
	refgenome/Gmax_508_v4.0_mit_chlp.fasta \
	reference_signals/bandillo2017_signals_curated.tsv \
	refgenome/lookup_v1_to_v4.rds \
	refgenome/soybase_genome_annotation_v4.0_04-20-2021.txt \
	refgenome/gmax_v4_genes.rds
	$(RSCRIPT) reference_signals/format_bandillo2017_signals.R

# Generating the reference signals obtained from Bandillo et al. (2015) for oil and protein
reference_signals/bandillo2015_signals.rds: reference_signals/format_bandillo2015_signals.R \
	refgenome/glyma.Wm82.gnm1.FCtY.genome_main.fna \
	refgenome/Gmax_508_v4.0_mit_chlp.fasta \
	reference_signals/bandillo2015_table1_curated.csv
	$(RSCRIPT) reference_signals/format_bandillo2015_signals.R

# Generating the reference signal obtained from de Ronne et al. (2022) for resistance to P. sojae
reference_signals/deronne2022_signal.rds: reference_signals/format_deronne2022_signal.R
	$(RSCRIPT) reference_signals/format_deronne2022_signal.R

# Generating the customs signals observed from our analyses for various signals
reference_signals/custom_signals.rds: reference_signals/format_custom_signals.R
	$(RSCRIPT) reference_signals/format_custom_signals.R

# Putting the k-mer 5% thresholds on the same scale as the thresholds for the other programs
gwas_results/kmers/%_threshold_5per.txt: gwas_results/scale_kmer_thresholds.R gwas_results/kmers/KMER_GWAS
	$(RSCRIPT) gwas_results/scale_kmer_thresholds.R $*

# Generating the .csv file with data regarding which reference signals were found and their p-values
tables/signals_table.csv: tables/signals_table.R \
	$(signals_gr) \
	$(shell cat utilities/trait_names.txt | xargs -I {} echo gwas_results/platypus/{}_signal.rds) \
	$(shell cat utilities/trait_names.txt | xargs -I {} echo gwas_results/paragraph/{}_signal.rds) \
	$(shell cat utilities/trait_names.txt | xargs -I {} echo gwas_results/kmers/{}_signal.rds)
	$(RSCRIPT) tables/signals_table.R

# Generating the .csv file with data regarding the performance of various approaches on loci for which genes are known
tables/loci_table.csv: tables/loci_table.R \
	utilities/all_signals.rds \
	refgenome/gmax_v4_genes.rds \
	cnv_analysis/i_cnv_range.rds \
	cnv_analysis/hps_cnv_range.rds \
	cnv_analysis/ps_cnv_range.rds \
	$(signalplots)
	$(RSCRIPT) $<

# PHENOTYPIC DATA --------------------------------------------------

# Creating the table with the numeric coding for a given trait
tables/%_gwas_table.csv: tables/gwas_table.R \
	phenotypic_data/trait_names.rds \
	phenotypic_data/pheno_names_lookup.rds \
	phenotypic_data/phenotypic_data.csv \
	phenotypic_data/lookup_tables.rds
	$(RSCRIPT) tables/gwas_table.R $*

# Creating the lookup tables used to translate the phenotypes into numeric codes
phenotypic_data/lookup_tables.rds: phenotypic_data/lookup_tables.R
	$(RSCRIPT) phenotypic_data/lookup_tables.R

# Creating the lookup table used to translate the GRIN trait names into the ones used in this analysis
phenotypic_data/trait_names.rds: phenotypic_data/trait_names.R
	$(RSCRIPT) phenotypic_data/trait_names.R

# KMER HAPLOTYPE PLOTS --------------------------------------------------
# Generating the k-mer plot from the consensus sequences (the k-mer p-values are missing from the list of dependencies)
figures/%_kmers.png gwas_results/kmers/%_phenodata.rds gwas_results/kmer_consensus/%_plotting_data.rds gwas_results/kmer_consensus/%_difflist.rds gwas_results/kmer_consensus/%_causal_gene.rds: figures/kmer_plot.R \
	$(signals_gr) \
	$(txdb) \
	phenotypic_data/pheno_names_lookup.rds \
	utilities/kmer_plot_ranges.txt \
	phenotypic_data/trait_names.rds \
	phenotypic_data/phenotypic_data.csv \
	gwas_results/kmer_consensus/%_sequences.fa
	$(RSCRIPT) figures/kmer_plot.R $*

# Extracting the consensus sequences from significant k-mer assemblies for a given locus
# The BWA alignment results are missing from the list of dependencies
gwas_results/kmer_consensus/%_sequences.fa: gwas_results/gather_consensus.sh \
	utilities/kmer_plot_ranges.txt \
	refgenome/Gmax_508_v4.0_mit_chlp.fasta
	gwas_results/gather_consensus.sh $*

# Creating lookup tables to simplify phenotype names
phenotypic_data/pheno_names_lookup.rds: phenotypic_data/pheno_names_lookup.R
	$(RSCRIPT) $<

# KMER LD PLOTS --------------------------------------------------

# Generating the LD plot from the clustered matrix
figures/%_ld.png: figures/ld_plot.R \
	gwas_results/kmers/%_clustered_ld.txt \
	gwas_results/kmers/%_gwas.rds
	$(RSCRIPT) $< $*

# Generating the clustered LD matrix from the GWAS results and the full k-mers table
gwas_results/kmers/%_clustered_ld.txt: gwas_results/ld_analysis.R \
	kmers_table/kmers_table.table \
	kmers_table/kmers_table.names \
	gwas_results/kmers/%_gwas.rds
	$(RSCRIPT) $< $*

# GWAS RESULTS --------------------------------------------------

# Preparing the GWAS results for platypus paragraph
$(foreach prog,platypus paragraph,$(eval gwas_results/$(prog)/%_gwas.rds: gwas_results/format_gwas_results.R \
	$(refgen) \
	filtered_variants/$(prog)/filtered_variants.vcf.gz \
	gwas_results/$(prog)/%_gwas.csv \
	gwas_results/$(prog)/%_threshold_5per.txt ; \
	$(RSCRIPT) gwas_results/format_gwas_results.R $$* $(prog)))

# kmers need special treatment
gwas_results/kmers/%_gwas.rds: gwas_results/format_gwas_results.R \
	$(refgen) \
	gwas_results/kmers/%_kmer_positions.rds \
	gwas_results/kmers/%_threshold_5per.txt
	$(RSCRIPT) gwas_results/format_gwas_results.R $* kmers

# SIGNALS --------------------------------------------------
# Creating a GRanges object containing the signal(s) for all traits for each program but platypus
# Also creating subsets of GWAS results that only overlap signals
$(foreach prog,paragraph kmers,$(eval gwas_results/$(prog)/%_signal.rds gwas_results/$(prog)/%_gwas_subset.rds : gwas_results/find_signals.R \
	gwas_results/$(prog)/%_gwas.rds \
	gwas_results/$(prog)/%_threshold_5per.txt ; \
	$(RSCRIPT) gwas_results/find_signals.R $$* $(prog)))

# Platypus requires special treatment because of the pruned markers
gwas_results/platypus/%_signal.rds gwas_results/platypus/%_gwas_subset.rds : gwas_results/find_signals.R \
	gwas_results/platypus/%_gwas.rds \
	gwas_results/platypus/%_threshold_5per.txt \
	filtered_variants/platypus_full.vcf.gz \
	filtered_variants/platypus_gapit_kinship.rds \
	filtered_variants/platypus_gapit_pca.rds \
	phenotypic_data/phenotypic_data.csv
	$(RSCRIPT) gwas_results/find_signals.R $* platypus

# GENE ANALYSIS  --------------------------------------------------
# Generating a GRanges object and tsv files for each locus
$(foreach prog,platypus paragraph kmers,$(eval \
	gwas_results/$(prog)/%_top_markers.rds gwas_results/$(prog)/%_all_genes.tsv \
	gwas_results/$(prog)/%_top_genes.tsv gwas_results/$(prog)/%_nearest_gene.tsv : \
	gwas_results/gene_analysis.R \
	$(signals_gr) \
	refgenome/gmax_v4_genes.rds \
	refgenome/soybase_gmax_v4_annotations.rds \
	gwas_results/$(prog)/%_gwas_locus.rds \
	gwas_results/$(prog)/%_signal_locus.rds ; \
	$(RSCRIPT) gwas_results/gene_analysis.R $$* $(prog)))

# MANHATTAN PLOTS --------------------------------------------------

# This block assembles the three manhattan subplots for a given trait
figures/%_manhattan.png: figures/manhattan_plot.R \
	$(grobdir)/platypus_%_manhattan.rds \
	$(grobdir)/paragraph_%_manhattan.rds \
	$(grobdir)/kmers_%_manhattan.rds
	$(RSCRIPT) figures/manhattan_plot.R $*

# Preparing the manhattan subplots for platypus, paragraph and k-mers
$(foreach prog,platypus paragraph kmers,$(eval $(grobdir)/$(prog)_%_manhattan.rds: figures/manhattan_subplot.R \
	$(signals_gr) \
	gwas_results/$(prog)/%_gwas.rds \
	gwas_results/$(prog)/%_signal.rds \
	gwas_results/$(prog)/%_threshold_5per.txt ; \
	$(RSCRIPT) figures/manhattan_subplot.R $$* $(prog)))

# SIGNAL PLOTS --------------------------------------------------

# This block assembles the three signal subplots for a given locus
figures/%_signal.png: figures/signal_plot.R \
	$(signals_gr) \
	$(txdb) \
	$(grobdir)/platypus_%_signal.rds \
	$(grobdir)/paragraph_%_signal.rds \
	$(grobdir)/kmers_%_signal.rds
	$(RSCRIPT) figures/signal_plot.R $*

# Preparing the signal subplots for each of platypus, k-mers and paragraph
$(foreach prog,platypus paragraph kmers,$(eval $(grobdir)/$(prog)_%_signal.rds: figures/signal_subplot.R \
	$(signals_gr) \
	utilities/cnv_genes.txt \
	gwas_results/$(prog)/%_top_markers.rds \
	refgenome/gmax_v4_genes.rds \
	gwas_results/$(prog)/%_gwas_locus.rds \
	gwas_results/$(prog)/%_signal_locus.rds ; \
	$(RSCRIPT) figures/signal_subplot.R $$* $(prog)))

# GENE PLOTS --------------------------------------------------

# This block assembles the three gene subplots for a given locus
figures/%_gene.png: figures/gene_plot.R \
	$(signals_gr) \
	$(txdb) \
	$(grobdir)/platypus_%_gene.rds \
	$(grobdir)/paragraph_%_gene.rds \
	$(grobdir)/kmers_%_gene.rds
	$(RSCRIPT) figures/gene_plot.R $*

# Preparing the gene subplots for each of platypus, k-mers and paragraph
$(foreach prog,platypus paragraph kmers,$(eval $(grobdir)/$(prog)_%_gene.rds: figures/gene_subplot.R \
	$(signals_gr) \
	utilities/cnv_genes.txt \
	utilities/variant_ranges.txt \
	$(shell cut -f2 utilities/cnv_genes.txt | uniq) \
	refgenome/gmax_v4_genes.rds \
	gwas_results/$(prog)/%_gwas_locus.rds \
	gwas_results/$(prog)/%_signal_locus.rds ; \
	$(RSCRIPT) figures/gene_subplot.R $$* $(prog)))

# CNV ANALYSES --------------------------------------------------
# Preparing the CNV ranges for the loci that harbor CNV (only Ps and Hps so far)
cnv_analysis/%_cnv_range.rds: cnv_analysis/%_analysis.R \
	cnv_analysis/cnv_functions.R \
	utilities/srr_id_correspondence.txt \
	$(shell cut -d ' ' -f1 utilities/srr_id_correspondence.txt | xargs -I {} echo illumina_data/merged_bams/{}_merged.bam) \
	$(shell cut -d ' ' -f1 utilities/srr_id_correspondence.txt | xargs -I {} echo illumina_data/merged_bams/{}_merged.bam.bai) \
	$(shell cut -d ' ' -f1 utilities/srr_id_correspondence.txt | xargs -I {} echo sv_genotyping/paragraph/manifest_files/{}_manifest.txt) \
	refgenome/gmax_v4_genes.rds
	$(RSCRIPT) $<

# The I locus is a special case because it does not have as many prerequisites
cnv_analysis/i_cnv_range.rds: cnv_analysis/i_analysis.R \
	cnv_analysis/i_locus.bam
	$(RSCRIPT) $<

# The script that creates the alignment of the I locus on the reference assembly
cnv_analysis/i_locus.bam: cnv_analysis/extract_I_locus.sh \
	external_data/BAC77G7-a.fasta \
	external_data/BAC77G7-a.fasta.fai \
	$(refgen)
	$<

# SYMLINKS --------------------------------------------------

# Creating locus-specific symlinks to signal files of the relevant phenotype
$(foreach signal,$(shell cut -d "," -f1 utilities/signal_ids.txt | xargs -I {} echo gwas_results/%/{}_signal_locus.rds),$(eval $(signal) : \
	$(shell echo $(signal) | sed -E 's/[^_]+_signal_locus/signal/') \
	gwas_results/create_symlink.R ; \
	$(RSCRIPT) gwas_results/create_symlink.R $$< $$@))

# Creating locus-specific symlinks to GWAS results files of the relevant phenotype
$(foreach signal,$(shell cut -d "," -f1 utilities/signal_ids.txt | xargs -I {} echo gwas_results/%/{}_gwas_locus.rds),$(eval $(signal) : \
	$(shell echo $(signal) | sed -E 's/[^_]+_gwas_locus/gwas_subset/') \
	gwas_results/create_symlink.R ; \
	$(RSCRIPT) gwas_results/create_symlink.R $$< $$@))


# SUPPLEMENTAL FIGURE - HISTOGRAM of concordance results
figures/concordance_histogram.png: figures/concordance_histogram.R \
	illumina_data/soysnp50k_genotyping/gtcheck_results.tsv
	$(RSCRIPT) $<

# COMPARING THE GENOTYPES DERIVED FROM WGS TO SOYSNP50K
illumina_data/soysnp50k_genotyping/gtcheck_results.tsv: illumina_data/soysnp50k_genotyping/gtcheck_analysis.R \
	illumina_data/soysnp50k_genotyping/pi_ids.txt \
	illumina_data/soysnp50k_genotyping/GTCHECK
	$(RSCRIPT) $<

illumina_data/soysnp50k_genotyping/pi_ids.txt illumina_data/soysnp50k_genotyping/read_groups.txt illumina_data/soysnp50k_genotyping/alleles.tsv.gz: \
	illumina_data/soysnp50k_genotyping/prepare_files.sh \
	phenotypic_data/phenotypic_data.csv \
	illumina_data/soysnp50k_genotyping/soysnp50K_gmax_v4.vcf
	$<

illumina_data/soysnp50k_genotyping/soysnp50K_gmax_v4.vcf: illumina_data/soysnp50k_genotyping/convert_positions.R \
	refgenome/Gmax_nuclv2_mit_chlp.fasta \
	refgenome/Gmax_508_v4.0_mit_chlp.fasta \
	illumina_data/soysnp50k_genotyping/soysnp50K_maf10.vcf
	$(RSCRIPT) $<

illumina_data/soysnp50k_genotyping/soysnp50K_maf10.vcf: illumina_data/soysnp50k_genotyping/filter_soysnp50k.sh \
	phenotypic_data/phenotypic_data.csv \
	external_data/soysnp50k_wm82.a2_41317.vcf.gz 
	$<

illumina_data/soysnp50k_genotyping/GTCHECK: illumina_data/soysnp50k_genotyping/mpileup_all.sh \
	refgenome/Gmax_508_v4.0_mit_chlp.fasta \
	illumina_data/merged_bams/ILLUMINA_BAM_MERGING \
	illumina_data/soysnp50k_genotyping/alleles.tsv.gz \
	illumina_data/soysnp50k_genotyping/read_groups.txt \
	illumina_data/soysnp50k_genotyping/pi_ids.txt \
	illumina_data/soysnp50k_genotyping/soysnp50K_gmax_v4.vcf
	$<

#
# ILLUMINA DATA PROCESSING --------------------------------------------------
#
illumina_data/merged_bams/ILLUMINA_BAM_MERGING: illumina_data/merge_bams.sh \
	illumina_data/index_bams.sh \
	utilities/srr_id_correspondence.txt \
	illumina_data/BWA_ILLUMINA_MAPPING
	illumina_data/merge_bams.sh ; illumina_data/index_bams.sh

illumina_data/BWA_ILLUMINA_MAPPING: illumina_data/bwa_mapping.sh \
	refgenome/Gmax_508_v4.0_mit_chlp.fasta \
	utilities/correct_sra_metadata.csv \
	illumina_data/BBDUK_TRIMMING
	$<

illumina_data/BBDUK_TRIMMING: illumina_data/bbduk_trimming.sh \
	utilities/correct_sra_metadata.csv \
	external_data/adapters.fa \
	$(shell cut -d, -f1 utilities/correct_sra_metadata.csv | tail -n +2 | sed 's/"//g' | xargs -I {} echo "illumina_data/raw_fastq/{}_1.fastq.gz illumina_data/raw_fastq/{}_2.fastq.gz")
	$<

# COMPUTING GWAS ANALYSES --------------------------------------------------

# Computing the GWAS analyses for all traits for Platypus
gwas_results/platypus/%_gwas.csv: gwas_results/platypus/platypus_gwas.R \
	phenotypic_data/phenotypic_data.csv \
	variant_calling/platypus/platypus_formatted.hmp.txt
	$(RSCRIPT) $< $*

# Getting the thresholds for all traits for Platypus
gwas_results/platypus/%_threshold_5per.txt: gwas_results/platypus/platypus_thresholds.R \
	phenotypic_data/phenotypic_data.csv \
	variant_calling/platypus/platypus_formatted.hmp.txt
	$(RSCRIPT) $< $*

# Computing the GWAS analyses for all traits for Paragraph
gwas_results/paragraph/%_gwas.csv: gwas_results/paragraph/paragraph_gwas.R \
	phenotypic_data/phenotypic_data.csv \
	sv_genotyping/paragraph/paragraph_formatted.hmp.txt
	$(RSCRIPT) $< $*

# Getting the thresholds for all traits for Paragraph
gwas_results/paragraph/%_threshold_5per.txt: gwas_results/paragraph/paragraph_thresholds.R \
	phenotypic_data/phenotypic_data.csv \
	sv_genotyping/paragraph/paragraph_formatted.hmp.txt
	$(RSCRIPT) $< $*

# SNP AND INDEL GENOTYPING AND FILTERING FOR GWAS --------------------------------------------------

variant_calling/platypus/platypus_formatted.hmp.txt: variant_calling/platypus/platypus_filtering.sh \
	refgenome/Gmax_508_v4.0_mit_chlp.fasta.fai \
	variant_calling/platypus/platypus_all.vcf
	$<

variant_calling/platypus/platypus_all.vcf: variant_calling/platypus/platypus_call.sh \
	refgenome/Gmax_508_v4.0_mit_chlp.fasta \
	illumina_data/merged_bams/ILLUMINA_BAM_MERGING
	$<

# SV GENOTYPING USING PARAGRAPH AND FILTERING

# Filtering the Paragraph genotype calls and preparing them for input to GAPIT
sv_genotyping/paragraph/paragraph_formatted.hmp.txt: sv_genotyping/paragraph/paragraph_filtering.sh \
	utilities/srr_id_correspondence.txt \
	sv_genotyping/paragraph/PARAGRAPH_GENOTYPING
	$<

sv_genotyping/paragraph/PARAGRAPH_GENOTYPING: sv_genotyping/paragraph/paragraph_genotyping.sh \
	refgenome/Gmax_508_v4.0_mit_chlp.fasta \
	sv_genotyping/paragraph/MANIFEST_FILES \
	illumina_data/merged_bams/ILLUMINA_BAM_MERGING \
	sv_genotyping/paragraph/all_svs_padded.vcf
	$<

sv_genotyping/paragraph/all_svs_padded.vcf: sv_genotyping/paragraph/prepare_variants.sh \
	variant_calling/merging/candidate_svs.vcf \
	refgenome/Gmax_508_v4.0_mit_chlp.fasta \
	sv_genotyping/paragraph/truncate_variants.awk \
	sv_genotyping/paragraph/addMissingPaddingGmax4.py
	$<

sv_genotyping/paragraph/MANIFEST_FILES: sv_genotyping/paragraph/manifest_files.R \
	utilities/srr_id_correspondence.txt \
	illumina_data/SAMTOOLS_COVERAGE \
	illumina_data/SAMTOOLS_STATS
	$(RSCRIPT) $<

# Computing some coverage statistics to create the manifest files
illumina_data/SAMTOOLS_COVERAGE: illumina_data/samtools_coverage.sh \
	utilities/srr_id_correspondence.txt \
	illumina_data/merged_bams/ILLUMINA_BAM_MERGING
	$<

# Computing some statistics to create the manifest files
illumina_data/SAMTOOLS_STATS: illumina_data/samtools_stats.sh \
	utilities/srr_id_correspondence.txt \
	illumina_data/merged_bams/ILLUMINA_BAM_MERGING
	$<

# MERGING OF SVS CALLED USING VARIOUS TOOLS --------------------------------------------------

variant_calling/merging/candidate_svs.vcf: variant_calling/merging/merge_assembly_svs.sh \
	refgenome/Gmax_508_v4.0_mit_chlp.fasta \
	variant_calling/merging/assembly_merging_files.txt \
	variant_calling/merging/usda_svs.vcf \
	variant_calling/assemblies/assembly_svs.vcf
	$<

variant_calling/merging/usda_svs.vcf: variant_calling/merging/merge_nanopore_svs.sh \
	refgenome/Gmax_508_v4.0_mit_chlp.fasta \
	variant_calling/merging/illumina.clustered.vcf \
	external_data/nanopore_svs/nanopore_svs.vcf \
	variant_calling/merging/select_svs.R \
	variant_calling/merging/nanopore_merging_files.txt
	$<

external_data/nanopore_svs/nanopore_svs.vcf: external_data/nanopore_svs/merge_svs.sh \
	refgenome/Gmax_508_v4.0_mit_chlp.fasta \
	external_data/nanopore_svs/files.txt \
	external_data/nanopore_svs/merge_realigned.R \
	$(shell cat external_data/nanopore_svs/files.txt | xargs -I {} echo external_data/nanopore_svs/{}) \
	external_data/nanopore_svs/header_lines.txt
	$<

variant_calling/merging/illumina.clustered.vcf: variant_calling/merging/merge_illumina_svs.sh \
	refgenome/Gmax_508_v4.0_mit_chlp.fasta \
	variant_calling/asmvar/asmvar_svmerged.clustered.vcf \
	variant_calling/manta/manta_svmerged.clustered.vcf \
	variant_calling/smoove/smoove_svmerged.clustered.vcf \
	variant_calling/svaba/svaba_svmerged.clustered.vcf
	$<

# CALLING SVS WITH ASMVAR --------------------------------------------------

# Merging the variants called on all samples with AsmVar
variant_calling/asmvar/asmvar_svmerged.clustered.vcf: variant_calling/asmvar/asmvar_svmerge.sh \
	refgenome/Gmax_508_v4.0_mit_chlp.fasta \
	variant_calling/asmvar/asmvar_filtered.vcf
	$<

# Filtering the variants called from each sample with AsmVar and performing a first merge
variant_calling/asmvar/asmvar_filtered.vcf: variant_calling/asmvar/asmvar_filter.sh \
	refgenome/Gmax_508_v4.0_mit_chlp.fasta \
	refgenome/Gmax_508_v4.0_mit_chlp.fasta.fai \
	scripts/add_svtype.awk \
	scripts/extract_svs_50.awk \
	variant_calling/asmvar/svtype_header_line.txt \
	variant_calling/asmvar/ASMVAR_CALLING
	$<

# Assembling genomes with SOAPdenovo2, aligning them with LAST and calling variants with AsmVar
variant_calling/asmvar/ASMVAR_CALLING: variant_calling/asmvar/asmvar_call.sh \
	utilities/srr_id_correspondence.txt \
	refgenome/Gmax_508_v4.0_mit_chlp.fasta \
	$(gmax_db) \
	illumina_data/BBDUK_TRIMMING
	$<

# CALLING SVS WITH SMOOVE --------------------------------------------------

# Merging the variants called with smoove together
variant_calling/smoove/smoove_svmerged.clustered.vcf: variant_calling/smoove/smoove_svmerge.sh \
	refgenome/Gmax_508_v4.0_mit_chlp.fasta \
	variant_calling/smoove/smoove_filtered.vcf
	$<

variant_calling/smoove/smoove_filtered.vcf: variant_calling/smoove/smoove_filter.sh \
	refgenome/Gmax_508_v4.0_mit_chlp.fasta \
	refgenome/Gmax_508_v4.0_mit_chlp.fasta.fai \
	variant_calling/smoove/all_samples.smoove.square.vcf.gz \
	variant_calling/smoove/ACO_header_line.txt \
	scripts/extract_svs_50.awk
	$<

# Calling variants and genotyping them on the whole population with smoove
variant_calling/smoove/all_samples.smoove.square.vcf.gz: variant_calling/smoove/smoove_call.sh \
	refgenome/Gmax_508_v4.0_mit_chlp.fasta \
	utilities/srr_id_correspondence.txt \
	illumina_data/merged_bams/ILLUMINA_BAM_MERGING
	$<

# CALLING SVS WITH MANTA --------------------------------------------------

# Merging the filtered variants called with Manta using SVmerge
variant_calling/manta/manta_svmerged.clustered.vcf: variant_calling/manta/manta_svmerge.sh \
	refgenome/Gmax_508_v4.0_mit_chlp.fasta \
	variant_calling/manta/manta_filtered.vcf
	$<

variant_calling/manta/manta_filtered.vcf: variant_calling/manta/manta_filter.sh \
	variant_calling/manta/MANTA_CALLING \
	refgenome/Gmax_508_v4.0_mit_chlp.fasta \
	variant_calling/manta/ACO_header_line.txt \
	scripts/extract_svs_50.awk
	$<

variant_calling/manta/MANTA_CALLING: variant_calling/manta/manta_call.sh \
	variant_calling/manta/manta_config.txt \
	variant_calling/manta/config_commands.sh \
	refgenome/Gmax_508_v4.0_mit_chlp.fasta \
	illumina_data/merged_bams/ILLUMINA_BAM_MERGING
	$<

# CALLING SVS WITH SVABA --------------------------------------------------

variant_calling/svaba/svaba_svmerged.clustered.vcf: variant_calling/svaba/svaba_svmerge.sh \
	refgenome/Gmax_508_v4.0_mit_chlp.fasta \
	variant_calling/svaba/svaba_filtered.vcf
	$<

variant_calling/svaba/svaba_filtered.vcf: variant_calling/svaba/svaba_filter.sh \
	refgenome/Gmax_508_v4.0_mit_chlp.fasta \
	refgenome/Gmax_508_v4.0_mit_chlp.fasta.fai \
	utilities/srr_id_correspondence.txt \
	variant_calling/svaba/convert_svaba.R \
	scripts/svaba_process.R \
	variant_calling/svaba/SVABA_CALLING \
	variant_calling/svaba/annotate_svtype.awk \
	variant_calling/svaba/ACO_header_line.txt \
	scripts/extract_svs_50.awk
	$<

variant_calling/svaba/SVABA_CALLING: variant_calling/svaba/svaba_call.sh \
	refgenome/Gmax_508_v4.0_mit_chlp.fasta \
	utilities/srr_id_correspondence.txt \
	illumina_data/merged_bams/ILLUMINA_BAM_MERGING
	$<

# CALLING SVS WITH WITH SVMU BASED ON DE NOVO ASSEMBLIES --------------------------------------------------

# Merging the SVs found for all genome assemblies
variant_calling/assemblies/assembly_svs.vcf: variant_calling/assemblies/assembly_svmerge.sh \
	refgenome/Gmax_508_v4.0_mit_chlp.fasta \
	variant_calling/assemblies/merging_files.txt \
	variant_calling/assemblies/ASSEMBLY_FILTERING
	$<

# Converting from svmu's sv.txt to VCF format
variant_calling/assemblies/ASSEMBLY_FILTERING: variant_calling/assemblies/assembly_filter.sh \
	refgenome/Gmax_508_v4.0_mit_chlp.fasta \
	variant_calling/assemblies/assembly_samples.txt \
	$(shell cut -d " " -f1 variant_calling/assemblies/assembly_samples.txt | xargs -I {} echo external_data/genome_assemblies/{}) \
	$(shell cut -d " " -f1 variant_calling/assemblies/assembly_samples.txt | xargs -I {} echo external_data/genome_assemblies/{}.fai) \
	variant_calling/assemblies/svmu_to_vcf.R \
	variant_calling/assemblies/SVMU_CALLING
	$<

# Call variants from MuMMER alignments using svmu
variant_calling/assemblies/SVMU_CALLING: variant_calling/assemblies/svmu_call.sh \
	variant_calling/assemblies/assembly_samples.txt \
	variant_calling/assemblies/MUMMER_ALIGNMENT \
	$(shell cut -d " " -f1 variant_calling/assemblies/assembly_samples.txt | xargs -I {} echo external_data/genome_assemblies/{}) \
	$(shell cut -d " " -f1 variant_calling/assemblies/assembly_samples.txt | xargs -I {} echo external_data/genome_assemblies/{}.fai) \
	refgenome/Gmax_508_v4.0_mit_chlp.fasta
	$<

# Align whole-genome assemblies to the reference Williams82 genome v.4 using MuMMER
variant_calling/assemblies/MUMMER_ALIGNMENT: variant_calling/assemblies/mummer.sh \
	variant_calling/assemblies/assembly_samples.txt \
	$(shell cut -d " " -f1 variant_calling/assemblies/assembly_samples.txt | xargs -I {} echo external_data/genome_assemblies/{}) \
	$(shell cut -d " " -f1 variant_calling/assemblies/assembly_samples.txt | xargs -I {} echo external_data/genome_assemblies/{}.fai) \
	refgenome/Gmax_508_v4.0_mit_chlp.fasta
	$<

# GWAS ANALYSIS WITH K-MERS --------------------------------------------------

# Creating the .rds object with the filtered set of reads containing significant k-mers
gwas_results/kmers/%_kmer_positions.rds: gwas_results/kmers/combine_reads.R \
	gwas_results/kmers/%/katcher_results/KATCHER \
	refgenome/Gmax_508_v4.0_mit_chlp.fasta \
	refgenome/Gmax_508_v4.0_mit_chlp.fasta.fai
	$(RSCRIPT) $< $*

# Running katcher for each of the traits
gwas_results/kmers/%/katcher_results/KATCHER: gwas_results/kmers/katcher.sh \
	illumina_data/merged_bams/ILLUMINA_BAM_MERGING \
	utilities/srr_id_correspondence.txt \
	gwas_results/kmers/KMER_GWAS \
	kmers_table/kmers_table.names \
	kmers_table/kmers_table.table
	$< $*

# A target for running Voichek and Weigel's (2020) k-mer analysis
gwas_results/kmers/KMER_GWAS: gwas_results/kmers/kmer_gwas.R \
	utilities/trait_names.txt \
	phenotypic_data/phenotypic_data.csv \
	kmers_table/kmers_table.names \
	kmers_table/kmers_table.table \
	kmers_table/kmers_table.kinship
	$(RSCRIPT) $<

# kmers_table/kmers_table.table kmers_table/kmers_table.names:
# kmers_table/kmers_table.kinship:






# filtered_variants/$(prog)/filtered_variants.vcf.gz:


# PREPARING THE PHENOTYPIC DATA FOR GWAS
# phenotypic_data/phenotypic_data.csv:
