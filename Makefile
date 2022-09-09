RSCRIPT := ~/.local/bin/Rscript --quiet
PDFLATEX := ~/.local/texlive/2022/bin/x86_64-linux/pdflatex

SDIR := additional_files

# Creating a few variables for improving the readability of rules
txdb := refgenome/gmax_v4_genes.rds refgenome/gmax_v4_transcripts.rds refgenome/gmax_v4_exons.rds refgenome/gmax_v4_cds.rds
refgen := refgenome/Gmax_508_v4.0_mit_chlp.fasta
signals_gr := utilities/all_signals.rds
grobdir := figures/grobs

suptables := tables/FLOWER.COLOR_gwas_table.csv \
	tables/PUBESCENCE.COLOR_gwas_table.csv \
	tables/SEED.COAT.COLOR_gwas_table.csv \
	tables/STEM.TERMINATION.TYPE_gwas_table.csv \
	tables/HILUM.COLOR_gwas_table.csv \
	tables/POD.COLOR_gwas_table.csv \
	tables/PUBESCENCE.FORM_gwas_table.csv \
	tables/PUBESCENCE.DENSITY_gwas_table.csv \
	tables/SEED.COAT.LUSTER_gwas_table.csv \
	tables/MATURITY.GROUP_gwas_table.csv

supfigures := $(shell cat utilities/trait_names.txt | xargs -I {} echo figures/{}_manhattan.png) \
	$(shell cut -d "," -f1 utilities/signal_ids.txt | xargs -I {} echo figures/{}_signal.png) \
	$(shell cat utilities/gene_signal_ids.txt | xargs -I {} echo figures/{}_gene.png) \
	$(shell cut -f1 utilities/kmer_plot_ranges.txt | xargs -I {} echo figures/{}_kmers.png)

all: $(SDIR)/additional_file_1.pdf

SUPTABLES: $(suptables)
SUPFIGURES: $(supfigures)

# Making sure that intermediate files do not get deleted
.PRECIOUS: $(grobdir)/platypus_%_manhattan.rds $(grobdir)/vg_%_manhattan.rds \
	$(grobdir)/paragraph_%_manhattan.rds $(grobdir)/kmers_%_manhattan.rds \
	$(grobdir)/platypus_%_signal.rds $(grobdir)/vg_%_signal.rds \
	$(grobdir)/paragraph_%_signal.rds $(grobdir)/kmers_%_signal.rds \
	$(grobdir)/platypus_%_gene.rds $(grobdir)/vg_%_gene.rds \
	$(grobdir)/paragraph_%_gene.rds $(grobdir)/kmers_%_gene.rds \
	gwas_results/kmers/%_threshold_5per.txt \
	gwas_results/kmer_consensus/%_sequences.fa

# Adding some more intermediate files to .PRECIOUS
$(foreach signal,$(shell cut -d "," -f1 utilities/signal_ids.txt | xargs -I {} echo gwas_results/%/{}_gwas_locus.rds),$(eval .PRECIOUS: $(signal)))
$(foreach signal,$(shell cut -d "," -f1 utilities/signal_ids.txt | xargs -I {} echo gwas_results/%/{}_signal_locus.rds),$(eval .PRECIOUS: $(signal)))
$(foreach prog,vg platypus paragraph kmers,$(eval .PRECIOUS : gwas_results/$(prog)/%_signal.rds))
$(foreach prog,vg platypus paragraph kmers,$(eval .PRECIOUS : gwas_results/$(prog)/%_gwas.rds))

# Compiling the Supplemental Data file from the .tex file
$(SDIR)/additional_file_1.pdf: $(SDIR)/additional_file_1.tex $(supfigures) $(suptables)
	cd $(SDIR) ; $(PDFLATEX) additional_file_1.tex ; $(PDFLATEX) additional_file_1.tex

# SIGNALS --------------------------------------------------

# This script prepares the reference signals GRanges object
$(signals_gr) utilities/signal_ids.txt utilities/gene_signal_ids.txt: utilities/make_signals_granges.R \
	reference_signals/bandillo2017_signals.rds \
	reference_signals/bandillo2015_signals.rds \
	reference_signals/deronne2022_signal.rds \
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

# Putting the k-mer 5% thresholds on the same scale as the thresholds for the other programs
gwas_results/kmers/%_threshold_5per.txt: gwas_results/scale_kmer_thresholds.R gwas_results/kmers/%_threshold_5per
	$(RSCRIPT) gwas_results/scale_kmer_thresholds.R $*

# PHENOTYPIC DATA --------------------------------------------------

# Creating the table with the numeric coding for a given trait
tables/%_gwas_table.csv: tables/gwas_table.R \
	phenotypic_data/trait_names.rds \
	phenotypic_data/phenotypic_data.csv \
	phenotypic_data/lookup_tables.rds
	$(RSCRIPT) tables/gwas_table.R $*

# Creating the lookup tables used to translate the phenotypes into numeric codes
phenotypic_data/lookup_tables.rds: phenotypic_data/lookup_tables.R
	$(RSCRIPT) phenotypic_data/lookup_tables.R

# Creating the lookup table used to translate the GRIN trait names into the ones used in this analysis
phenotypic_data/trait_names.rds: phenotypic_data/trait_names.R
	$(RSCRIPT) phenotypic_data/trait_names.R

# KMER PLOTS --------------------------------------------------
# Generating the k-mer plot from the consensus sequences
figures/%_kmers.png: figures/kmer_plot.R \
	gwas_results/kmer_consensus/%_sequences.fa
	$(RSCRIPT) figures/kmer_plot.R $*

# Extracting the consensus sequences from significant k-mer assemblies for a given locus
gwas_results/kmer_consensus/%_sequences.fa: gwas_results/gather_consensus.sh \
	utilities/kmer_plot_ranges.txt \
	refgenome/Gmax_508_v4.0_mit_chlp.fasta
	gwas_results/gather_consensus.sh $*

# GWAS RESULTS --------------------------------------------------

# Preparing the GWAS results for each of platypus, vg and paragraph
$(foreach prog,platypus vg paragraph,$(eval gwas_results/$(prog)/%_gwas.rds: gwas_results/format_gwas_results.R \
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
# Creating a GRanges object containing the signal(s) for all traits for each program but platypus; also creating subsets of GWAS results that only overlap signals
$(foreach prog,vg paragraph kmers,$(eval gwas_results/$(prog)/%_signal.rds gwas_results/$(prog)/%_gwas_subset.rds : gwas_results/find_signals.R \
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

# MANHATTAN PLOTS --------------------------------------------------

# This block assembles the four manhattan subplots for a given trait
figures/%_manhattan.png: figures/manhattan_plot.R \
	$(grobdir)/platypus_%_manhattan.rds \
	$(grobdir)/vg_%_manhattan.rds \
	$(grobdir)/paragraph_%_manhattan.rds \
	$(grobdir)/kmers_%_manhattan.rds
	$(RSCRIPT) figures/manhattan_plot.R $*

# Preparing the manhattan subplots for each of platypus, vg and paragraph
$(foreach prog,platypus vg paragraph kmers,$(eval $(grobdir)/$(prog)_%_manhattan.rds: figures/manhattan_subplot.R \
	$(signals_gr) \
	gwas_results/$(prog)/%_gwas.rds \
	gwas_results/$(prog)/%_signal.rds \
	gwas_results/$(prog)/%_threshold_5per.txt ; \
	$(RSCRIPT) figures/manhattan_subplot.R $$* $(prog)))

# SIGNAL PLOTS --------------------------------------------------

# This block assembles the four signal subplots for a given locus
figures/%_signal.png: figures/signal_plot.R \
	$(signals_gr) \
	$(txdb) \
	$(grobdir)/platypus_%_signal.rds \
	$(grobdir)/vg_%_signal.rds \
	$(grobdir)/paragraph_%_signal.rds \
	$(grobdir)/kmers_%_signal.rds
	$(RSCRIPT) figures/signal_plot.R $*

# Preparing the signal subplots for each of platypus, vg and paragraph
$(foreach prog,platypus vg paragraph kmers,$(eval $(grobdir)/$(prog)_%_signal.rds: figures/signal_subplot.R \
	$(signals_gr) \
	refgenome/gmax_v4_genes.rds \
	gwas_results/$(prog)/%_gwas_locus.rds \
	gwas_results/$(prog)/%_signal_locus.rds ; \
	$(RSCRIPT) figures/signal_subplot.R $$* $(prog)))

# GENE PLOTS --------------------------------------------------

# This block assembles the four gene subplots for a given locus
figures/%_gene.png: figures/gene_plot.R \
	$(signals_gr) \
	$(txdb) \
	$(grobdir)/platypus_%_gene.rds \
	$(grobdir)/vg_%_gene.rds \
	$(grobdir)/paragraph_%_gene.rds \
	$(grobdir)/kmers_%_gene.rds
	$(RSCRIPT) figures/gene_plot.R $*

# Preparing the gene subplots for each of platypus, vg and paragraph
$(foreach prog,platypus vg paragraph kmers,$(eval $(grobdir)/$(prog)_%_gene.rds: figures/gene_subplot.R \
	$(signals_gr) \
	refgenome/gmax_v4_genes.rds \
	gwas_results/$(prog)/%_gwas_locus.rds \
	gwas_results/$(prog)/%_signal_locus.rds ; \
	$(RSCRIPT) figures/gene_subplot.R $$* $(prog)))

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

