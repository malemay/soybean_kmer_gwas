RSCRIPT := ~/.local/bin/Rscript --quiet
PDFLATEX := ~/.local/texlive/2022/bin/x86_64-linux/pdflatex

SDIR := additional_files

# Creating a few variables for improving the readability of rules
txdb := refgenome/gmax_v4_genes.rds refgenome/gmax_v4_transcripts.rds refgenome/gmax_v4_exons.rds refgenome/gmax_v4_cds.rds
refgen := refgenome/Gmax_508_v4.0_mit_chlp.fasta
signals_gr := utilities/all_signals.rds
grobdir := figures/grobs
addfile := additional_files/additional_file_1.tex

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

# Getting the LD png figures to generate dynamically from additional_files/additional_file_1.tex
manhattanplots := $(shell grep '^\\manhattanplot' $(addfile) | grep -v scaffolds | sed -E 's/\\manhattanplot\{(.*)\}/\1/' | xargs -I {} echo figures/{}_manhattan.png)
scaffoldplots := $(shell grep '^\\manhattanplot' $(addfile) | grep scaffolds | sed -E 's/\\manhattanplot\{(.*)\}/\1/' | xargs -I {} echo figures/{}_manhattan.png)
ldplots := $(shell grep '^\\ldplot' $(addfile) | sed -E 's/\\ldplot\{(.*)\}/\1/' | xargs -I {} echo figures/{}_ld.png)
geneplots := $(shell grep '^\\geneplot' $(addfile) | sed -E 's/\\geneplot\{(.*)\}/\1/' | xargs -I {} echo figures/{}_gene.png)
kmerplots := $(shell grep '^\\kmerplot' $(addfile) | sed -E 's/\\kmerplot\{(.*)\}/\1/' | xargs -I {} echo figures/{}_kmers.png)
signalplots := $(shell grep '^\\signalplot' $(addfile) | sed -E 's/\\signalplot\{(.*)\}/\1/' | xargs -I {} echo figures/{}_signal.png)

supfigures := $(manhattanplots) \
	$(signalplots) \
	$(geneplots) \
	$(kmerplots) \
	$(scaffoldplots) \
	$(ldplots)


topgranges := $(foreach prog,platypus vg paragraph kmers,$(shell cut -d "," -f1 utilities/signal_ids.txt | xargs -I {} echo gwas_results/$(prog)/{}_top_markers.rds))

allgenes := $(foreach prog,platypus vg paragraph kmers,$(shell cut -d "," -f1 utilities/signal_ids.txt | xargs -I {} echo gwas_results/$(prog)/{}_all_genes.tsv))

topgenes := $(foreach prog,platypus vg paragraph kmers,$(shell cut -d "," -f1 utilities/signal_ids.txt | xargs -I {} echo gwas_results/$(prog)/{}_top_genes.tsv))

nearestgene := $(foreach prog,platypus vg paragraph kmers,$(shell cut -d "," -f1 utilities/signal_ids.txt | xargs -I {} echo gwas_results/$(prog)/{}_nearest_gene.tsv))

genetables := $(allgenes) $(topgenes) $(nearestgene)

phenodata := $(shell cut -f1 utilities/kmer_plot_ranges.txt | xargs -I {} echo gwas_results/kmers/{}_phenodata.rds)


all: $(SDIR)/additional_file_1.pdf $(genetables) $(phenodata)

GENETABLES : $(genetables)
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
	gwas_results/kmer_consensus/%_sequences.fa \
	gwas_results/kmers/%_clustered_ld.txt \
	cnv_analysis/hps_cnv_range.rds \
	cnv_analysis/ps_cnv_range.rds

# Adding some more intermediate files to .PRECIOUS
$(foreach signal,$(shell cut -d "," -f1 utilities/signal_ids.txt | xargs -I {} echo gwas_results/%/{}_gwas_locus.rds),$(eval .PRECIOUS: $(signal)))
$(foreach signal,$(shell cut -d "," -f1 utilities/signal_ids.txt | xargs -I {} echo gwas_results/%/{}_signal_locus.rds),$(eval .PRECIOUS: $(signal)))
$(foreach prog,vg platypus paragraph kmers,$(eval .PRECIOUS : gwas_results/$(prog)/%_signal.rds))
$(foreach prog,vg platypus paragraph kmers,$(eval .PRECIOUS : gwas_results/$(prog)/%_gwas.rds))

.PRECIOUS : $(topgranges)

# Compiling the Supplemental Data file from the .tex file
$(SDIR)/additional_file_1.pdf: $(SDIR)/additional_file_1.tex $(supfigures) $(suptables)
	cd $(SDIR) ; $(PDFLATEX) additional_file_1.tex ; $(PDFLATEX) additional_file_1.tex

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
gwas_results/kmers/%_threshold_5per.txt: gwas_results/scale_kmer_thresholds.R gwas_results/kmers/%_threshold_5per
	$(RSCRIPT) gwas_results/scale_kmer_thresholds.R $*

# Generating the .csv file with data regarding which reference signals were found and their p-values
tables/signals_table.csv: tables/signals_table.R \
	$(signals_gr) \
	$(shell cat utilities/trait_names.txt | xargs -I {} echo gwas_results/platypus/{}_signal.rds) \
	$(shell cat utilities/trait_names.txt | xargs -I {} echo gwas_results/vg/{}_signal.rds) \
	$(shell cat utilities/trait_names.txt | xargs -I {} echo gwas_results/paragraph/{}_signal.rds) \
	$(shell cat utilities/trait_names.txt | xargs -I {} echo gwas_results/kmers/{}_signal.rds)
	$(RSCRIPT) tables/signals_table.R

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

# KMER HAPLOTYPE PLOTS --------------------------------------------------
# Generating the k-mer plot from the consensus sequences (the k-mer p-values are missing from the list of dependencies)
figures/%_kmers.png gwas_results/kmers/%_phenodata.rds: figures/kmer_plot.R \
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
# Creating a GRanges object containing the signal(s) for all traits for each program but platypus
# Also creating subsets of GWAS results that only overlap signals
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

# GENE ANALYSIS  --------------------------------------------------
# Generating a GRanges object and tsv files for each locus
$(foreach prog,platypus vg paragraph kmers,$(eval \
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

# Manhattan plots for unanchored scaffolds (for the k-mers only)
figures/%_scaffolds_manhattan.png: figures/scaffolds_manhattan.R \
	gwas_results/kmers/%_gwas.rds \
	gwas_results/kmers/%_threshold_5per.txt
	$(RSCRIPT) figures/scaffolds_manhattan.R $* kmers

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
	gwas_results/$(prog)/%_top_markers.rds \
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
	utilities/cnv_genes.txt \
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

