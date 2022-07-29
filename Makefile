RSCRIPT := ~/.local/bin/Rscript
PDFLATEX := ~/.local/texlive/2022/bin/x86_64-linux/pdflatex

SDIR := additional_files

SUPFIGURES := figures/stem_termination_all_manhattan.png \
	figures/stem_termination_sn_manhattan.png \
	figures/pubescence_color_all_manhattan.png \
	figures/pubescence_color_nogray_manhattan.png \
	figures/protein_manhattan.png \
	figures/pubescence_form_all_manhattan.png \
	figures/pubescence_form_noerect_manhattan.png \
	figures/oil_manhattan.png \
	figures/pubescence_density_manhattan.png \
	figures/seed_coat_luster_all_manhattan.png \
	figures/seed_coat_luster_nointermediate_manhattan.png \
	figures/seed_coat_luster_dullshiny_manhattan.png \
	figures/pod_color_all_manhattan.png \
	figures/pod_color_blbr_manhattan.png \
	figures/flower_color_manhattan.png \
	figures/seed_coat_color_all_manhattan.png \
	figures/seed_coat_color_greenyellow_manhattan.png \
	figures/hilum_color_all_manhattan.png \
	figures/hilum_color_blackbrown_manhattan.png \
	figures/hilum_color_rbr_manhattan.png \
	figures/maturity_group_manhattan.png \
	figures/corrected_dry_weight_manhattan.png \
	figures/fc_w_signal.png \
	figures/pca_t_signal.png \
	figures/sda_dt1_signal.png

figures/%_manhattan.png: figures/manhattan_plot.R \
	figures/ggplots/platypus_%_manhattan.rds \
	figures/ggplots/vg_%_manhattan.rds \
	figures/ggplots/paragraph_%_manhattan.rds \
	figures/ggplots/kmers_%_manhattan.rds
	$(RSCRIPT) figures/manhattan_plot.R $*

.PRECIOUS: figures/ggplots/platypus_%_manhattan.rds figures/ggplots/vg_%_manhattan.rds \
	figures/ggplots/paragraph_%_manhattan.rds figures/ggplots/kmers_%_manhattan.rds

figures/ggplots/platypus_%_manhattan.rds: figures/manhattan_subplot.R \
	refgenome/Gmax_508_v4.0_mit_chlp.fasta \
	gwas_results/platypus/%_gwas.csv \
	gwas_results/platypus/%_threshold_5per.txt
	$(RSCRIPT) figures/manhattan_subplot.R $* platypus

figures/ggplots/vg_%_manhattan.rds: figures/manhattan_subplot.R \
	refgenome/Gmax_508_v4.0_mit_chlp.fasta \
	gwas_results/vg/%_gwas.csv \
	gwas_results/vg/%_threshold_5per.txt
	$(RSCRIPT) figures/manhattan_subplot.R $* vg

figures/ggplots/paragraph_%_manhattan.rds: figures/manhattan_subplot.R \
	refgenome/Gmax_508_v4.0_mit_chlp.fasta \
	gwas_results/paragraph/%_gwas.csv \
	gwas_results/paragraph/%_threshold_5per.txt
	$(RSCRIPT) figures/manhattan_subplot.R $* paragraph

figures/ggplots/kmers_%_manhattan.rds: figures/manhattan_subplot.R \
	gwas_results/kmers/%_kmer_positions.rds \
	gwas_results/kmers/%_threshold_5per \
	refgenome/Gmax_508_v4.0_mit_chlp.fasta
	$(RSCRIPT) figures/manhattan_subplot.R $* kmers

# Compiling the Supplemental Data file from the .tex file
$(SDIR)/additional_file_1.pdf: $(SDIR)/additional_file_1.tex $(SUPFIGURES)
	cd $(SDIR) ; $(PDFLATEX) additional_file_1.tex

figures/%_signal.png: figures/signal_plot.R \
	utilities/loci.txt \
	refgenome/Gmax_508_v4.0_mit_chlp.fasta \
	refgenome/lookup_v1_to_v4.rds \
	refgenome/gmax_v4_genes.rds \
	refgenome/gmax_v4_transcripts.rds \
	refgenome/gmax_v4_exons.rds \
	refgenome/gmax_v4_cds.rds
	$(RSCRIPT) figures/signal_plot.R $*

