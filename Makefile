RSCRIPT := ~/.local/bin/Rscript --quiet
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
	figures/stem_termination_all_Dt2_signal.png \
	figures/stem_termination_all_Dt1_signal.png \
	figures/stem_termination_all_E3_signal.png \
	figures/stem_termination_sn_Dt2_signal.png \
	figures/stem_termination_sn_Dt1_signal.png \
	figures/stem_termination_sn_E3_signal.png \
	figures/pubescence_color_all_Td_signal.png \
	figures/pubescence_color_all_T_signal.png \
	figures/pubescence_color_nogray_Td_signal.png \
	figures/pubescence_color_nogray_T_signal.png \
	figures/pubescence_form_all_Pa1_signal.png \
	figures/pubescence_form_all_Pa2_signal.png \
	figures/pubescence_form_noerect_Pa1_signal.png \
	figures/pubescence_form_noerect_Pa2_signal.png \
	figures/pubescence_density_Pd1_signal.png \
	figures/pubescence_density_Ps_signal.png \
	figures/pubescence_density_P1_signal.png \
	figures/seed_coat_luster_all_I_signal.png \
	figures/seed_coat_luster_all_Bq_signal.png \
	figures/seed_coat_luster_all_B1_signal.png \
	figures/seed_coat_luster_all_Hps_signal.png \
	figures/seed_coat_luster_all_B_signal.png \
	figures/seed_coat_luster_nointermediate_I_signal.png \
	figures/seed_coat_luster_nointermediate_Bq_signal.png \
	figures/seed_coat_luster_nointermediate_B1_signal.png \
	figures/seed_coat_luster_nointermediate_Hps_signal.png \
	figures/seed_coat_luster_nointermediate_B_signal.png \
	figures/seed_coat_luster_dullshiny_I_signal.png \
	figures/seed_coat_luster_dullshiny_Bq_signal.png \
	figures/seed_coat_luster_dullshiny_B1_signal.png \
	figures/seed_coat_luster_dullshiny_Hps_signal.png \
	figures/seed_coat_luster_dullshiny_B_signal.png \
	figures/pod_color_all_L2_signal.png \
	figures/pod_color_all_L1_signal.png \
	figures/pod_color_blbr_L2_signal.png \
	figures/pod_color_blbr_L1_signal.png \
	figures/flower_color_W1_signal.png \
	figures/flower_color_L1_signal.png \
	figures/seed_coat_color_all_G_signal.png \
	figures/seed_coat_color_all_T_signal.png \
	figures/seed_coat_color_all_O_signal.png \
	figures/seed_coat_color_all_I_signal.png \
	figures/seed_coat_color_all_R_signal.png \
	figures/seed_coat_color_greenyellow_G_signal.png \
	figures/seed_coat_color_greenyellow_T_signal.png \
	figures/seed_coat_color_greenyellow_O_signal.png \
	figures/seed_coat_color_greenyellow_I_signal.png \
	figures/seed_coat_color_greenyellow_R_signal.png \
	figures/hilum_color_all_T_signal.png \
	figures/hilum_color_all_O_signal.png \
	figures/hilum_color_all_I_signal.png \
	figures/hilum_color_all_R_signal.png \
	figures/hilum_color_all_W1_signal.png \
	figures/hilum_color_blackbrown_T_signal.png \
	figures/hilum_color_blackbrown_O_signal.png \
	figures/hilum_color_blackbrown_I_signal.png \
	figures/hilum_color_blackbrown_R_signal.png \
	figures/hilum_color_blackbrown_W1_signal.png \
	figures/hilum_color_rbr_T_signal.png \
	figures/hilum_color_rbr_O_signal.png \
	figures/hilum_color_rbr_I_signal.png \
	figures/hilum_color_rbr_R_signal.png \
	figures/hilum_color_rbr_W1_signal.png \
	figures/maturity_group_E1_signal.png \
	figures/maturity_group_E2_signal.png \
	figures/maturity_group_E3_signal.png \
	figures/maturity_group_E4_signal.png \
	figures/stem_termination_all_Dt2_gene.png \
	figures/stem_termination_all_Dt1_gene.png \
	figures/stem_termination_all_E3_gene.png \
	figures/stem_termination_sn_Dt2_gene.png \
	figures/stem_termination_sn_Dt1_gene.png \
	figures/stem_termination_sn_E3_gene.png \
	figures/pubescence_color_all_T_gene.png \
	figures/pubescence_color_nogray_T_gene.png \
	figures/pubescence_density_P1_gene.png \
	figures/seed_coat_luster_all_I_gene.png \
	figures/seed_coat_luster_all_Hps_gene.png \
	figures/seed_coat_luster_nointermediate_I_gene.png \
	figures/seed_coat_luster_nointermediate_Hps_gene.png \
	figures/seed_coat_luster_dullshiny_I_gene.png \
	figures/seed_coat_luster_dullshiny_Hps_gene.png \
	figures/pod_color_all_L1_gene.png \
	figures/pod_color_blbr_L1_gene.png \
	figures/flower_color_W1_gene.png \
	figures/flower_color_L1_gene.png \
	figures/seed_coat_color_all_T_gene.png \
	figures/seed_coat_color_all_I_gene.png \
	figures/seed_coat_color_all_R_gene.png \
	figures/seed_coat_color_greenyellow_T_gene.png \
	figures/seed_coat_color_greenyellow_I_gene.png \
	figures/seed_coat_color_greenyellow_R_gene.png \
	figures/hilum_color_all_T_gene.png \
	figures/hilum_color_all_I_gene.png \
	figures/hilum_color_all_R_gene.png \
	figures/hilum_color_all_W1_gene.png \
	figures/hilum_color_blackbrown_T_gene.png \
	figures/hilum_color_blackbrown_I_gene.png \
	figures/hilum_color_blackbrown_R_gene.png \
	figures/hilum_color_blackbrown_W1_gene.png \
	figures/hilum_color_rbr_T_gene.png \
	figures/hilum_color_rbr_I_gene.png \
	figures/hilum_color_rbr_R_gene.png \
	figures/hilum_color_rbr_W1_gene.png \
	figures/maturity_group_E1_gene.png \
	figures/maturity_group_E2_gene.png \
	figures/maturity_group_E3_gene.png \
	figures/maturity_group_E4_gene.png

all: $(SDIR)/additional_file_1.pdf

.PRECIOUS: figures/grobs/platypus_%_manhattan.rds figures/grobs/vg_%_manhattan.rds \
	figures/grobs/paragraph_%_manhattan.rds figures/grobs/kmers_%_manhattan.rds \
	figures/grobs/platypus_%_signal.rds figures/grobs/vg_%_signal.rds \
	figures/grobs/paragraph_%_signal.rds figures/grobs/kmers_%_signal.rds \
	figures/grobs/platypus_%_gene.rds figures/grobs/vg_%_gene.rds \
	figures/grobs/paragraph_%_gene.rds figures/grobs/kmers_%_gene.rds \
	gwas_results/platypus/%_gwas.rds gwas_results/vg/%_gwas.rds \
	gwas_results/paragraph/%_gwas.rds gwas_results/kmers/%_gwas.rds

# Compiling the Supplemental Data file from the .tex file
$(SDIR)/additional_file_1.pdf: $(SDIR)/additional_file_1.tex $(SUPFIGURES)
	cd $(SDIR) ; $(PDFLATEX) additional_file_1.tex

# The next four blocks compute the GWAS results and output them to .rds format for each SV genotyping approach
# This creates a unified framework for downstream analyses
gwas_results/platypus/%_gwas.rds: gwas_results/format_gwas_results.R \
	refgenome/Gmax_508_v4.0_mit_chlp.fasta \
	utilities/signal_ids.txt \
	gwas_results/platypus/%_gwas.csv 
	$(RSCRIPT) gwas_results/format_gwas_results.R $* platypus

gwas_results/vg/%_gwas.rds: gwas_results/format_gwas_results.R \
	refgenome/Gmax_508_v4.0_mit_chlp.fasta \
	utilities/signal_ids.txt \
	gwas_results/vg/%_gwas.csv 
	$(RSCRIPT) gwas_results/format_gwas_results.R $* vg

gwas_results/paragraph/%_gwas.rds: gwas_results/format_gwas_results.R \
	refgenome/Gmax_508_v4.0_mit_chlp.fasta \
	utilities/signal_ids.txt \
	gwas_results/paragraph/%_gwas.csv 
	$(RSCRIPT) gwas_results/format_gwas_results.R $* paragraph

gwas_results/kmers/%_gwas.rds: gwas_results/format_gwas_results.R \
	utilities/signal_ids.txt \
	refgenome/Gmax_508_v4.0_mit_chlp.fasta \
	gwas_results/kmers/%_kmer_positions.rds 
	$(RSCRIPT) gwas_results/format_gwas_results.R $* kmers

# This block assembles the four manhattan subplots for a given trait
figures/%_manhattan.png: figures/manhattan_plot.R \
	figures/grobs/platypus_%_manhattan.rds \
	figures/grobs/vg_%_manhattan.rds \
	figures/grobs/paragraph_%_manhattan.rds \
	figures/grobs/kmers_%_manhattan.rds
	$(RSCRIPT) figures/manhattan_plot.R $*

# The next four blocks pertain to the sub-manhattan plots of each SV genotyping approach
figures/grobs/platypus_%_manhattan.rds: figures/manhattan_subplot.R \
	utilities/all_signals.rds \
	gwas_results/platypus/%_gwas.rds \
	gwas_results/platypus/%_threshold_5per.txt
	$(RSCRIPT) figures/manhattan_subplot.R $* platypus

figures/grobs/vg_%_manhattan.rds: figures/manhattan_subplot.R \
	utilities/all_signals.rds \
	gwas_results/vg/%_gwas.rds \
	gwas_results/vg/%_threshold_5per.txt
	$(RSCRIPT) figures/manhattan_subplot.R $* vg

figures/grobs/paragraph_%_manhattan.rds: figures/manhattan_subplot.R \
	utilities/all_signals.rds \
	gwas_results/paragraph/%_gwas.rds \
	gwas_results/paragraph/%_threshold_5per.txt
	$(RSCRIPT) figures/manhattan_subplot.R $* paragraph

figures/grobs/kmers_%_manhattan.rds: figures/manhattan_subplot.R \
	utilities/all_signals.rds \
	gwas_results/kmers/%_gwas.rds \
	gwas_results/kmers/%_threshold_5per
	$(RSCRIPT) figures/manhattan_subplot.R $* kmers

# This script prepares the reference signals GRanges object
utilities/all_signals.rds utilities/signal_ids.txt utilities/gene_signal_ids.txt: utilities/make_signals_granges.R \
	utilities/bandillo_signals.rds \
	utilities/trait_names.txt
	$(RSCRIPT) utilities/make_signals_granges.R

# This block assembles the four signal subplots for a given locus
figures/%_signal.png: figures/signal_plot.R \
	figures/grobs/platypus_%_signal.rds \
	figures/grobs/vg_%_signal.rds \
	figures/grobs/paragraph_%_signal.rds \
	figures/grobs/kmers_%_signal.rds
	$(RSCRIPT) figures/signal_plot.R $*

# The next four blocks prepare the sub-plots of the signal plots for each SV genotyping approach
figures/grobs/platypus_%_signal.rds: figures/signal_subplot.R \
	utilities/all_signals.rds \
	refgenome/gmax_v4_genes.rds \
	refgenome/gmax_v4_transcripts.rds \
	refgenome/gmax_v4_exons.rds \
	refgenome/gmax_v4_cds.rds \
	gwas_results/platypus/%_locus_gwas.rds \
	gwas_results/platypus/%_locus_threshold_5per.txt
	$(RSCRIPT) figures/signal_subplot.R $* platypus

figures/grobs/vg_%_signal.rds: figures/signal_subplot.R \
	utilities/all_signals.rds \
	refgenome/gmax_v4_genes.rds \
	refgenome/gmax_v4_transcripts.rds \
	refgenome/gmax_v4_exons.rds \
	refgenome/gmax_v4_cds.rds \
	gwas_results/vg/%_locus_gwas.rds \
	gwas_results/vg/%_locus_threshold_5per.txt
	$(RSCRIPT) figures/signal_subplot.R $* vg

figures/grobs/paragraph_%_signal.rds: figures/signal_subplot.R \
	utilities/all_signals.rds \
	refgenome/gmax_v4_genes.rds \
	refgenome/gmax_v4_transcripts.rds \
	refgenome/gmax_v4_exons.rds \
	refgenome/gmax_v4_cds.rds \
	gwas_results/paragraph/%_locus_gwas.rds \
	gwas_results/paragraph/%_locus_threshold_5per.txt
	$(RSCRIPT) figures/signal_subplot.R $* paragraph

figures/grobs/kmers_%_signal.rds: figures/signal_subplot.R \
	utilities/all_signals.rds \
	refgenome/gmax_v4_genes.rds \
	refgenome/gmax_v4_transcripts.rds \
	refgenome/gmax_v4_exons.rds \
	refgenome/gmax_v4_cds.rds \
	gwas_results/kmers/%_locus_gwas.rds \
	gwas_results/kmers/%_locus_threshold_5per
	$(RSCRIPT) figures/signal_subplot.R $* kmers




# This block assembles the four gene subplots for a given locus
figures/%_gene.png: figures/gene_plot.R \
	figures/grobs/platypus_%_gene.rds \
	figures/grobs/vg_%_gene.rds \
	figures/grobs/paragraph_%_gene.rds \
	figures/grobs/kmers_%_gene.rds
	$(RSCRIPT) figures/gene_plot.R $*

# The next four blocks prepare the sub-plots of the gene plots for each SV genotyping approach
figures/grobs/platypus_%_gene.rds: figures/gene_subplot.R \
	utilities/all_signals.rds \
	refgenome/gmax_v4_genes.rds \
	refgenome/gmax_v4_transcripts.rds \
	refgenome/gmax_v4_exons.rds \
	refgenome/gmax_v4_cds.rds \
	gwas_results/platypus/%_locus_gwas.rds \
	gwas_results/platypus/%_locus_threshold_5per.txt
	$(RSCRIPT) figures/gene_subplot.R $* platypus

figures/grobs/vg_%_gene.rds: figures/gene_subplot.R \
	utilities/all_signals.rds \
	refgenome/gmax_v4_genes.rds \
	refgenome/gmax_v4_transcripts.rds \
	refgenome/gmax_v4_exons.rds \
	refgenome/gmax_v4_cds.rds \
	gwas_results/vg/%_locus_gwas.rds \
	gwas_results/vg/%_locus_threshold_5per.txt
	$(RSCRIPT) figures/gene_subplot.R $* vg

figures/grobs/paragraph_%_gene.rds: figures/gene_subplot.R \
	utilities/all_signals.rds \
	refgenome/gmax_v4_genes.rds \
	refgenome/gmax_v4_transcripts.rds \
	refgenome/gmax_v4_exons.rds \
	refgenome/gmax_v4_cds.rds \
	gwas_results/paragraph/%_locus_gwas.rds \
	gwas_results/paragraph/%_locus_threshold_5per.txt
	$(RSCRIPT) figures/gene_subplot.R $* paragraph

figures/grobs/kmers_%_gene.rds: figures/gene_subplot.R \
	utilities/all_signals.rds \
	refgenome/gmax_v4_genes.rds \
	refgenome/gmax_v4_transcripts.rds \
	refgenome/gmax_v4_exons.rds \
	refgenome/gmax_v4_cds.rds \
	gwas_results/kmers/%_locus_gwas.rds \
	gwas_results/kmers/%_locus_threshold_5per
	$(RSCRIPT) figures/gene_subplot.R $* kmers






