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
	figures/pca_t_signal.png \
	figures/fc_w_signal.png \
	figures/sda_dt1_signal.png \
	figures/hca_t_signal.png \
	figures/pfa_pa1_signal.png \
	figures/pfne_pa1_signal.png \
	figures/pbcng_td_signal.png \
	figures/oil_gm15_signal.png \
	figures/pd_ps_signal.png \
	figures/scla_b_signal.png \
	figures/sclni_b_signal.png \
	figures/sclds_b_signal.png \
	figures/pdca_l1_signal.png \
	figures/pdcbb_l1_signal.png \
	figures/scca_i_signal.png \
	figures/sccgy_g_signal.png \
	figures/hcbb_r_signal.png \
	figures/pca_t_gene.png \
	figures/fc_w_gene.png \
	figures/sda_dt1_gene.png \
	figures/hca_t_gene.png \
	figures/hcbb_r_gene.png #\
	figures/pdca_l1_gene.png \
	figures/pdcbb_l1_gene.png \

all: $(SDIR)/additional_file_1.pdf

.PRECIOUS: figures/ggplots/platypus_%_manhattan.rds figures/ggplots/vg_%_manhattan.rds \
	figures/ggplots/paragraph_%_manhattan.rds figures/ggplots/kmers_%_manhattan.rds \
	figures/ggplots/platypus_%_signal.rds figures/ggplots/vg_%_signal.rds \
	figures/ggplots/paragraph_%_signal.rds figures/ggplots/kmers_%_signal.rds \
	figures/ggplots/platypus_%_gene.rds figures/ggplots/vg_%_gene.rds \
	figures/ggplots/paragraph_%_gene.rds figures/ggplots/kmers_%_gene.rds

# Compiling the Supplemental Data file from the .tex file
$(SDIR)/additional_file_1.pdf: $(SDIR)/additional_file_1.tex $(SUPFIGURES)
	cd $(SDIR) ; $(PDFLATEX) additional_file_1.tex

# This block assembles the four manhattan subplots for a given trait
figures/%_manhattan.png: figures/manhattan_plot.R \
	figures/ggplots/platypus_%_manhattan.rds \
	figures/ggplots/vg_%_manhattan.rds \
	figures/ggplots/paragraph_%_manhattan.rds \
	figures/ggplots/kmers_%_manhattan.rds
	$(RSCRIPT) figures/manhattan_plot.R $*

# The next four blocks pertain to the sub-manhattan plots of each SV genotyping approach
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

# This block assembles the four signal subplots for a given locus
figures/%_signal.png: figures/signal_plot.R \
	figures/ggplots/platypus_%_signal.rds \
	figures/ggplots/vg_%_signal.rds \
	figures/ggplots/paragraph_%_signal.rds \
	figures/ggplots/kmers_%_signal.rds
	$(RSCRIPT) figures/signal_plot.R $*

# The next four blocks prepare the sub-plots of the signal plots for each SV genotyping approach
figures/ggplots/platypus_%_signal.rds: figures/signal_subplot.R \
	utilities/loci.txt \
	refgenome/Gmax_508_v4.0_mit_chlp.fasta \
	refgenome/lookup_v1_to_v4.rds \
	refgenome/gmax_v4_genes.rds \
	refgenome/gmax_v4_transcripts.rds \
	refgenome/gmax_v4_exons.rds \
	refgenome/gmax_v4_cds.rds \
	gwas_results/platypus/%_gwas.csv \
	gwas_results/platypus/%_threshold_5per.txt
	$(RSCRIPT) figures/signal_subplot.R $* platypus

figures/ggplots/vg_%_signal.rds: figures/signal_subplot.R \
	utilities/loci.txt \
	refgenome/Gmax_508_v4.0_mit_chlp.fasta \
	refgenome/lookup_v1_to_v4.rds \
	refgenome/gmax_v4_genes.rds \
	refgenome/gmax_v4_transcripts.rds \
	refgenome/gmax_v4_exons.rds \
	refgenome/gmax_v4_cds.rds \
	gwas_results/vg/%_gwas.csv \
	gwas_results/vg/%_threshold_5per.txt
	$(RSCRIPT) figures/signal_subplot.R $* vg

figures/ggplots/paragraph_%_signal.rds: figures/signal_subplot.R \
	utilities/loci.txt \
	refgenome/Gmax_508_v4.0_mit_chlp.fasta \
	refgenome/lookup_v1_to_v4.rds \
	refgenome/gmax_v4_genes.rds \
	refgenome/gmax_v4_transcripts.rds \
	refgenome/gmax_v4_exons.rds \
	refgenome/gmax_v4_cds.rds \
	gwas_results/paragraph/%_gwas.csv \
	gwas_results/paragraph/%_threshold_5per.txt
	$(RSCRIPT) figures/signal_subplot.R $* paragraph

figures/ggplots/kmers_%_signal.rds: figures/signal_subplot.R \
	utilities/loci.txt \
	refgenome/Gmax_508_v4.0_mit_chlp.fasta \
	refgenome/lookup_v1_to_v4.rds \
	refgenome/gmax_v4_genes.rds \
	refgenome/gmax_v4_transcripts.rds \
	refgenome/gmax_v4_exons.rds \
	refgenome/gmax_v4_cds.rds \
	gwas_results/kmers/%_kmer_positions.rds \
	gwas_results/kmers/%_threshold_5per
	$(RSCRIPT) figures/signal_subplot.R $* kmers




# This block assembles the four gene subplots for a given locus
figures/%_gene.png: figures/gene_plot.R \
	figures/ggplots/platypus_%_gene.rds \
	figures/ggplots/vg_%_gene.rds \
	figures/ggplots/paragraph_%_gene.rds \
	figures/ggplots/kmers_%_gene.rds
	$(RSCRIPT) figures/gene_plot.R $*

# The next four blocks prepare the sub-plots of the gene plots for each SV genotyping approach
figures/ggplots/platypus_%_gene.rds: figures/gene_subplot.R \
	utilities/loci.txt \
	refgenome/Gmax_508_v4.0_mit_chlp.fasta \
	refgenome/lookup_v1_to_v4.rds \
	refgenome/gmax_v4_genes.rds \
	refgenome/gmax_v4_transcripts.rds \
	refgenome/gmax_v4_exons.rds \
	refgenome/gmax_v4_cds.rds \
	gwas_results/platypus/%_gwas.csv \
	gwas_results/platypus/%_threshold_5per.txt
	$(RSCRIPT) figures/gene_subplot.R $* platypus

figures/ggplots/vg_%_gene.rds: figures/gene_subplot.R \
	utilities/loci.txt \
	refgenome/Gmax_508_v4.0_mit_chlp.fasta \
	refgenome/lookup_v1_to_v4.rds \
	refgenome/gmax_v4_genes.rds \
	refgenome/gmax_v4_transcripts.rds \
	refgenome/gmax_v4_exons.rds \
	refgenome/gmax_v4_cds.rds \
	gwas_results/vg/%_gwas.csv \
	gwas_results/vg/%_threshold_5per.txt
	$(RSCRIPT) figures/gene_subplot.R $* vg

figures/ggplots/paragraph_%_gene.rds: figures/gene_subplot.R \
	utilities/loci.txt \
	refgenome/Gmax_508_v4.0_mit_chlp.fasta \
	refgenome/lookup_v1_to_v4.rds \
	refgenome/gmax_v4_genes.rds \
	refgenome/gmax_v4_transcripts.rds \
	refgenome/gmax_v4_exons.rds \
	refgenome/gmax_v4_cds.rds \
	gwas_results/paragraph/%_gwas.csv \
	gwas_results/paragraph/%_threshold_5per.txt
	$(RSCRIPT) figures/gene_subplot.R $* paragraph

figures/ggplots/kmers_%_gene.rds: figures/gene_subplot.R \
	utilities/loci.txt \
	refgenome/Gmax_508_v4.0_mit_chlp.fasta \
	refgenome/lookup_v1_to_v4.rds \
	refgenome/gmax_v4_genes.rds \
	refgenome/gmax_v4_transcripts.rds \
	refgenome/gmax_v4_exons.rds \
	refgenome/gmax_v4_cds.rds \
	gwas_results/kmers/%_kmer_positions.rds \
	gwas_results/kmers/%_threshold_5per
	$(RSCRIPT) figures/gene_subplot.R $* kmers






