# This script carries out an association analysis for a single
# variant using GAPIT: a CNV at the Ps locus genotyped by smoove

# Loading the required packages
library(VariantAnnotation)
library(GAPIT3)
library(gwask)

# Loading the variant from the VCF file
smoove_vcf <- VcfFile("~/sv_gwas/usda_lines/sv_calling/smoove/all_samples.smoove.square.vcf.gz")
ps_cnv <- readVcf(smoove_vcf, param = GRanges(seqnames = "Gm12", ranges = IRanges(start = 36.2 * 10^6, end = 36.3 * 10^6)))
ps_cnv <- ps_cnv[info(ps_cnv)$SVTYPE != "BND"]
ps_cnv <- ps_cnv[abs(unlist(info(ps_cnv)$SVLEN)) < 100000 & abs(unlist(info(ps_cnv)$SVLEN)) > 10000]

# Making sure that we end up with a single variant
stopifnot(length(ps_cnv) == 1)

# Artificially doubling because GAPIT bugs otherwise
ps_cnv <- ps_cnv[c(1, 1)]

# Reading the kinship and pca results
kinship <- readRDS("filtered_variants/platypus_gapit_kinship.rds")
pca <- readRDS("filtered_variants/platypus_gapit_pca.rds")

gapit_results <- gapit_vcf(vcf = ps_cnv,
			   kinship = kinship,
			   pca = pca,
			   phenodata = "phenotypic_data/phenotypic_data.csv",
			   trait = "pubescence_density",
			   N.sig = min(length(ps_cnv), 20),
			   id_column = "bayer_id",
			   tmproot = "tmpdir")

