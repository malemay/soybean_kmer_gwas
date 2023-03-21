# Loading the required libraries
library(GenomicFeatures)

# Reading the gene annotations from a GFF file
# DEPENDENCY: refgenome/Gmax_508_Wm82.a4.v1.gene_exons.gff3
txdb <- makeTxDbFromGFF("refgenome/Gmax_508_Wm82.a4.v1.gene_exons.gff3")
gmax_v4_genes <- genes(txdb)

# Standardizing the gene names to lower case and saving to file
names(gmax_v4_genes) <- tolower(names(gmax_v4_genes))
saveRDS(gmax_v4_genes, file = "refgenome/gmax_v4_genes.rds")

# Extracting the transcripts by gene from the database and saving to file
gmax_v4_transcripts <- transcriptsBy(txdb, "gene")
names(gmax_v4_transcripts) <- tolower(names(gmax_v4_transcripts))
saveRDS(gmax_v4_transcripts, file = "refgenome/gmax_v4_transcripts.rds")

# Extracting the coding sequences by transcript
gmax_v4_cds <- cdsBy(txdb, "tx", use.names = TRUE)
saveRDS(gmax_v4_cds, file = "refgenome/gmax_v4_cds.rds")

# Extracting the exons by transcript
gmax_v4_exons <- exonsBy(txdb, "tx", use.names = TRUE)
saveRDS(gmax_v4_exons, file = "refgenome/gmax_v4_exons.rds")

