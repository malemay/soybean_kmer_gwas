touch refgenome/glyma.Wm82.gnm1.FCtY.genome_main.fna
touch refgenome/Gmax_508_v4.0_mit_chlp.fasta
touch reference_signals/bandillo2017_signals_curated.tsv
touch refgenome/lookup_v1_to_v4.rds
touch refgenome/soybase_genome_annotation_v4.0_04-20-2021.txt
touch refgenome/Gmax_508_Wm82.a4.v1.gene_exons.gff3
touch reference_signals/bandillo2015_table1_curated.csv
touch refgenome/Gmax_508_v4.0_mit_chlp.fasta.fai
touch external_data/adapters.fa
touch refgenome/Gmax4_db.bck
touch refgenome/Gmax4_db.des
touch refgenome/Gmax4_db.prj
touch refgenome/Gmax4_db.sds
touch refgenome/Gmax4_db.ssp
touch refgenome/Gmax4_db.suf
touch refgenome/Gmax4_db.tis

mkdir -p illumina_data/raw_fastq

for i in $(cut -d , -f1 utilities/correct_sra_metadata.csv | tail -n +2 | sed 's/"//g')
do
	touch illumina_data/raw_fastq/${i}_1.fastq.gz
	touch illumina_data/raw_fastq/${i}_2.fastq.gz
done

for i in $(cat external_data/nanopore_svs/files.txt)
do
	touch external_data/nanopore_svs/${i}
done

mkdir -p external_data/genome_assemblies/
for i in $(cut -d " " -f1 variant_calling/assemblies/assembly_samples.txt)
do
	touch external_data/genome_assemblies/${i}
	touch external_data/genome_assemblies/${i}.fai
done

for program in platypus kmers paragraph
do
	for i in $(cut -d, -f1 utilities/signal_ids.txt)
	do
		#touch gwas_results/${program}/${i}_locus.rds
		touch gwas_results/${program}/${i}_top_markers.rds
		touch gwas_results/${program}/${i}_all_genes.tsv
		touch gwas_results/${program}/${i}_top_genes.tsv
		touch gwas_results/${program}/${i}_nearest_gene.tsv
	done
done

mkdir -p illumina_data/merged_bams/
mkdir -p sv_genotyping/paragraph/manifest_files/

for i in $(cut -d ' ' -f1 utilities/srr_id_correspondence.txt)
do
	touch illumina_data/merged_bams/${i}_merged.bam
	touch illumina_data/merged_bams/${i}_merged.bam.bai
	touch sv_genotyping/paragraph/manifest_files/${i}_manifest.txt
done

touch external_data/BAC77G7-a.fasta
touch external_data/BAC77G7-a.fasta.fai
touch refgenome/Gmax_nuclv2_mit_chlp.fasta
touch external_data/soysnp50k_wm82.a2_41317.vcf.gz
