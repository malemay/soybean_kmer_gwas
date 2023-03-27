touch refgenome/glyma.Wm82.gnm1.FCtY.genome_main.fna
touch refgenome/Gmax_508_v4.0_mit_chlp.fasta
touch reference_signals/bandillo2017_signals_curated.tsv
touch refgenome/lookup_v1_to_v4.rds
touch refgenome/soybase_genome_annotation_v4.0_04-20-2021.txt
touch refgenome/Gmax_508_Wm82.a4.v1.gene_exons.gff3
touch reference_signals/bandillo2015_table1_curated.csv
touch refgenome/Gmax_508_v4.0_mit_chlp.fasta.fai
touch external_data/adapters.fa

mkdir -p illumina_data/raw_fastq

for i in $(cut -d , -f1 utilities/correct_sra_metadata.csv | tail -n +2 | sed 's/"//g')
do
	touch illumina_data/raw_fastq/${i}_1.fastq.gz
	touch illumina_data/raw_fastq/${i}_2.fastq.gz
done

