#!/bin/bash

# Setting the current directory
cd variant_calling/asmvar/

# This variable contains the path to the index database for LAST
# DEPENDENCY: the following files comprising the database
# refgenome/Gmax4_db.bck
# refgenome/Gmax4_db.des
# refgenome/Gmax4_db.prj
# refgenome/Gmax4_db.sds
# refgenome/Gmax4_db.ssp
# refgenome/Gmax4_db.suf
# refgenome/Gmax4_db.tis
refdb=../../../../refgenome/Gmax4_db

# Creating the directory structure and configuration files to run SOAPdenovo2 on all samples
mkdir -p soapdenovo_assembly
cd soapdenovo_assembly

# DEPENDENCY: utilities/srr_id_correspondence.txt
# DEPENDENCY: illumina_data/BBDUK_TRIMMING
# First creating a directory for all samples
cat ../../../utilities/srr_id_correspondence.txt | while read i
do
	bayer_id=$(echo $i | cut -d " " -f1)
	mkdir -p $bayer_id

	srr_ids=$(echo $i | cut -d " " -f2)
	IFS=';' read -r -a srr_array <<< "$srr_ids"

	# Outputting the maximum read length to the config file
	echo "max_rd_len=250" > ${bayer_id}/${bayer_id}_config.txt

	for id in ${srr_array[@]}
	do
		echo "
[LIB]
avg_ins=500
reverse_seq=0
asm_flags=1
q1=../../../../illumina_data/trimmed_fastq/${id}/${id}_R1_trimmed.fastq.gz
q2=../../../../illumina_data/trimmed_fastq/${id}/${id}_R1_trimmed.fastq.gz

[LIB]
asm_flags=1
q=../../../../illumina_data/trimmed_fastq/${id}/${id}_sing_trimmed.fastq.gz
" >> ${bayer_id}/${bayer_id}_config.txt
	done
done

# Moving back to the main AsmVar directory
cd ..

# Looping over all samples
for i in $(cut -d " " -f1 ../../../utilities/srr_id_correspondence.txt)
do

	# Going to the directory for the assembly of that sample
	cd soapdenovo_assembly/${i}

	# Performing the assembly with SOAPdenovo2
	SOAPdenovo-63mer sparse_pregraph -p 2 -s ${i}_config.txt -K 49 -R -z 1100000000 -o ${i} 1>pregraph.log 2>pregraph.err
	SOAPdenovo-63mer contig -p 2 -g ${i} -R 1>contig.log 2>contig.err

	if [ $? -eq 0 ]
	then
		rm ${i}.path ${i}.ht_content ${i}.sparse.edge ${i}.updated.edge ${i}.edge.gz ${i}.preArc ${i}.Arc ${i}.vertex ${i}.ht_idx
	fi

	# Aligning with LAST
	# -m 0.01 to keep the previous default setting of last-926
	lastal -D1000 -P2 -Q0 -e20 -j4 -v $refdb ${i}.contig | last-split -m 0.01 -s30 -v > ${i}.maf

	# Moving back to the main AsmVar directory directory
	cd ../..

  	# Creating output directory for AsmVar calls and moving to it
  	mkdir -p ${i}
  	cd ${i}

  # Creating some variables to make scripting clearer
  query_dir=../soapdenovo_assembly/${i}
  # DEPENDENCY: refgenome/Gmax_508_v4.0_mit_chlp.fasta
  ref=../../../refgenome/Gmax_508_v4.0_mit_chlp.fasta
  chromosomes=$(seq -w 1 20 | xargs -I XX echo GmXX)

  for j in $chromosomes
  do
    ASV_VariantDetector -s ${i} -i ${query_dir}/${i}.maf -t $ref -q ${query_dir}/${i}.contig -o asmvar_results_${j} -r $j > ${i}_${j}.age
  done

  # Moving back to the starting directory
  cd ..
done

touch ASMVAR_CALLING
