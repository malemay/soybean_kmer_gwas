# Code for the analysis of *k*-mer- and structural variant-based GWAS in soybean

## Overview

This repository contains all the code needed to reproduce the analyses
presented in the paper titled "*k*-mer-based GWAS enhances the discovery of
causal variants and candidate genes in soybean".

As a disclaimer, readers should be aware that most of the code was reorganized
and integrated into the Makefile only after analyses were performed.
Therefore, those trying to run the analyses might run into issues related to
paths or software version incompatibilities.  We encourage users who encounter
problems while trying to run this code to open a GitHub issue or contact the
repo maintainer directly. We believe that the code in this repository and the
associated Makefile should still be useful to help those interested in
understanding the analyses that were performed.

## Software dependencies

The following software should be installed to reproduce the analyses.  Some of
these programs may themselves have additional dependencies.  It is assumed that
all these programs are found in your `$PATH` for the analyses to run properly.

* [AsmVar](https://github.com/bioinformatics-centre/AsmVar)
* [BBDuk](https://jgi.doe.gov/data-and-tools/bbtools/bb-tools-user-guide/bbduk-guide/)
* [bamaddrg](https://github.com/ekg/bamaddrg)
* [BayesTyperTools](https://github.com/bioinformatics-centre/BayesTyper)
* [bwa](https://github.com/lh3/bwa)
* [bcftools](https://github.com/samtools/bcftools)
* [gwask](https://github.com/malemay/gwask) R package
* [htslib](https://github.com/samtools/htslib)
* [katcher](https://github.com/malemay/katcher) and associated binaries
* [KMC](https://github.com/refresh-bio/KMC)
* [*k*-mer GWAS](https://github.com/voichek/kmersGWAS) binaries
* [kmers_ld](https://github.com/malemay/kmers_ld)
* [LAST](https://gitlab.com/mcfrith/last)
* [MAFFT](https://mafft.cbrc.jp/alignment/software/)
* [manta](https://github.com/Illumina/manta)
* [MUMmer](https://github.com/mummer4/mummer)
* [paragraph](https://github.com/Illumina/paragraph)
* [Platypus](https://www.well.ox.ac.uk/research/research-groups/lunter-group/lunter-group/platypus-a-haplotype-based-variant-caller-for-next-generation-sequence-data)
* [PLINK](https://www.cog-genomics.org/plink2)
* [R programming language](https://cran.r-project.org/)
* [samtools](https://github.com/samtools/samtools)
* [smoove](https://github.com/brentp/smoove)
* [SOAPdenovo2](https://github.com/aquaskyline/SOAPdenovo2)
* [SPAdes](https://github.com/ablab/spades)
* [SvABA](https://github.com/walaj/svaba)
* [SVanalyzer](https://github.com/nhansen/SVanalyzer)
* [svmutools](https://github.com/malemay/svmutools) R package
* [TASSEL command-line tools](https://bitbucket.org/tasseladmin/tassel-5-source/wiki/Home)
* [TeX Live](https://www.tug.org/texlive/)
* [vcftools](https://github.com/vcftools/vcftools)

Some programs needed for reproducing analyses were modified from existing software:

* We forked [svmu](https://github.com/mahulchak/svmu) and slightly modified it
  to improve memory usage and execution speed. This version can be installed
  from [our fork](https://github.com/malemay/svmu) by using the commit 378719b on
  branch malemay-fork.

* The script
  [scripts/addMissingPaddingGmax4.py](https://github.com/malemay/soybean_sv_paper/blob/master/scripts/addMissingPaddingGmax4.py)
was adapted from
[addMissingPaddingHg38.py](https://github.com/vgteam/sv-genotyping-paper/blob/master/human/misc-scripts/addMissingPaddingHg38.py)
to use the soybean reference genome instead of the human reference genome.

## Data availability

### Sequencing data

* The Illumina data used in this project is available from the
  [SRA](https://www.ncbi.nlm.nih.gov/sra) using the BioProject accession
  numbers [PRJNA257011](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA257011),
  [PRJNA289660](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA289660),
  and [PRJNA639876](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA639876).
  This data should be placed under `illumina_data/raw_fastq/` in compressed
  fastq format to reproduce the analyses.

### High-quality assemblies for SV calling

* The assemblies of [Liu et al. (2020)](https://doi.org/10.1016/j.cell.2020.05.023)
  are available on the Genome Warehouse through Accession Number
  [PRJCA002030](https: //ngdc.cncb.ac.cn/search/?dbId=gsa&q=PRJCA002030).

* The assemblies of ZH13, W05 and Lee are available on
  [SoyBase](https://soybase.org/ GlycineBlastPages/blast_descriptions.php).

All these assemblies should be placed in external_data/genome_assemblies/

### Structural variants called from Oxford Nanopore data

The SVs identified from Oxford Nanopore data by [Lemay et al. (2022)](https://doi.org/10.1186/s12915-022-01255-w)
are [available on figshare](https://doi.org/10.6084/ m9.figshare.15127730.v1).
These should be placed under `external_data/nanopore_svs/`.

### SoySNP50K calls

SoySNP50K genotype calls are available from [SoyBase](https://soybase.org/snps/).
These should be placed under external_data/soysnp50k_wm82.a2_41317.vcf.gz

### Reference data

The following datasets are available from the Web and should be added to the
repository to reproduce the analyses:

* The reference genome sequence and annotation of soybean cultivar Williams82,
  assembly version 4 can be downloaded from
  [Phytozome](https://phytozome-next.jgi.doe.gov/).  The files needed
  (`Gmax_508_v4.0.fa`, `Gmax_508_Wm82.a4.v1.gene_exons.gff3`) should be placed
  under `refgenome/`.

* The soybean chloroplast and mitochondrion genome sequences can be downloaded
  from
[SoyBase](https://www.soybase.org/GlycineBlastPages/blast_descriptions.php).
They should also be concatenated to `Gmax_508_v4.0.fa` and placed under the
name `refgenome/Gmax_508_v4.0_mit_chlp.fasta`

### Data generated by the analysis

Some of the intermediate datasets generated as part of this analysis are
 [available on figshare](https://doi.org/10.6084/m9.figshare.21699464.v3).

### Other datasets

* Reference Illumina adapters are distributed with BBDUK and should be placed
  under external_data/adapters.fa.

* The *I* locus contig assembled by [Tuteja and Vodkin (2008)](https://doi.org/10.2135/cropsci2007.10.0542tpg)
  is available from [NCBI](https://www.ncbi.nlm.nih.gov/nuccore/EF623854)
  and should be placed under `external_data/BAC77G7-a.fasta`.

* The signals discovered by [Bandillo et al. (2017)](https://doi.org/10.3835/plantgenome2016.06.0054)
	and [Bandillo et al. (2015)](https://doi.org/10.3835/plantgenome2015.04.0024) are available from their publications
	and were placed under `reference_signals/bandillo2017_signals_curated.tsv` and
	`reference_signals/bandillo2015_table1_curated.csv`, respectively.

## Querying the Makefile

Please refer to the
[corresponding section](https://github.com/malemay/soybean_sv_paper#querying-the-makefile) in
our previous work for instructions on how to use the Makefile to understand the
sequence of the analyses performed as part of this work.

## Citation

The reference to this software will be updated shortly.

