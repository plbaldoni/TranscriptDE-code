#!/bin/bash

#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G

module load picard-tools

annodir=../../../data/annotation/
destdir=../../../data/lung/index/

mkdir $destdir

# Adding sequin decoy chromosome to GRCh38.p13

NormalizeFasta LINE_LENGTH=60 I=${annodir}sequins/rnasequin_decoychr_2.4.fa O=${destdir}normalized_rnasequin_decoychr_2.4.fa

gzip ${destdir}normalized_rnasequin_decoychr_2.4.fa

cat ${annodir}/hg38/GRCh38.p13.genome.fa.gz ${destdir}normalized_rnasequin_decoychr_2.4.fa.gz > ${destdir}GRCh38.p13.genome_sequins.fa.gz

rm -rf ${destdir}normalized_rnasequin_decoychr_2.4.fa.gz

# Combining Gencode GTF with sequin GTF

gzip -c ${annodir}sequins/rnasequin_annotation_2.4.gtf > ${destdir}rnasequin_annotation_2.4.gtf.gz

cat ${annodir}/hg38/gencode.v33.annotation.gtf.gz ${destdir}rnasequin_annotation_2.4.gtf.gz > ${destdir}gencode.v33.annotation_sequins.gtf.gz

rm -rf ${destdir}rnasequin_annotation_2.4.gtf.gz

# Combining Gencode transcriptome with sequin transcriptome

gzip -c ${annodir}sequins/rnasequin_sequences_2.4.fa > ${destdir}rnasequin_sequences_2.4.fa.gz

cat ${annodir}/hg38/gencode.v33.transcripts.fa.gz ${destdir}rnasequin_sequences_2.4.fa.gz > ${destdir}gencode.v33.transcripts_sequins.fa.gz

gunzip ${destdir}gencode.v33.transcripts_sequins.fa.gz

NormalizeFasta LINE_LENGTH=60 I=${destdir}gencode.v33.transcripts_sequins.fa O=${destdir}out_gencode.v33.transcripts_sequins.fa

rm -rf ${destdir}gencode.v33.transcripts_sequins.fa

mv ${destdir}out_gencode.v33.transcripts_sequins.fa ${destdir}gencode.v33.transcripts_sequins.fa

gzip ${destdir}gencode.v33.transcripts_sequins.fa

# Removing

rm -rf ${destdir}rnasequin_sequences_2.4.fa.gz


