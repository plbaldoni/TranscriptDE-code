#!/bin/bash

#SBATCH --mem=200g
#SBATCH --cpus-per-task=20
#SBATCH --time=48:00:00

# Creating decoy-aware transcriptome as descripted in the first bullet point of
# https://salmon.readthedocs.io/en/latest/salmon.html#preparing-transcriptome-indices-mapping-based-mode

destdir=../../../data/lung/index/
salmon=../../../../../../software/salmon-1.9.0_linux_x86_64/bin/salmon
salmondecoy=../../../../../../software/SalmonTools/scripts/generateDecoyTranscriptome.sh
mashmap=../../../../../../software/MashMap/mashmap
bedtools=/stornext/System/data/apps/bedtools/bedtools-2.26.0/bin/bedtools

cp ${destdir}GRCh38.p13.genome_sequins.fa.gz \
  ${destdir}gencode.v33.annotation_sequins.gtf.gz \
  ${destdir}gencode.v33.transcripts_sequins.fa.gz ./

gunzip GRCh38.p13.genome_sequins.fa.gz
gunzip gencode.v33.annotation_sequins.gtf.gz
gunzip gencode.v33.transcripts_sequins.fa.gz

$salmondecoy \
-j 18 \
-g GRCh38.p13.genome_sequins.fa \
-t gencode.v33.transcripts_sequins.fa \
-a gencode.v33.annotation_sequins.gtf \
-m $mashmap \
-b $bedtools \
-o $destdir

# Building index

$salmon index \
-t ${destdir}gentrome.fa \
-i ${destdir}transcripts_index \
-d ${destdir}decoys.txt \
-k 31 \
-p 18

# Removing extra files

rm -rf GRCh38.p13.genome_sequins.fa gencode.v33.transcripts_sequins.fa gencode.v33.annotation_sequins.gtf
