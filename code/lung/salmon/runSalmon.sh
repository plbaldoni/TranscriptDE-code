#!/bin/bash

#SBATCH --mem=120g
#SBATCH --cpus-per-task=18
#SBATCH --time=6:00:00

homedir=/stornext/General/data/academic/lab_smyth/baldoni.p
salmon=$homedir/software/salmon-1.9.0_linux_x86_64/bin/salmon
salmonindex=../../../data/lung/index/transcripts_index/

for filename in ../../../data/lung/fastq/*_R1.fastq.gz; do

outname=${filename%_R1.fastq.gz}
outname=${outname##*/}

$salmon quant \
-i $salmonindex \
-l A \
-1 $filename \
-2 ${filename/_R1.fastq.gz/_R2.fastq.gz} \
-p 16 \
--numBootstraps 100 \
--validateMappings \
-o ../../../output/lung/salmon/$outname

done
