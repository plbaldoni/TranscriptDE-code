#!/bin/bash

#SBATCH --mem=120g 
#SBATCH --cpus-per-task=18
#SBATCH --time=6:00:00 

homedir=/stornext/General/data/academic/lab_smyth/baldoni.p
r1tail=_R1.fastq.gz
r2tail=_R2.fastq.gz

for filename in ../../../data/mouse/fastq/*$r1tail; do

outname=${filename%$r1tail}
outname=${outname##*/}

$homedir/software/salmon-1.9.0_linux_x86_64/bin/salmon quant \
-i $homedir/software/SalmonIndex/mm39/transcripts_index/ \
-l A \
-1 $filename \
-2 ${filename/$r1tail/$r2tail} \
-p 16 \
--numBootstraps 100 \
--validateMappings \
-o ../../../output/mouse/salmon/$outname

done
