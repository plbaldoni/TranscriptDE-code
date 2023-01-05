#!/bin/bash

#SBATCH --mem=120g
#SBATCH --cpus-per-task=18
#SBATCH --time=6:00:00

homedir=/stornext/General/data/academic/lab_smyth/baldoni.p

for filename in ../../../../data/mouse/single-end/fastq/*.fastq.gz; do

outname=${filename%.fastq.gz}
outname=${outname##*/}

$homedir/software/salmon-1.9.0_linux_x86_64/bin/salmon quant \
-i $homedir/software/SalmonIndex/mm39/transcripts_index/ \
-l A \
-r $filename \
-p 16 \
--numBootstraps 100 \
--validateMappings \
-o ../../../../output/mouse/single-end/salmon/$outname

done
