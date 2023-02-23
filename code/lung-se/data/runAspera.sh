#!/bin/bash
#
#SBATCH --job-name=get_ena
#SBATCH --time=48:00:00
#SBATCH --mem=32g
#SBATCH --cpus-per-task=4

module load aspera

tmpdir=.

key=/stornext/System/data/tools/aspera/aspera-3.9.6/etc/asperaweb_id_dsa.openssh

prj=../../../data/lung-se/misc/filereport_read_run_PRJNA341465_tsv.txt
dest=../../../data/lung-se/fastq/

path1=( $(awk -v FS='\t' -v OFS='\t' '{ print $31 }' $prj | tail -n +2 | awk -v FS=';' -v OFS=';' '{print $1}') )
# path2=( $(awk -v FS='\t' -v OFS='\t' '{ print $31 }' $prj | tail -n +2 | awk -v FS=';' -v OFS=';' '{print $2}') )

alias=( $(awk -v FS='\t' -v OFS='\t' '{ print $26 }' $prj | tail -n +2) )
len=${#alias[@]}

for i in $(seq 1 $len) ; do
  declare -i j=$i-1

  from1=${path1[$j]}
  # from2=${path2[$j]}

  aliasgz1=${alias[$j]}_R1.fastq.gz
  # aliasgz2=${alias[$j]}_R2.fastq.gz

  to1=${dest}${aliasgz1}
  # to2=${dest}${aliasgz2}

  echo "Downloading $from1..."
  ascp -QT -l 300m -P33001 -i $key era-fasp@$from1 $to1

  # echo "Downloading $from2..."
  # ascp -QT -l 300m -P33001 -i $key era-fasp@$from2 $to2

done

