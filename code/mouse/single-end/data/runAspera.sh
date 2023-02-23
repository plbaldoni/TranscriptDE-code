#!/bin/bash
#
#SBATCH --job-name=get_ena
#SBATCH --time=48:00:00
#SBATCH --mem=32g
#SBATCH --cpus-per-task=4

module load aspera

tmpdir=.

key=/stornext/System/data/tools/aspera/aspera-3.9.6/etc/asperaweb_id_dsa.openssh

prj=../../../../data/mouse/single-end/misc/filereport_read_run_PRJNA258286_tsv.txt
dest=../../../../data/mouse/single-end/fastq/
path=( $(awk -v FS='\t' -v OFS='\t' '{ print $31 }' $prj | tail -n +2 | awk -v FS=';' -v OFS=';' '{print $1}') )
alias=( $(awk -v FS='\t' -v OFS='\t' '{ print $26 }' $prj | tail -n +2) )
len=${#alias[@]}

for i in $(seq 1 $len) ; do
  declare -i j=$i-1
  from=${path[$j]}

  aliasgz=${alias[$j]}.fastq.gz

  to=${dest}${aliasgz}

  echo "Downloading $from..."
  ascp -QT -l 300m -P33001 -i $key era-fasp@$from $to
done

