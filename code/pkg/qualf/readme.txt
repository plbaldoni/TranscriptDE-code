# Files here are the quality scores from two experiments with 50bp, 125bp, and  150bp
# reads. There are 20,000 in each file and these are used to simulate RNA-seq
# reads with Rsubread::simReads with the appropriate read length. These
# files were created with the following command calls:
# gunzip -c /vast/scratch/users/baldoni.p/tmp/ENCFF713MNU.fastq.gz | sed -n '0~4p' | sed -r '/^.{,49}$/d' | sed '/^.\{50\}./d' | shuf -n 20000 > ref-quality-strings-20k-50bp-ENCFF713MNU.txt
# gunzip -c /vast/scratch/users/baldoni.p/tmp/ENCFF102BXZ.fastq.gz | sed -n '0~4p' | sed -r '/^.{,150}$/d' | sed '/^.\{151\}./d' | sed 's/.$//' | shuf -n 20000 > ref-quality-strings-20k-150bp-ENCFF102BXZ.txt
# gunzip -c /vast/scratch/users/baldoni.p/ENCFF126GLV.fastq.gz | sed -n '0~4p' | sed -r '/^.{,124}$/d' | sed '/^.\{125\}./d' | shuf -n 20000 > ref-quality-strings-20k-125bp-ENCFF126GLV.txt
