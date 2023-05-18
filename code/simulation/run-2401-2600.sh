#!/bin/bash 
#SBATCH --job-name=tx-pipe
#SBATCH --array=401-600
#SBATCH --mem=66G
#SBATCH --cpus-per-task=11
#SBATCH --time=48:00:00
#SBATCH --qos=bonus
#SBATCH --output=slurm-%A_%a.out

workers=10
CurrentDir=$(pwd)
MaxJobNum=$(tail -n +2 parameters.txt | wc -l)

# Turn on below for job numbers greater than 1000 (there is an array limit in the HPC)
StepSize=2000
SLURM_ARRAY_TASK_ID=$((SLURM_ARRAY_TASK_ID+StepSize))

# If array job is greater than the maximum number of jobs, exit
if [[ $SLURM_ARRAY_TASK_ID -gt $MaxJobNum ]]; then
  exit 1
fi

# Reading parameters
ParVec=$(tail -n +2 parameters.txt | sed -n "$SLURM_ARRAY_TASK_ID"p)
dest=$(echo $ParVec | awk '{ print $1 }')
genome=$(echo $ParVec | awk '{ print $2 }')
rlen=$(echo $ParVec | awk '{ print $3 }')
fc=$(echo $ParVec | awk '{ print $4 }')
pe=$(echo $ParVec | awk '{ print $5 }')
mtx=$(echo $ParVec | awk '{ print $6 }')
scenario=$(echo $ParVec | awk '{ print $7 }')
libs=$(echo $ParVec | awk '{ print $8 }')
simulation=$(echo $ParVec | awk '{ print $9 }')

# Check if job has been completed
if [[ -f $dest/dte-salmon/time.tsv && -f $dest/dte-kallisto/time.tsv ]]; then
  exit 1
fi

# Submitting Rscript
echo $dest
cd $dest
Rscript --no-save --no-restore --verbose $CurrentDir/main.R \
  --dest=$dest \
  --genome=$genome \
  --rlen=$rlen \
  --fc=$fc \
  --pe=$pe \
  --mtx=$mtx \
  --scenario=$scenario \
  --libs=$libs \
  --simulation=$simulation \
  --workers=$workers > output.log 2>&1
