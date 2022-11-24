# This directory contains scripts used to run the complete simulation pipeline,
# including (1) sequencing files generation in FASTQ format, (2) quantification
# with Salmon and kallisto, and (3) differential transcript expression.
#
# Steps (1-3) are run on a per-experiment basis, and each experiment is run
# in a SLURM batch job triggered with array batch jobs.
# To run the simulation pipeline, do:
# 1 - run the file parameters.R with, for example, the command
#     Rscript parameters.R.
# 2 - run the SLURM bash script run.sh with, for example, sbatch run.sh on your
#     HPC running SLURM scheduler (change this accordingly).
#
# Step 1 above with create a parameter.txt file with all combinations of
# scenarios for the simulation. This file parameter.txt is read in Step 2 when
# you run run.sh. It provides the necessary parameters for each simulation.
#
# 3 - Once all simulations have been completed, run the script summarize.R with,
#     for example, the command Rscript summarize.R to summarize all the results.
#
# Computing time per step:
# 1 - approximately a minute
# 2 - approximately six days*
# 3 - approximately six hours*
#
# *: do note that these steps use parallel computing in several parts of the
#     scripts, please adjust accordingly.
