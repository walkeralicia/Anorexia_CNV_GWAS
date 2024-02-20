#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --time=3:00:00
#SBATCH --mem=100G
#SBATCH --job-name=devCNVs
#SBATCH --partition=general
#SBATCH --account=a_mcrae
#SBATCH -o /Anorexia/UKB/burden_analysis/devCNVs/devCNVs.stdout
#SBATCH -e /Anorexia/UKB/burden_analysis/devCNVs/devCNVs.stderr

# This bash script submits the S3_CNV_burden.R script as a job

module load r
export R_LIBS=/rlib_4.2.1

R --file=/Anorexia/UKB/burden_analysis/devCNVs/S3_devCNVs_burden.R