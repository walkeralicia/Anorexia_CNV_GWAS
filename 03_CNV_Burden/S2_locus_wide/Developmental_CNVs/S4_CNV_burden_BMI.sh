#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --time=3:00:00
#SBATCH --mem=100G
#SBATCH --job-name=devCNVs_BMI
#SBATCH --partition=general
#SBATCH --account=a_mcrae
#SBATCH -o /Anorexia/UKB/burden_analysis/devCNVs/BMI/devCNVs_BMI.stdout
#SBATCH -e /Anorexia/UKB/burden_analysis/devCNVs/BMI/devCNVs_BMI.stderr


module load r
export R_LIBS=/rlib_4.2.1

cd /Anorexia/UKB/burden_analysis/devCNVs/BMI

R --file=/Anorexia/UKB/burden_analysis/devCNVs/BMI/S4_devCNVs_burden_BMI.R