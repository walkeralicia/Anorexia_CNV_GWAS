#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --time=3:00:00
#SBATCH --mem=100G
#SBATCH --job-name=devCNVs
#SBATCH --partition=general
#SBATCH --account=a_mcrae
#SBATCH -o /QRISdata/Q4399/Anorexia/UKB/burden_analysis/devCNVs/devCNVs.stdout
#SBATCH -e /QRISdata/Q4399/Anorexia/UKB/burden_analysis/devCNVs/devCNVs.stderr


module load r
export R_LIBS=/home/uqawal15/R_libraries/rlib_4.2.1

cd /QRISdata/Q4399/Anorexia/UKB/burden_analysis/devCNVs

R --file=/QRISdata/Q4399/Anorexia/UKB/burden_analysis/devCNVs/S3_devCNVs_burden.R