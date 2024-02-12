#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --time=3:00:00
#SBATCH --mem=100G
#SBATCH --job-name=devCNVs_BMI
#SBATCH --partition=general
#SBATCH --account=a_mcrae
#SBATCH -o /QRISdata/Q4399/Anorexia/UKB/burden_analysis/devCNVs/BMI/devCNVs_BMI.stdout
#SBATCH -e /QRISdata/Q4399/Anorexia/UKB/burden_analysis/devCNVs/BMI/devCNVs_BMI.stderr


module load r
export R_LIBS=/home/uqawal15/R_libraries/rlib_4.2.1

cd /QRISdata/Q4399/Anorexia/UKB/burden_analysis/devCNVs/BMI

R --file=/QRISdata/Q4399/Anorexia/UKB/burden_analysis/devCNVs/BMI/S4_devCNVs_burden_BMI.R