
#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --time=3:00:00
#SBATCH --mem=100G
#SBATCH --job-name=drCNVs
#SBATCH --partition=general
#SBATCH --account=a_mcrae
#SBATCH -o /Anorexia/UKB/burden_analysis/drCNVs/burden.stdout
#SBATCH -e /Anorexia/UKB/burden_analysis/drCNVs/burden.stderr

# Submit this job script to run S3_CNV_burden.R

module load r
export R_LIBS=/rlib_4.2.1

R --file=/Anorexia/UKB/burden_analysis/drCNVs/S3_CNV_burden.R