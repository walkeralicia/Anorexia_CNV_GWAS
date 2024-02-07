
#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --time=3:00:00
#SBATCH --mem=100G
#SBATCH --job-name=drCNVs
#SBATCH --partition=general
#SBATCH --account=a_mcrae
#SBATCH -o /QRISdata/Q4399/Anorexia/UKB/burden_analysis/drCNVs/burden.stdout
#SBATCH -e /QRISdata/Q4399/Anorexia/UKB/burden_analysis/drCNVs/burden.stderr



module load r
export R_LIBS=/home/uqawal15/R_libraries/rlib_4.2.1

cd /QRISdata/Q4399/Anorexia/UKB/burden_analysis/drCNVs

R --file=/QRISdata/Q4399/Anorexia/UKB/burden_analysis/drCNVs/S3_CNV_burden.R