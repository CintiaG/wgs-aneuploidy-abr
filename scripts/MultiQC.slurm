#!/bin/bash
#SBATCH --job-name=multiqc
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=10GB
#SBATCH --time=0-1:00
#SBATCH --account=mps_chrom_cp_nb
#SBATCH --output=jobs/multiqc.out.%j
#SBATCH --error=jobs/multiqc.err.%j

module load multiqc/1.13

echo -e "START\n"

date

INPUT=$1

echo "Script ran with the following parameters:"

echo "multiqc $INPUT"

multiqc $INPUT

echo -e "END\n"

date
