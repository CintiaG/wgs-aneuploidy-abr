#!/bin/bash
#SBATCH --job-name=fastqc
#SBATCH --nodes=1
#SBATCH --cpus-per-task=20
#SBATCH --mem=20GB
#SBATCH --time=0-1:00
#SBATCH --account=mps_chrom_cp_nb
#SBATCH --output=jobs/fastqc.out.%j
#SBATCH --error=jobs/fastqc.err.%j

module load fastqc/0.12.1

echo -e "START\n"

date

INPUT=$1

OUTPUT=$2

echo "Script ran with the following parameters:"

echo "fastqc -o $OUTPUT -t 20 --noextract "$INPUT"/*.fq.gz"

fastqc -o "$OUTPUT" -t 20 --noextract "$INPUT"/*.fq.gz

echo -e "END\n"

date
