#!/bin/bash
#SBATCH --job-name=trimming
#SBATCH --nodes=1
#SBATCH --time=0-8:00
#SBATCH --cpus-per-task=32
#SBATCH --mem=20GB
#SBATCH --account=mps_chrom_cp_nb
#SBATCH --output=jobs/trimming.out.%j
#SBATCH --error=jobs/trimming.err.%j

module load trimmomatic/0.39
module load parallel/20190322

echo -e "\nSTART\n"

date

echo "Running Trimmomatic"

ls reads/raw/ | sed 's/_.*//' | sort -u | parallel -j 4 "echo '{}' ; trimmomatic PE -validatePairs -threads 8 reads/raw/{}_1.fastq.gz reads/raw/{}_2.fastq.gz reads/paired/{}_1.fastq.gz reads/unpaired/{}_unp_1.fastq.gz reads/paired/{}_2.fastq.gz reads/unpaired/{}_unp_2.fastq.gz ILLUMINACLIP:scripts/NexteraPE-PE.fa:2:30:10 HEADCROP:15 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36"

echo "Done"

echo -e "\nEND\n"

date
