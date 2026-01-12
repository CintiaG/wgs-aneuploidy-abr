#!/bin/bash
#SBATCH --job-name=bedtools
#SBATCH --nodes=1
#SBATCH --time=0-8:00
#SBATCH --cpus-per-task=32
#SBATCH --mem=10GB
#SBATCH --account=mps_chrom_cp_nb
#SBATCH --output=jobs/bedtools.out.%j
#SBATCH --error=jobs/bedtools.err.%j
#SBATCH --mail-user=cintia.gomez_munoz@sorbonne-universite.fr
#SBATCH --mail-type=END

module load bedtools/2.29.2
module load parallel/20190322

INPUT_DIR=$1
OUTPUT_DIR=$2

echo -e "\nSTART\n"

date

echo "Running bedtools"

ls $INPUT_DIR/*.bam | sed 's/.*\/// ; s/_.*//' | \
parallel -j 32 "echo '{}' ; \
bedtools genomecov -d -ibam $INPUT_DIR/{}_sorted_marked.bam > $OUTPUT_DIR/{}_depth_per_base.txt ; \
bash scripts/doawk.sh $OUTPUT_DIR/{}_depth_per_base.txt $OUTPUT_DIR/{}_depth_per_base.bed ; \
bedtools map -b $OUTPUT_DIR/{}_depth_per_base.bed -a poly_CCN_and_VC/resources/ref/S288C_reference_sequence_R64-4-1_20230830_30Kb.bed -c 4 -o median > $OUTPUT_DIR/{}_depth_per_30Kb.bed ; \
bedtools genomecov -ibam $INPUT_DIR/{}_sorted_marked.bam > $OUTPUT_DIR/{}_histograms.bed"

echo "Done"

echo -e "\nEND\n"

date
