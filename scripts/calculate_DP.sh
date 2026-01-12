#!/bin/bash
# A script to calculate DP
# run: screen -SL bcf_tools -Logfile logs/calculate_DP.log bash scripts/calculate_DP.sh results/vcf/sequencing_strains.vcf.gz out_data/

INPUT=$1
OUT_DIR=$2
PREFIX=$(basename $INPUT)
OUTPUT=$(echo "$PREFIX" | sed 's/\..*//')_DP.txt

echo "bcftools query -l $INPUT | sed ':a;N;$!ba;s/\n/\t/g; s/^/CHROM\tPOS\tREF\tALT\t/' > $OUT_DIR/$OUTPUT"
bcftools query -l $INPUT | sed ':a;N;$!ba;s/\n/\t/g; s/^/CHROM\tPOS\tREF\tALT\t/' > $OUT_DIR/$OUTPUT

echo "bcftools query -f '%CHROM\t%POS\t%REF\t%ALT[\t%DP]\n' $INPUT >> $OUT_DIR/$OUTPUT"
bcftools query -f '%CHROM\t%POS\t%REF\t%ALT[\t%DP]\n' $INPUT >> $OUT_DIR/$OUTPUT
