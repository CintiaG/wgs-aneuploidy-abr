#!/usr/bin/env Rscript
# Convert AD/DP tables (text) to FST
# Usage:
#   Rscript scripts/convert2fst.R out_data/sequencing_strains_AD.txt
library(fst)

args <- commandArgs(trailingOnly = TRUE)

FileIn <- args[1]

FileOut <- sub(".txt", ".fst", FileIn)

Df <- read.table(FileIn, header = TRUE)

write_fst(Df, FileOut)
