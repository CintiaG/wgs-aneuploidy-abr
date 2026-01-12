# Compile sequencing, quality control, and alignment statistics
#
# This script combines Trimmomatic trimming statistics, post-trimming
# FastQC/MultiQC metrics, and alignment statistics into a single
# per-sample summary table.
#
# Samples are ordered according to the original trimming statistics
# to preserve the initial sequencing order.
#
# Usage:
#   Rscript scripts/get_sequencing_stats.R \
#     out_data/trimming_stats.txt \
#     quals/multiqc_rnd2/multiqc_fastqc.txt \
#     out_data/alignment_stats.txt \
#     out_data/sequencing_stats.txt

library(dplyr)

# Read input arguments
args <- commandArgs(trailingOnly = TRUE)

TrimFile <- args[1]
QualFile <- args[2]
AlnFile <- args[3]
OutFile <- args[4]

# Trimmomatic stats file
Df1 <- read.table(TrimFile, sep = "\t", header = TRUE)

# MultiQC file after trimming
Df2 <- read.table(QualFile, sep = "\t", header = TRUE)

# Alignment stats file
Df3 <- read.table(AlnFile, sep = "\t", header = TRUE)

# Separate Total Bases info in the MultiQC file
Df2$Total.Bases.Units <- sub(".* Gbp", "Gbp", sub(".* Mbp", "Mbp", Df2$Total.Bases))
Df2$Total.Bases <- as.numeric(sub(" Gbp", "", sub(" Mbp", "", Df2$Total.Bases)))
Df2$Total.Bases.Mbp <- ifelse(
  Df2$Total.Bases.Units == "Mbp",
  Df2$Total.Bases,
  Df2$Total.Bases * 1e3
)

# Rename sample
Df2$Sample <- sub("_.*", "", Df2$Sample)

# Summarize data frame by sample
Df2 <- Df2 %>% group_by(Sample) %>%
  summarize(
    GC_mean = mean(X.GC),
    Total_reads = first(Total.Sequences),
    Reads_length_mean = round(mean(avg_sequence_length), 1),
    Total_bases_Mbp = sum(Total.Bases.Mbp)
  )

# Join data frames information
Df4 <- left_join(Df1, Df2, by = "Sample")
Df4 <- left_join(Df4, Df3, by = "Sample")

# Save overall stats
write.table(Df4, OutFile, row.names = FALSE, quote = FALSE, sep = "\t")
