# Plot sequencing depth histograms and estimate modal genome-wide depth
#
# This script reads per-sample depth histograms (BED-like table) and generates
# one depth-frequency plot per sample, colored by chromosome. The modal depth
# is estimated from the "genome" histogram (excluding Depth <= 10 to reduce noise).
#
# Inputs:
#   1) Histogram file (e.g., out_data/sequencing_histograms.bed)
#      Columns: Sample, Chr, Depth, Freq, Tot, Fract
#   2) Chromosome index file (e.g., files/chrom_index_cen.txt)
#      Must include: Chr_NC, Chr_rom, Colors
#
# Output:
#   - figures/histograms_<Prefix>/*_histogram.tiff   (one plot per sample)
#   - out_data/<Prefix>_depth_mode.csv              (Sample, Mode)
#
# Usage:
#   Rscript scripts/plot_depth.R \
#     out_data/sequencing_histograms.bed \
#     files/chrom_index_cen.txt \
#     sequencing \
#     600

# Load libraries
library(ggplot2)

# Read input arguments
args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 4) {
  stop("Usage: Rscript scripts/plot_depth.R <hist_file> <index_file> <prefix> <max_depth>")
}

HistFile <- args[1]
IndexFile <- args[2]
Prefix <- args[3]
Max <- as.numeric(args[4])
OutDir <- paste0("figures/histograms_", Prefix, "/")

dir.create(OutDir, recursive = TRUE, showWarnings = FALSE)
dir.create("out_data", recursive = TRUE, showWarnings = FALSE)

DfIndex <- read.table(IndexFile, header = TRUE, sep = "\t", quote = "\"",
                      stringsAsFactors = FALSE)
DfIndex <- DfIndex[DfIndex$Chr != "mitochondrial",]

# Build named palette keyed on histogram chromosome IDs (Chr_NC), plus genome
Pal <- c("darkgray", DfIndex$Colors)
names(Pal) <- c("genome", DfIndex$Chr_NC)

# Read all reads histogram data frame
HistDf <- read.table(HistFile, header = FALSE, sep = "\t",
                     col.names = c("Sample", "Chr", "Depth", "Freq", "Tot", "Fract"),
                     stringsAsFactors = FALSE)
HistDf$Depth <- as.numeric(HistDf$Depth)
HistDf$Freq  <- as.numeric(HistDf$Freq)

# Filter mitochondrial DNA
HistDf <- HistDf[HistDf$Chr != "ref|NC_001224|",]

Samples <- unique(HistDf$Sample)

Labeler <- setNames(DfIndex$Chr_rom, DfIndex$Chr_NC)
Labeler["genome"] <- "genome"

DfMode <- NULL

for (i in 1:length(Samples)){
  Sample <- Samples[i]
  #Code <- Codes[i]
  
  # Subset one sample
  WkHistDf <- HistDf[HistDf$Sample == Sample,]
  
  # Filter depth above 10 and below Max
  WkHistDf <- WkHistDf[WkHistDf$Depth > 10,]
  WkHistDf <- WkHistDf[WkHistDf$Depth < Max,]
  if (nrow(WkHistDf) == 0) next
  
  # Find genome reads mode (exclude the first 10 positions to filter noise) 
  WkHistDfGenome <- WkHistDf[WkHistDf$Chr == "genome", ]
  if (nrow(WkHistDfGenome) == 0) next
  
  ReadsMode <- WkHistDfGenome$Depth[which.max(WkHistDfGenome$Freq)]
  
  # Keep track of depth bin
  DfMode <- rbind(DfMode, data.frame("Sample" = Sample, "Mode" = ReadsMode))
  
  Pl <- ggplot(WkHistDf, aes(Depth, Freq, color = Chr)) +
    geom_line() +
    xlim(c(0, Max)) +
    geom_vline(xintercept = ReadsMode, linetype = "dashed", color = "red") +
    theme_bw() +
    scale_y_continuous(labels = scales::label_number()) +
    ggtitle(paste0(i, " ", Sample, " (Mode ", ReadsMode, ")")) +
    theme(plot.title = element_text(hjust = 0.5)) +
    scale_color_manual(values = Pal, breaks = names(Labeler), labels = Labeler)
  
  tiff(file = paste0(OutDir, i, "_", Sample, "_histogram.tiff"), compression = "lzw", width = 200, height = 150, units = "mm", res = 300)
  print(Pl)
  dev.off()
}

# Save
write.csv(DfMode, paste0("out_data/", Prefix, "_depth_mode.csv"), row.names = FALSE)
