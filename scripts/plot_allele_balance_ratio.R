#!/usr/bin/env Rscript
# Plot combined chromosome copy number profiles and allele balance ratios (ABR)
#
# This script generates per-sample multi-panel figures combining:
#   1) Genome-wide copy number histograms
#   2) Copy number profiles in 30 kb windows with detected aneuploid segments
#   3) Genome-wide allele balance ratio (ABR) distributions
#   4) ABR values plotted along the genome
#
# The script is robust to optional metadata columns and will automatically
# adapt plot titles and labels depending on availability of:
#   - Code
#   - Mean_ploidy
#   - Rounded_ploidy
#   - Ploidy_bioinfo_exact
#
# In addition, the script computes and exports the number of heterozygous sites
# per sample (after ABR filtering) to:
#   out_data/<prefix>_het_sites.csv
#
# Inputs:
#   1) Sample information CSV (e.g. files/sequencing.csv)
#   2) Allele balance ratio table in fst format (e.g. out_data/sequencing_ABR.fst)
#   3) Copy-number histogram table (30 kb windows)
#   4) Copy-number depth table (30 kb windows)
#   5) Chromosome index file with centromere positions and colors
#   6) Per-chromosome copy-number summary table
#   7) Detected aneuploid regions table
#   8) Output prefix (e.g. sequencing)
#   9) Maximum copy number to display in plots
#
# Outputs:
#   - Per-sample ABR and copy-number figures:
#       figures/allele_balance_ratio_<prefix>/*_ABR.tiff
#   - Per-sample heterozygous site counts:
#       out_data/<prefix>_het_sites.csv
#
# Usage:
#   Rscript scripts/plot_allele_balance_ratio.R \
#     files/sequencing.csv \
#     out_data/sequencing_ABR.fst \
#     out_data/sequencing_histograms_30Kb_aneu.csv \
#     out_data/sequencing_depth_per_30Kb_aneu.csv \
#     files/chrom_index_cen.txt \
#     out_data/sequencing_copy_number_30Kb.csv \
#     out_data/sequencing_aneuploid_regions.csv \
#     sequencing \
#     8
library(fst)
library(ggplot2)
library(ggpubr)
library(ggrepel)

# Read input arguments
args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 9) {
  stop("Usage: Rscript scripts/plot_allele_balance_ratio.R <strains.csv> <abr.fst> <hist.csv> <depth.csv> <chrom_index.txt> <copies.csv> <regions.csv> <prefix> <max_copies>")
}

StrainsFile <- args[1]
AlleleFile <- args[2]
HistFile <- args[3]
DepthFile <- args[4]
IndexFile <- args[5]
CopiesFile <- args[6]
RegionsFile <- args[7]
Prefix <- args[8]
OutDir <- paste0("figures/allele_balance_ratio_", Prefix, "/")

dir.create("out_data", recursive = TRUE, showWarnings = FALSE)
dir.create(OutDir, recursive = TRUE, showWarnings = FALSE)

Top <- as.numeric(args[9])

# Load index data frame to get mids and lines information
DfIndex <- read.table(IndexFile, header = TRUE, sep = "\t", quote = "\"",
                      stringsAsFactors = FALSE)
DfIndex <- DfIndex[!DfIndex$Chr == "mitochondrial",]

StrainsDf <- read.csv(StrainsFile)

AlleleDf <- read_fst(AlleleFile)

HistDf <- read.csv(HistFile)

DepthDf <- read.csv(DepthFile)

CopiesDf <- read.csv(CopiesFile)

RegionsDf <- read.csv(RegionsFile)

# Filter aberrant number of alleles
AlleleDf <- AlleleDf[AlleleDf$Num > 1,]
AlleleDf <- AlleleDf[AlleleDf$Num < 4,]

AlleleDf <- AlleleDf[AlleleDf$Ref_Alt == "ALT",]

# Remove samples that were not sequenced
StrainsDf <- StrainsDf[!is.na(StrainsDf$Sample),]
DepthDf$Sample <- gsub("-", ".", DepthDf$Sample)
HistDf$Sample <- gsub("-", ".", HistDf$Sample)

Samples <- gsub("-", ".", StrainsDf$Sample)
HasCode <- "Code" %in% colnames(StrainsDf)
Strains <- if (HasCode) StrainsDf$Code else StrainsDf$Sample
Ploidies <- StrainsDf$Base_ploidy

for (i in 1:length(Strains)){
  Strain <- Strains[i]
  Sample <- Samples[i]
  Ploidy <- Ploidies[i]
  
  # Use Code if available, otherwise fallback to Sample
  LabelID <- if (HasCode && !is.na(StrainsDf$Code[i]) && StrainsDf$Code[i] != "") StrainsDf$Code[i] else Samples[i]
  
  # Safe access to optional ploidy columns
  HasMeanPloidy <- "Mean_ploidy" %in% colnames(StrainsDf)
  HasRoundedPloidy <- "Rounded_ploidy" %in% colnames(StrainsDf)
  HasBioinfoPloidy <- "Ploidy_bioinfo_exact" %in% colnames(StrainsDf)
  
  MeanP <- if (HasMeanPloidy) StrainsDf$Mean_ploidy[i] else NA
  RoundP <- if (HasRoundedPloidy) StrainsDf$Rounded_ploidy[i] else NA
  BaseP <- StrainsDf$Base_ploidy[i]
  BioP  <- if (HasBioinfoPloidy) StrainsDf$Ploidy_bioinfo_exact[i] else NA
  
  UpThres <- (1 + round(0.7 * (1 / Ploidy), 2)) * Ploidy
  DownThres <- (1 - round(0.7 * (1 / Ploidy), 2)) * Ploidy
  
  # Working data frames
  WkAlleleDf <- AlleleDf[AlleleDf$Strain == Sample,]
  if (nrow(WkAlleleDf) < 1){
    WkAlleleDf <- AlleleDf[AlleleDf$Strain == Strain,]
  }
  
  WkHistDf <- HistDf[HistDf$Sample == Sample,]
  WkDepthDf <- DepthDf[DepthDf$Sample == Sample,]
  
  HasCodeCopies <- "Code" %in% colnames(CopiesDf)
  HasCodeRegions <- "Code" %in% colnames(RegionsDf)
  
  WkCopiesDf <- if (HasCodeCopies) CopiesDf[CopiesDf$Code == Strain,] else CopiesDf[CopiesDf$Sample == Sample,]
  WkRegionDf <- if (HasCodeRegions) RegionsDf[RegionsDf$Code == Strain,] else RegionsDf[RegionsDf$Sample == Sample,]
  
  if (nrow(WkCopiesDf) == 0) {
    warning("No copy-number rows found for: ", LabelID, " (skipping)")
    next
  }
  if (nrow(WkDepthDf) == 0 || nrow(WkHistDf) == 0) {
    warning("Missing depth/hist data for: ", LabelID, " (skipping)")
    next
  }
  if (nrow(WkAlleleDf) == 0) {
    warning("No ABR rows found for: ", LabelID, " (skipping)")
    next
  }
  
  ChrCols <- setNames(DfIndex$Colors, DfIndex$Chr_NC)
  
  PlA <- ggplot(WkHistDf, aes(Copies, Freq, color = Chr)) +
    geom_line(show.legend = FALSE, linewidth = 0.25) +
    geom_vline(xintercept = StrainsDf$Base_ploidy[i], linetype = "dashed", color = "red", linewidth = 0.25) +
    theme_bw() +
    scale_y_continuous(labels = function(x) x / 10000) +
    scale_x_continuous(limits = c(0, Top), labels = function(x){
      paste0("  ", x, ".0")
    }) +
    scale_color_manual(values = ChrCols) +
    coord_flip() +
    ylab("Freq (1e4)")
  
  PlB <- ggplot(WkDepthDf, aes(Index, CopiesRaw, color = Chr_new)) +
    geom_line(show.legend = FALSE, linewidth = 0.25) +
    scale_x_continuous(breaks = DfIndex$Mids_30000,
                       labels = as.roman(1:16),
                       expand = c(0, 0)) +
    geom_vline(xintercept = DfIndex$Lines_30000,
               linetype = "dashed", color = "gray30") +
    theme_bw() +
    scale_color_manual(values = WkCopiesDf$Colors) +
    scale_y_continuous(limits = c(0, Top), labels = function(x){
      paste0("  ", x, ".0")
    }) +
    xlab("Chromosome") +
    ylab("Copies") +
    geom_text(data = WkCopiesDf,
              aes(x = Mids,
                  y = 5,
                  label = Label),
              color = ifelse(WkCopiesDf$Mean_Copies == WkCopiesDf$Label,
                             "black",
                             "red"),
              show.legend = FALSE) +
  # Upper threshold
  geom_hline(yintercept = UpThres, linetype="dashed", 
             color = "red") +
  # Lower threshold
  geom_hline(yintercept = DownThres, linetype="dashed", 
             color = "red")
  
  # Make sure that regions were obtained
  if (nrow(WkRegionDf) > 0){
    # Add aneuploidies
    PlB <- PlB + geom_rect(data = WkRegionDf, aes(xmin = Start_index - 1, xmax = End_index + 1, ymin = Ymin - 0.01, ymax = Ymax + 0.01), alpha = 0.3, inherit.aes = FALSE, color = "red") +
      geom_text_repel(show.legend = FALSE, data = WkRegionDf, aes(x = Mid_index, y = Mid_y, label = Aneu), color = "black", nudge_x = 10, nudge_y = -1)
  }
  
  # Initialize centromere dataframe
  WkCentroDf <- DfIndex[,c(2,16)]
  
  # Find copies intersection
  WkCentroDf$Copies <- sapply(1:nrow(WkCentroDf), function(i){
    CEN_Idx <- NULL
    Chr <- WkCentroDf$Chr_NC[i]
    CEN <- WkCentroDf$CEN_30000[i]
    DfChr <- WkDepthDf[WkDepthDf$Chr == Chr,]
    
    # If centromere does not have and index, search down
    while(length(CEN_Idx) == 0){
      CEN <- CEN - 1
      CEN_Idx <- DfChr[DfChr$Index == CEN,]$CopiesRaw
    }
    CEN_Idx
  })
  
  PlB <- PlB + geom_point(data = WkCentroDf, aes(x = CEN_30000, y = Copies), inherit.aes = FALSE, color = "black")
  
  PlC <- ggplot(WkAlleleDf, aes(Strain, Ratio)) +
    geom_violin(fill = "lightblue") +
    ylim(c(0, 1)) +
    theme_bw() +
    scale_x_discrete(labels = Strain)
  
  PlD <- ggplot(WkAlleleDf, aes(x = Index, y = Ratio)) +
    geom_point(shape = 1, size = 0.1, alpha = 0.5, color = "darkgray") +
    ylim(c(0, 1)) +
    geom_vline(xintercept = DfIndex$Lines,
               linetype = "dashed", color = "gray30") +
    scale_x_continuous(
      breaks = DfIndex$Mids,
      labels = as.roman(1:16),
      expand = c(0, 0)) +
    theme_bw() +
    xlab("Chromosome") +
    theme(plot.title = element_text(hjust = 0.5))
  # Make plot
  PlAll <- ggarrange(PlA, PlB, PlC, PlD, ncol = 2, nrow = 2, widths = c(2.5, 17.5)) + theme(plot.margin = margin(1, 0, 0, 0, unit = "cm"))
  
  TitleText <- paste0(i, " ", LabelID, " ",
                      if (!is.na(MeanP) && !is.na(RoundP)) {
                        paste0(MeanP, "n (", RoundP, "n)")
                      } else {
                        paste0(BaseP, "n")
                      },
                      if (!is.na(RoundP) && RoundP != BaseP) paste0(" [", BaseP, "n]") else "",
                      if (!is.na(BioP)) paste0(" ", BioP, "n") else "")
  
  PlAll <- annotate_figure(PlAll, top = text_grob(TitleText))
  
  FigureName <- paste0(OutDir, i, "_", LabelID, "_ABR.tiff")
  
  tiff(file = FigureName, compression = "lzw", width = 300, height = 150, units = "mm", res = 300)
  print(PlAll)
  dev.off()
}

# Calculate number of heterozygous sites
AlleleDf <- AlleleDf[,c(1,4)]
AlleleDf <- unique(AlleleDf)

# Add number of heterozygous sites after filtering
StrainsHet <- data.frame("Sample" = Samples,
                         "Heterozygous_sites_nb" = unlist(lapply(Samples, function(Strain){
                           nrow(AlleleDf[AlleleDf$Strain == Strain,])
                         })))

# Save
write.csv(StrainsHet, paste0("out_data/", Prefix, "_het_sites.csv"), row.names = FALSE)