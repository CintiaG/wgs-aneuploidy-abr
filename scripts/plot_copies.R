# Detect segmental aneuploidies and plot copy number profiles (30 kb windows)
#
# This script integrates:
#   (1) per-sample depth histograms,
#   (2) per-sample sequencing depth computed in 30 kb windows,
#   (3) a sample metadata table containing baseline ploidy and depth mode,
#   (4) a chromosome index file with centromere positions and chromosome colors,
# to estimate chromosome and segmental copy number changes.
#
# Copy number is inferred by normalizing depth values by the genome-wide depth
# mode and baseline ploidy. Consecutive windows with deviating copy number are
# grouped into segments to identify chromosome-scale and segmental aneuploidies.
#
# Inputs:
#   1) Depth histogram file (no header):
#        Sample, Chr, Depth, Freq, Tot, Fract
#      e.g. out_data/sequencing_histograms.bed
#
#   2) Depth per 30 kb windows file (no header):
#        Sample, Chr, Start, End, Depth
#      e.g. out_data/sequencing_depth_per_30Kb.bed
#
#   3) Sample metadata CSV file:
#        Required columns:
#          - Sample
#          - Base_ploidy
#          - Depth_mode
#        Optional columns:
#          - Code (short sample identifier for plotting)
#          - Miss_Chr (manual copy number corrections)
#          - Mean_ploidy, Rounded_ploidy (used only for plot titles, if present)
#
#   4) Chromosome index file (tab-delimited, header required):
#        Must include: Chr_NC, Size, Colors, Lines_30000, Mids_30000, CEN_30000
#
# Outputs:
#   - figures/copies_<Prefix>/*_copies.tiff
#       Per-sample copy number and depth plots.
#
#   - out_data/<Prefix>_copy_number_aneu.csv
#       Copy number per chromosome and segmented region.
#
#   - out_data/<Prefix>_copy_number.csv
#       Genome-wide copy number and aneuploidy summary per sample.
#
#   - out_data/<Prefix>_aneuploid_regions.csv
#       Coordinates and properties of all detected copy number segments.
#
#   - out_data/<Prefix>_depth_per_aneu.csv
#       Depth and inferred copy number per 30 kb window.
#
#   - out_data/<Prefix>_histograms_aneu.csv
#       Histogram data restricted to regions used in aneuploidy detection.
#
# Usage:
#   Rscript used_scripts/plot_copies.R \
#     out_data/sequencing_histograms.bed \
#     out_data/sequencing_depth_per_30Kb.bed \
#     files/sequencing.csv \
#     files/chrom_index_cen.txt \
#     sequencing \
#     600 \
#     8
#
# Notes:
#   - The last argument sets the maximum copy number displayed in plots.
#   - If Code or ploidy summary columns are not provided, the script
#     automatically falls back to sample identifiers and baseline ploidy.


library(ggpubr)
library(dplyr)
library(ggrepel)

GetRows <- function(WkDf, Change, Thres){
  # Function to group windows by depth
  # Set empty variables
  Rows <- data.frame()
  AllRows <- data.frame()
  PrevChr <- ""
  PrevEnd <- 0
  Counter <- 1
  
  # Filter according to the hard threshold, that has to be estimated depending on ploidy
  if (Change == "increase"){
    WkDf <- WkDf[WkDf$Copies > Thres,]
  }
  if (Change == "decrease"){
    WkDf <- WkDf[WkDf$Copies < Thres,]
  }
  if (Change == "normal"){
    WkDf <- WkDf[WkDf$Copies < Thres[1],]
    WkDf <- WkDf[WkDf$Copies > Thres[2],]
  }
  
  # Make sure to execute this part only if the filtered data frame is not empty
  if (nrow(WkDf) > 0){
    for (j in 1:nrow(WkDf)){
      # Calculate the overlap margin with the previous end, which allows to account for gaps the size of the slide
      Margin <- PrevEnd + Slide
      Start <- WkDf[j,]$Start
      Index <- WkDf[j,]$Index
      Chr <- WkDf[j,]$Chr
      WkRow <- WkDf[j,]
      
      # This happens when there are no or small gaps
      if (Start <= Margin && Chr == PrevChr){
        # Group the row to the block
        Rows <- rbind(Rows, WkRow)
        # Re-start the block search if the start is after the slide threshold or it is in a different chromosome
      } else {
        # Filter extreme depth values
        Rows <- Rows[(Rows$Copies < 8 & Rows$Copies > 0),]
        # Keep only blocks larger than 5 windows
        if (nrow(Rows) > 5){
          # Get the first and last row information as the detected block
          if (Chr == PrevChr){
            Rows$Chr <- paste0(Rows$Chr, "block", Counter)
            # Keep track of the number of detected regions
            Counter <- Counter + 1
          } else {
            # Restart counter
            Counter <- 1
          }
          AllRows <- rbind(AllRows, Rows)
        }
        # Assign current value as the first value of the next detected block
        Rows <- WkRow
      }
      # Save previous values for next iteration
      PrevChr <- Chr
      PrevEnd <- WkDf[j,]$End
    }
    # Dump information of the last analyzed row
    Rows <- Rows[(Rows$Copies < 8 & Rows$Copies > 0),]
    if (nrow(Rows) > 5){
      AllRows <- rbind(AllRows, Rows)
    }
  }
   
  AllRows
}

GetRegions <- function(RowsDf, Change){
  # Function that returns regions data frame
  # Empty data frame
  RegionsDf <- data.frame()
  
  # Execute this code only if input data frame is not empty
  if (nrow(RowsDf) > 0){
    # Gather each chromosome information
    for (Chr in unique(RowsDf$Chr)){
      WkDf <- RowsDf[RowsDf$Chr == Chr,]
      Start <- head(WkDf$Start, 1)
      End <- tail(WkDf$End, 1)
      StartIdx <- head(WkDf$Index, 1)
      EndIdx <- tail(WkDf$Index, 1)
      Copies <- WkDf$Copies
      
      RegionsDf <- rbind(RegionsDf, data.frame("Chr" = Chr,
                                               "Start" = Start,
                                               "End" = End,
                                               "Size" = End - Start,
                                               "Start_index" = StartIdx,
                                               "End_index" = EndIdx,
                                               "Mid_index" = round(mean(c(StartIdx, EndIdx))),
                                               "Ymin" = round(min(Copies), 2),
                                               "Ymax" = round(max(Copies), 2),
                                               "Mid_y" = mean(Copies)))
    }
    # Add change
    RegionsDf$Change <- Change
    # Add aneuploidy size proportion according to chromosome size
    RegionsDf$Proportion <- unlist(lapply(1:length(RegionsDf$Chr), function(i){
      Chr <- sub("block[0-9]", "", RegionsDf$Chr[i])
      round(RegionsDf[i,]$Size / DfIndex[DfIndex$Chr_NC == Chr,]$Size, 2)
    }))
  }
  RegionsDf
}

RenameRows <- function(RowsDf, Change){
  # A function to rename the rows data frame according to the new established names
  unlist(lapply(RowsDf$Chr, function(Chr){
    RegionsDf[(RegionsDf$Block == Chr & RegionsDf$Change == Change),]$Chr_new
  }))
}

PostFilter <- function(Df){
  # A function to further filter incoming depth regions, which allows to avoid overlappings
  if (nrow(Df) > 0){
    Df <- Df %>%
      group_by(Chr) %>%
      slice(-c(1, n())) %>%
      ungroup()
  }
  Df
}

CorrectCopies <- function(Correction, Df, Column, Type){
  if (Correction != ""){
    CorrSep <- unlist(strsplit(Correction, ";"))
    for (Corr in CorrSep){
      CorrIdx <- as.numeric(sub(".*\\*", "", Corr))
      CorrNmb <- as.numeric(sub("\\*.*", "", Corr))
      WkChr <- DfIndex$Chr_NC[CorrIdx]
      # All rows equal to that chromosome
      Idx <- Df$Chr == WkChr
      if (Type == "First"){
        # If it is the first time we perform the operation
        Df[[Column]][Idx] <- Df[[Column]][Idx] + CorrNmb
      }
      if (Type == "Second"){
        # If it is the second time, return to original value
        Df[[Column]][Idx] <- Df[[Column]][Idx] - CorrNmb
      }
      if (Type == "Label"){
        Sign <- ifelse(CorrNmb > 0,
                        "\u2191", "\u2193")
        Df[[Column]][Idx] <- paste0(Sign, Df[[Column]][Idx])
      }
    }
  }
  Df
}

# Read input arguments
args <- commandArgs(trailingOnly = TRUE)

HistFile <- args[1]
DepthFile <- args[2]
DfFile <- args[3]
IndexFile <- args[4]
Prefix <- args[5]
OutDir <- paste0("figures/copies_", Prefix, "/")

dir.create(OutDir, recursive = TRUE, showWarnings = FALSE)
dir.create("out_data", recursive = TRUE, showWarnings = FALSE)

Max <- as.numeric(args[6])
Top <- as.numeric(args[7])
Slide <- 30000

# Read all reads histogram data frame
HistDf <- read.table(HistFile, col.names = c("Sample", "Chr", "Depth", "Freq", "Tot", "Fract"))
# Filter mitochondrial DNA
HistDf <- HistDf[HistDf$Chr != "ref|NC_001224|",]

# Read all reads depth data frame
DepthDf <- read.table(DepthFile, col.names = c("Sample", "Chr", "Start", "End", "Depth"))
# Remove mitochondrial DNA
DepthDf <- DepthDf[DepthDf$Chr != "ref|NC_001224|",]

# Read sample data frame to plot by order of sample
Df <- read.csv(DfFile)
if (!("Miss_Chr" %in% colnames(Df))) Df$Miss_Chr <- ""
Df$Miss_Chr[is.na(Df$Miss_Chr)] <- ""

Samples <- Df$Sample
Codes <- if ("Code" %in% colnames(Df)) Df$Code else rep(NA_character_, nrow(Df))
Ploidies <- Df$Base_ploidy
Corrections <- Df$Miss_Chr

# Read index and lines file
DfIndex <- read.table(IndexFile, header = TRUE, sep = "\t", quote = "\"",
                      stringsAsFactors = FALSE)
DfIndex <- DfIndex[DfIndex$Chr != "mitochondrial",]

# Normalize depth according to ploidy
Df$Norm <- Df$Depth_mode / Df$Base_ploidy

# Get genome size
GenomeSize <- sum(DfIndex$Size)

# Set empty variables
CopiesDf <- NULL
CopiesDfGenome <- NULL
AllRegionsDf <- NULL
AllDepths <- NULL
AllHists <- NULL

for (i in 1:length(Samples)){
  Sample <- Samples[i]
  Code <- Codes[i]
  Ploidy <- Ploidies[i]
  Correction <- Corrections[i]
  
  # Use Code if available, otherwise fallback to Sample
  LabelID <- if (!is.null(Code) && !is.na(Code) && Code != "") Code else Sample
  
  # Safe access to optional ploidy columns
  HasMeanPloidy <- "Mean_ploidy" %in% colnames(Df)
  HasRoundedPloidy <- "Rounded_ploidy" %in% colnames(Df)
  
  MeanP <- if (HasMeanPloidy) Df$Mean_ploidy[i] else NA
  RoundP <- if (HasRoundedPloidy) Df$Rounded_ploidy[i] else NA
  BaseP <- Df$Base_ploidy[i]
  
  # Subset one sample histogram
  WkHistDf <- HistDf[HistDf$Sample == Sample,]
  
  # Filter depth below 10 and frequency below user's max depth
  WkHistDf <- WkHistDf[WkHistDf$Depth > 10,]
  WkHistDf <- WkHistDf[WkHistDf$Depth < Max,]
  WkHistDf <- WkHistDf[! WkHistDf$Chr == "genome",]
  
  # Subset one sample depth
  WkDepthDf <- DepthDf[DepthDf$Sample == Sample,]
  
  # Add an index to the working depth data frame
  WkDepthDf$Index <- 1:nrow(WkDepthDf)
  
  # Remove the first and last three windows of each chromosome, which correspond to telomeres
  WkDepthDf <- WkDepthDf %>%
    group_by(Chr) %>%
    slice(-c(1, (n()-2):n())) %>%
    ungroup()
  
  # Get number of copies
  WkHistDf$Copies <- WkHistDf$Depth/Df$Norm[i]
  WkDepthDf$Copies <- WkDepthDf$Depth/Df$Norm[i]
  # Keep original number of copies
  WkDepthDf$CopiesRaw <- WkDepthDf$Depth/Df$Norm[i]
  
  # Fix the mis-detected chromosomes
  WkDepthDf <- CorrectCopies(Correction, WkDepthDf, "Copies", "First")
  
  # Calculate the threshold according to the ploidy
  UpThres <- (1 + round(0.7 * (1 / Ploidy), 2)) * Ploidy
  DownThres <- (1 - round(0.7 * (1 / Ploidy), 2)) * Ploidy
  
  # Group overlapping windows with changed depth
  UpRows <- GetRows(WkDepthDf, "increase", UpThres)
  DownRows <- GetRows(WkDepthDf, "decrease", DownThres)
  
  # Post filter regions to avoid overlapping in different change depth regions
  UpRows <- PostFilter(UpRows)
  DownRows <- PostFilter(DownRows)
  
  # Detect increased and decreased deviating depth regions
  UpRegions <- GetRegions(UpRows, "increase")
  DownRegions <- GetRegions(DownRows, "decrease")
  
  # Combine increase and decrease regions
  RegionsDf <- rbind(UpRegions, DownRegions)
  
  # Detect the chromosome indexes that have aneuploidies
  Idx <- WkDepthDf$Chr %in% unique(sub("block[0-9]", "", RegionsDf$Chr))
  
  # Filter depth data frame by only the chromosomes that have aneuploidies and group depth of their normal regions
  NormalRows <- GetRows(WkDepthDf[Idx,], "normal", c(UpThres, DownThres))
  NormalRows <- PostFilter(NormalRows)
  
  # Get normal regions only for the chromosomes that present aneuploidies
  NormalRegions <- GetRegions(NormalRows, "normal")
  
  # Combine also with the regions data frame
  RegionsDf <- rbind(RegionsDf, NormalRegions)
  
  # Copy detected block information to a new column to later get the general chromosome name
  RegionsDf$Block <- RegionsDf$Chr
  RegionsDf$Chr <- sub("block[0-9]", "", RegionsDf$Chr)
  
  # Sort data frame by chromosome and start position
  if (nrow(RegionsDf) > 0){
    RegionsDf$Sample <- Sample
    RegionsDf$Code <- Code
    RegionsDf <- RegionsDf[with(RegionsDf, order(Chr, Start)),]
  }
  
  # Modify name of chromosome according to the number of regions
  RegionsDf$Chr_new <- unlist(lapply(unique(RegionsDf$Chr), function(Chr){
    WkDf <- RegionsDf[RegionsDf$Chr == Chr,]
    Len <- nrow(WkDf)
    Len <- rep(Len, Len)
    # Left, center, right
    AllSuffixes <- c("L", "C", "R")
    Suffixes <- ifelse(Len == 3, AllSuffixes, ifelse(Len == 2, AllSuffixes[c(1,3)], ""))
    paste0(WkDf$Chr, Suffixes)
  }))
  
  # Add type of aneuploidy
  RegionsDf$Type <- unlist(lapply(RegionsDf$Chr, function(Chr){
    WkDf <- RegionsDf[RegionsDf$Chr == Chr,]
    Len <- nrow(WkDf)
    Type <- ifelse(Len > 1, "Complex", "Simple")
    Type
  }))
  
  # Fix y coordinates of mis-detected chromosomes
  RegionsDf <- CorrectCopies(Correction, RegionsDf, "Ymin", "Second")
  RegionsDf <- CorrectCopies(Correction, RegionsDf, "Ymax", "Second")
  RegionsDf <- CorrectCopies(Correction, RegionsDf, "Mid_y", "Second")
  
  # Rename the chromosomes. This allow to treat them as if they were separated chromosomes
  UpRows$Chr_new <- RenameRows(UpRows, "increase")
  DownRows$Chr_new <- RenameRows(DownRows, "decrease")
  NormalRows$Chr_new <- RenameRows(NormalRows, "normal")
  # Add new name column to chromosomes without aneuploidies, which is the same
  WkDepthDf$Chr_new <- WkDepthDf$Chr
  
  # Separate the data frame into chromosomes with and without aneuploidies with previous index
  # and combine with the increase, decrease and normal rows data frame
  WkDepthDf <- rbind(WkDepthDf[!Idx,], UpRows, DownRows, NormalRows)
  
  # Remove any remaining blocks name in data frame
  WkDepthDf$Chr <- sub("block[0-9]", "", WkDepthDf$Chr)
  
  # Calculate the copy number of each complete and segmented chromosome
  WkCopiesDf <- WkDepthDf %>%
    group_by(Chr, Chr_new) %>%
    summarize(Mean_Copies_exact = round(median(Copies, na.rm = TRUE), 2))
  
  # If mean copies exact is within the hard threshold use ploidy, if not, calculate
  WkCopiesDf$Mean_Copies <- ifelse((WkCopiesDf$Mean_Copies_exact < UpThres & WkCopiesDf$Mean_Copies_exact > DownThres),
                                   Ploidy,
                                   round(WkCopiesDf$Mean_Copies_exact))
  
  # Add the code
  WkCopiesDf$Code <- Code
  
  # Add chromosome size
  WkCopiesDf$Weight <- unlist(lapply(1:nrow(WkCopiesDf), function(i){
    Chr <- WkCopiesDf$Chr_new[i]
    Copies <- WkCopiesDf$Mean_Copies[i]
    # If chromosome has a unmodified name, it does not have aneuploidies
    ifelse(Chr %in% DfIndex$Chr_NC,
           # Use chromosome normal size
           Copies * DfIndex[DfIndex$Chr_NC == Chr,]$Size,
           # Use grouped window size
           Copies * RegionsDf[RegionsDf$Chr_new == Chr,]$Size)
  }))
  
  # Add mids for complex aneuploidies regions
  WkCopiesDf$Mids <- unlist(lapply(1:nrow(WkCopiesDf), function(i){
    Chr <- WkCopiesDf$Chr_new[i]
    ifelse(Chr %in% DfIndex$Chr_NC,
           DfIndex[DfIndex$Chr_NC == Chr,]$Mids_30000,
           RegionsDf[RegionsDf$Chr_new == Chr,]$Mid_index)
  }))
  
  # Add colors as well
  WkCopiesDf$Colors <- unlist(lapply(WkCopiesDf$Chr, function(Chr){
    DfIndex[DfIndex$Chr_NC == Chr,]$Colors
  }))
  
  # Labels for plotting
  WkCopiesDf$Label <- as.character(WkCopiesDf$Mean_Copies)
  WkCopiesDf <- CorrectCopies(Correction, WkCopiesDf, "Label", "Label")
  
  # Re-order data frame
  WkCopiesDf <- WkCopiesDf[,c(5,1:4,6:9)]
  
  # Estimate bioinformatic ploidy
  Ploidy_base <- round(Df$Base_ploidy[i])
  Ploidy_exact <- round(sum(WkCopiesDf$Weight) / GenomeSize, 2)
  
  # Determine aneuploidies
  Aneu <- unlist(lapply(1:nrow(WkCopiesDf), function(j){
    Diff <- WkCopiesDf$Mean_Copies[j] - Ploidy_base
    NameIdx <- match(WkCopiesDf$Chr[j], DfIndex$Chr_NC)
    Suffix <- sub(".*\\|", "", WkCopiesDf$Chr_new)[j]
    FullName <- paste0(NameIdx, Suffix)
    ifelse(Diff == 0, "", paste0(ifelse(Diff > 0, paste0("+", Diff), Diff), "*", FullName))
  }))
  
  # Add number of copies to copies data frame
  WkCopiesDf$Aneu <- Aneu
  RegionsDf$Aneu <- unlist(lapply(RegionsDf$Chr_new, function(Chr){
    WkCopiesDf[WkCopiesDf$Chr_new == Chr,]$Aneu
  }))
  
  # Filter aneuploidies
  Aneu <- Aneu[!Aneu == ""]
  
  # Add genome information 
  CopiesDfGenome <- rbind(CopiesDfGenome, data.frame("Sample" = Sample,
                                                     "Code" = Code,
                                                     "Sample_number" = i,
                                                     "Base_ploidy" = Ploidy_base,
                                                     "Bioinfo_ploidy" = Ploidy_exact,
                                                     "Bioinfo_ploidy_rounded" = round(sum(WkCopiesDf$Weight) / GenomeSize),
                                                     "Genome_Size" = round(sum(WkCopiesDf$Weight) / 1e6, 2),
                                                     "Aneuploidies" = Aneu))
  
  CopiesDf <- rbind(CopiesDf, WkCopiesDf)
  
  AllDepths <- rbind(AllDepths, WkDepthDf)
  
  AllHists <- rbind(AllHists, WkHistDf)
  
  AllRegionsDf <- rbind(AllRegionsDf, RegionsDf)
  
  PlA <- ggplot(WkHistDf, aes(Copies, Freq, color = Chr)) +
    geom_line(show.legend = FALSE, linewidth = 0.25) +
    geom_vline(xintercept = Df$Base_ploidy[i], linetype = "dashed", color = "red", linewidth = 0.25) +
    theme_bw() +
    scale_y_continuous(labels = function(x) x / 10000) +
    scale_x_continuous(limits = c(0, Top), labels = function(x){
      paste0(" ", x, ".0")
    }) +
    ggtitle("") +
    scale_color_manual(values = setNames(DfIndex$Colors, DfIndex$Chr_NC)) +
    coord_flip() +
    ylab("Freq (1e4)")
  
  PlB <- ggplot(WkDepthDf, aes(Index, CopiesRaw, color = Chr_new)) +
    geom_line(show.legend = FALSE, linewidth = 0.25) +
    scale_x_continuous(breaks = DfIndex$Mids_30000,
                       labels = paste0("Chr", 1:16),
                       expand = c(0, 0)) +
    geom_vline(xintercept = DfIndex$Lines_30000,
               linetype = "dashed", color = "darkgray") +
    theme_bw() +
    scale_color_manual(values = setNames(WkCopiesDf$Colors, WkCopiesDf$Chr_new)) +
    scale_y_continuous(limits = c(0, Top), labels = function(x){
      paste0(" ", x, ".0")
    }) +
    theme(plot.title = element_text(hjust = 0.5)) +
    ggtitle(paste0(
      i, " ", LabelID, " ",
      
      # Mean / rounded ploidy if available
      if (!is.na(MeanP) && !is.na(RoundP)) {
        paste0(MeanP, "n (", RoundP, "n)")
      } else {
        paste0(BaseP, "n")
      },
      
      # Base ploidy annotation if different
      if (!is.na(RoundP) && RoundP != BaseP) {
        paste0(" [", BaseP, "n]")
      } else {
        ""
      },
      
      # Bioinformatic ploidy
      " ", Ploidy_exact, "n"
    )) +
    xlab("Chromosome") +
    ylab("Copies") +
    # Plot the label of working copies
    geom_text(data = WkCopiesDf,
              aes(x = Mids,
                  y = 5,
                  label = Label),
              color = ifelse(WkCopiesDf$Mean_Copies == WkCopiesDf$Label,
                             "black",
                             "red"),
              show.legend = FALSE) +
    # Upper threshold
    geom_hline(yintercept = UpThres, linetype = "dashed", 
               color = "red") +
    # Lower threshold
    geom_hline(yintercept = DownThres, linetype = "dashed", 
               color = "red")
  
  # Make sure that regions were obtained
  if (nrow(RegionsDf) > 0){
    # Add aneuploidies
    PlB <- PlB + geom_rect(data = RegionsDf, aes(xmin = Start_index - 1, xmax = End_index + 1, ymin = Ymin - 0.01, ymax = Ymax + 0.01), alpha = 0.3, inherit.aes = FALSE, color = "red") +
      geom_text_repel(show.legend = FALSE, data = RegionsDf, aes(x = Mid_index, y = Mid_y, label = Aneu), color = "black", nudge_x = 10, nudge_y = -1)
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
  
  PlC <- ggarrange(PlA, PlB, ncol = 2, widths = c(2, 18))
  
  tiff(file = paste0(OutDir, i, "_", LabelID, "_copies.tiff"), compression = "lzw", width = 300, height = 75, units = "mm", res = 300)
  print(PlC)
  dev.off()
}

# Save
write.csv(CopiesDf, paste0("out_data/", Prefix, "_copy_number_aneu.csv"), row.names = FALSE)

write.csv(CopiesDfGenome, paste0("out_data/", Prefix, "_copy_number.csv"), row.names = FALSE)

write.csv(AllDepths, paste0("out_data/", Prefix, "_depth_per_aneu.csv"), row.names = FALSE)

write.csv(AllRegionsDf, paste0("out_data/", Prefix, "_aneuploid_regions.csv"), row.names = FALSE)

write.csv(AllHists, paste0("out_data/", Prefix, "_histograms_aneu.csv"), row.names = FALSE)