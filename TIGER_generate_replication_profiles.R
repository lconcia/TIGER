# TIGER_generate_replication_profiles.R
# R translation of TIGER_generate_replication_profiles.m
# Note: This script depends on R packages: data.table, dplyr, stringr, parallel, ggplot2, and R.matlab (for .mat files).
# Some custom MATLAB functions (e.g., TIGER_segment_filt, chrnum, TIGER_last_autosome) must be reimplemented or ported to R.

library(data.table)
library(dplyr)
library(stringr)
library(parallel)
library(ggplot2)
library(R.matlab)

# ---- User parameters ----

TIGER_folder <- "/TIGER" # change as appropriate
genome_build <- "hg19"   # hg19, hg38, mm10, dm6
read_length <- 100       # refers to alignability filter file
win_size <- 10000        # size of the windows in uniquely-alignable bps

# ---- Load alignability and window files ----

coords_file <- sprintf("%s/Alignability_and_GC_filters/%s_%dbp/Coordinates_to_remove_%s_%dbp.mat", 
                       TIGER_folder, genome_build, read_length, genome_build, read_length)
wins_file <- sprintf("%s/Alignability_and_GC_filters/%s_%dbp/%s_%dbp_wins_%d.mat",
                     TIGER_folder, genome_build, read_length, genome_build, read_length, win_size)

Coordinates_to_remove <- readMat(coords_file)$Coordinates_to_remove
Wins <- readMat(wins_file)$Wins

if (genome_build == "mm10") {
  Wins <- Wins[1:20] # no unique Y
}

# ---- Find all sam_* directories ----

dirs <- list.dirs(path = ".", full.names = FALSE, recursive = FALSE)
names <- dirs[grepl("^sam_", dirs)]

# ---- Helper functions (to be implemented) ----

chrnum <- function(chr_idx, genome_build) {
  # Should return chromosome name as string for the genome build
  # Implement as needed
  return(as.character(chr_idx))
}

TIGER_last_autosome <- function(genome_build) {
  # Return the number of autosomes for the genome build
  if (genome_build %in% c("hg19", "hg38")) return(22)
  if (genome_build == "mm10") return(19)
  if (genome_build == "dm6") return(3) # Example for Drosophila
  stop("Unknown genome build for TIGER_last_autosome")
}

# ---- Main Data Import Loop ----

Reads <- vector("list", length(Wins))
Windowed_data <- vector("list", length(Wins))
data_to_segment <- vector("list", length(Wins))

for (samp in seq_along(names)) {
  dir.create(names[samp], showWarnings = FALSE)
  setwd(names[samp])
  
  for (Chr in seq_along(Wins)) {
    message(sprintf("Importing, sample %d, Chr %d", samp, Chr))
    filename <- sprintf("readsChr%s.txt", chrnum(Chr, genome_build))
    if (!file.exists(filename)) next
    ReadsT <- scan(filename, what = numeric(), quiet = TRUE)
    if (is.list(Coordinates_to_remove)) {
      in_idx <- ReadsT %in% Coordinates_to_remove[[Chr]]
    } else {
      in_idx <- ReadsT %in% Coordinates_to_remove
    }
    ReadsT <- ReadsT[!in_idx]
    if (is.null(Reads[[Chr]])) Reads[[Chr]] <- list()
    Reads[[Chr]][[samp]] <- ReadsT
    
    # count reads in windows
    if (samp == 1) {
      Windowed_data[[Chr]] <- Wins[[Chr]][, 1:4, drop = FALSE]
    }
    win_coord <- c(Wins[[Chr]][,2], Wins[[Chr]][nrow(Wins[[Chr]]),3])
    read_count <- hist(Reads[[Chr]][[samp]], breaks = win_coord, plot = FALSE)$counts
    Windowed_data[[Chr]] <- cbind(Windowed_data[[Chr]], read_count)
  }
  
  # prepare data for segmentation
  last_auto <- TIGER_last_autosome(genome_build)
  mean_coverage <- sapply(Windowed_data[1:last_auto], function(x) mean(x[, ncol(x)], na.rm = TRUE))
  mean_coverage <- mean(mean_coverage, na.rm = TRUE)
  
  for (Chr in seq_along(Wins)) {
    data_to_segment[[Chr]] <- Windowed_data[[Chr]][, 1:4, drop = FALSE]
    seg_col <- Windowed_data[[Chr]][, ncol(Windowed_data[[Chr]])] / mean_coverage * 2
    data_to_segment[[Chr]] <- cbind(data_to_segment[[Chr]], seg_col)
  }
  
  names[samp] <- substr(names[samp], 5, nchar(names[samp]))
  setwd("..")
}

params <- list(
  genome_build = genome_build,
  read_length = read_length,
  win_size = win_size,
  TIGER_folder = TIGER_folder
)

# ---- Segmentation ----
# Placeholder: TIGER_segment_filt needs to be implemented in R!
# filt_in, data_in_segs_corrected, segment_params <- TIGER_segment_filt(data_to_segment, NULL, NULL, params$genome_build)
# params$segment_params <- segment_params

# ---- Save data ----
fold_name <- basename(getwd())
# save(list = c("names", "Reads", "Windowed_data", "params", "filt_in"), 
#      file = sprintf("%s_%s_%dbp_filtered_reads_%dKb.RData", fold_name, params$genome_build, params$read_length, params$win_size/1000))

# ---- Plotting and QC ----
# The plotting code can be written using ggplot2, using Data_filt_segments, Windowed_data, and names.
# The entire plotting part can be ported as needed.

# ---- Continue with other sections as needed ----

# Note: Many MATLAB functions (e.g., csaps, segment_filt_v2) don't have direct R equivalents.
# Smoothing can be done using smooth.spline or similar functions in R.
# Save and load MATLAB .mat files using R.matlab::writeMat() and readMat().
# If you need a full, ready-to-run pipeline in R, each MATLAB-specific function must be reimplemented in R.

# ---- Final notes ----
# This is a structural translation; you must fill in the details for segmentation, smoothing, and plotting as needed for your application.


## Translated with CoPilot
# 
# Important notes:
# 
# - The above code provides a structural translation of the MATLAB pipeline to R, with comments and placeholders where necessary.
# - Many custom MATLAB functions (e.g., TIGER_segment_filt, chrnum, TIGER_last_autosome, segment_filt_v2, etc.) must be reimplemented in R.
# - For reading .mat files, the R.matlab package is used. For saving intermediate results, use .RData or .rds as preferred.
# - Plotting is suggested with ggplot2, but the full code is not included due to complexity and space.
# - Further work is needed for a fully functional R pipeline; this translation aims to provide a clear starting point for that process.
