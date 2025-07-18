
# Load required libraries
library(R.matlab)

TIGER_merge_windows <- function(...) {
  # merges read number windows to larger windows. 
  # requires at least the first three of following input variables:
  # 1. Data (normalized data, used for the locations (four columns))
  # 2. Count_Expected_Data- this is what will actually be merged
  # 3. Windows_to_combine- number of windows to merge together
  # 4. Optional: genome_build- used for removal of windows spanning gaps. hg19 by default
  # Returns window coordinates (Loc_combined) and Count_expected_Data after window merging (Dat_combined)
  # Code can be modified to provide sliding windows as well
  
  # Parse variable arguments
  args <- list(...)
  
  Data <- args[[1]]
  Count_Expected_Data <- args[[2]]
  Windows_to_combine <- args[[3]]
  
  if (length(args) < 4) {
    genome_build <- 'hg19'
  } else {
    genome_build <- args[[4]]
  }
  
  # Load gap data based on genome build
  switch(genome_build,
         'hg19' = {
           gap_data <- readMat("Gap_hg19.mat")
           Gap <- gap_data$Gap
         },
         'hg38' = {
           gap_data <- readMat("Gap_hg38.mat")
           Gap <- gap_data$Gap
         },
         'mm10' = {
           gap_data <- readMat("Mouse_Gap.mat")
           Gap <- gap_data$Gap
         },
         'dm6' = {
           gap_data <- readMat("dm6_Gap.mat")
           Gap <- gap_data$Gap
         })
  
  # Initialize output lists
  Loc_combined <- vector("list", length(Data))
  Dat_combined <- vector("list", length(Data))
  for (i in 1:length(Data)) {
    Dat_combined[[i]] <- vector("list", ncol(Count_Expected_Data))
  }
  
  # combine- no overlapping windows
  for (Chr in 1:length(Data)) {
    
    Win_locs <- Data[[Chr]][, 2:3]
    if ((nrow(Win_locs) %% Windows_to_combine) != 0) {
      rows_to_remove <- (nrow(Win_locs) - (nrow(Win_locs) %% Windows_to_combine) + 1):nrow(Win_locs)
      Win_locs <- Win_locs[-rows_to_remove, ]
    }
    
    # Initialize Loc_combined for this chromosome
    num_combined_windows <- floor(nrow(Win_locs) / Windows_to_combine)
    Loc_combined[[Chr]] <- matrix(0, nrow = num_combined_windows, ncol = 4)
    
    Loc_combined[[Chr]][, 2] <- Win_locs[seq(1, nrow(Win_locs), by = Windows_to_combine), 1]
    Loc_combined[[Chr]][, 3] <- Win_locs[seq(Windows_to_combine, nrow(Win_locs), by = Windows_to_combine), 2]
    Loc_combined[[Chr]][, 1] <- Chr
    Loc_combined[[Chr]][, 4] <- round(rowMeans(Loc_combined[[Chr]][, 2:3]))
    
    for (c in 1:ncol(Count_Expected_Data)) {
      CN_data <- Count_Expected_Data[[Chr, c]]
      if ((nrow(CN_data) %% Windows_to_combine) != 0) {
        rows_to_remove <- (nrow(CN_data) - (nrow(CN_data) %% Windows_to_combine) + 1):nrow(CN_data)
        CN_data <- CN_data[-rows_to_remove, ]
      }
      
      # Initialize Dat_combined for this chromosome and column
      Dat_combined[[Chr]][[c]] <- matrix(0, nrow = num_combined_windows, ncol = 2)
      
      Dat_combined[[Chr]][[c]][, 1] <- CN_data[seq(1, nrow(CN_data), by = Windows_to_combine), 1]
      for (k in 2:Windows_to_combine) {
        Dat_combined[[Chr]][[c]][, 1] <- Dat_combined[[Chr]][[c]][, 1] + CN_data[seq(k, nrow(CN_data), by = Windows_to_combine), 1]
      }
      Dat_combined[[Chr]][[c]][, 2] <- CN_data[seq(1, nrow(CN_data), by = Windows_to_combine), 2]
      for (k in 2:Windows_to_combine) {
        Dat_combined[[Chr]][[c]][, 2] <- Dat_combined[[Chr]][[c]][, 2] + CN_data[seq(k, nrow(CN_data), by = Windows_to_combine), 2]
      }
    }
    
    # remove windows that span gaps
    G <- Gap[Gap[, 1] == Chr, 2:3]
    in_indices <- c()
    for (i in 1:nrow(G)) {
      in_indices <- c(in_indices, which(Loc_combined[[Chr]][, 2] <= G[i, 1] & Loc_combined[[Chr]][, 3] >= G[i, 2]))
    }
    
    if (length(in_indices) > 0) {
      Loc_combined[[Chr]] <- Loc_combined[[Chr]][-in_indices, ]
      for (c in 1:ncol(Count_Expected_Data)) {
        Dat_combined[[Chr]][[c]] <- Dat_combined[[Chr]][[c]][-in_indices, ]
      }
    }
  }
  
  # plot original and merged windows
  if (FALSE) {
    # Load required plotting library
    library(graphics)
    
    # Create plot
    Chr <- 2
    samp <- 1
    
    plot(Data[[Chr]][, 4] / 1e6, Data[[Chr]][, samp + 4], 
         pch = 16, cex = 0.5, 
         xlab = "Position (Mb)", ylab = "Value")
    points(Loc_combined[[Chr]][, 4] / 1e6, 
           Dat_combined[[Chr]][[samp]][, 1] / Dat_combined[[Chr]][[samp]][, 2], 
           pch = 16, cex = 0.5, col = "red")
    
    legend("topright", legend = c("Original", "Combined"), 
           pch = 16, col = c("black", "red"))
  }
  
  # Return results
  return(list(Loc_combined = Loc_combined, Dat_combined = Dat_combined))
}





# Key translation decisions and adaptations:
# 
# Function Arguments: Matlab's varargin and nargin were translated to R's ... (ellipsis) and list(...) for variable arguments handling.
# 
# Cell Arrays: Matlab cell arrays {} were translated to R lists using vector("list", length) and [[]] indexing.
# 
# Matrix Operations:
# 
# Matlab's matrix slicing (:,2:3) became R's [, 2:3]
# Matlab's mod function became R's %% operator
# Matlab's round(mean(...)') became R's round(rowMeans(...))
# File Loading: Matlab's load function was replaced with R's readMat() from the R.matlab library to read .mat files.
# 
# Switch Statement: Matlab's switch-case syntax was preserved but adapted to R's switch function syntax.
# 
# Indexing: Both Matlab and R use 1-based indexing, so no changes were needed for array indices.
# 
# Sequence Generation: Matlab's 1:Windows_to_combine:end became R's seq(1, nrow(data), by = Windows_to_combine).
# 
# Find Function: Matlab's find() function was replaced with R's which() function.
# 
# Matrix Initialization: Matlab's automatic matrix sizing was replaced with explicit matrix initialization using matrix(0, nrow = ..., ncol = ...).
# 
# Plotting: The optional plotting code was translated to use R's base graphics functions with plot(), points(), and legend().
# 
# Return Values: R functions return values using return(list(...)) instead of Matlab's multiple output variables.
# 
# Dependencies Required:
# 
# R.matlab package for reading .mat files
# Base R graphics for plotting (if plotting section is enabled)
# The translation preserves all functionality including:
# 
# Variable argument handling
# Genome build switching
# Window merging logic
# Gap removal logic
# Optional plotting functionality
# All mathematical operations and matrix manipulations
