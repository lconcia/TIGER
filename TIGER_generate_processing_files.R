# ---- Setup and Parameters ----
library(Biostrings)
library(data.table)
library(stringr)
library(R.utils)

genome_build <- "hg19"
read_length <- 100
TIGER_folder <- "/TIGER"

# ---- Helper: chromosome naming ----
chrnum <- function(chr, genome_build = "hg19") {
  # Map integer to chr string, e.g. 23 -> X
  if (chr == 23) return("X")
  if (chr == 24) return("Y")
  return(as.character(chr))
}

# ---- 1. Generate FASTQ with all possible reads ----
# Assumes you have: 
# - an RData file with a named list "Sequence" (chromosome strings) and "chr_length" (integer vector, by index)
# - a gap data.frame or matrix Gap: columns (chromosome, start, end)

setwd(file.path(TIGER_folder, "Alignability_and_GC_filters"))
outdir <- sprintf("%s_%dbp", genome_build, read_length)
if (!dir.exists(outdir)) dir.create(outdir)
setwd(outdir)

# Load gap file
gap_file <- switch(genome_build,
  "hg19" = "Gap_hg19.RData",
  "hg38" = "Gap_hg38.RData",
  "mm10" = "Mouse_Gap.RData",
  "dm6"  = "dm6_Gap.RData"
)
load(gap_file) # expects Gap variable: data.frame/matrix with columns: chromosome, start, end

# Load sequence/chr_length
load(sprintf("%s_sequence.RData", genome_build)) # expects Sequence (list of DNAStrings), chr_length (integer vector)

arbitrary_qual <- paste(rep("A", read_length), collapse = "")

for (Chr in unique(Gap[,1])) {
  cat("Chromosome", Chr, "\n")
  Use_sequence <- Sequence[[Chr]]
  G <- Gap[Gap[,1] == Chr, 2:3, drop=FALSE]
  
  # Remove gaps (set to N)
  for (i in seq(nrow(G),1)) {
    subseq(Use_sequence, G[i,1]+1, G[i,2]) <- DNAString(paste(rep("N", G[i,2]-G[i,1]+1), collapse=""))
  }
  seq_len <- length(Use_sequence)
  
  # Write FASTQ file
  fq_file <- sprintf("Chr%s_%dbp_R1.fastq", chrnum(Chr, genome_build), read_length)
  fq_conn <- file(fq_file, "w")
  nreads_total <- 0
  pb <- txtProgressBar(min=1, max=read_length, style=3)
  
  for (i in 1:read_length) {
    setTxtProgressBar(pb, i)
    # Make sure length is divisible by read_length
    seq_index <- i:(seq_len - ((seq_len-i+1) %% read_length))
    if (length(seq_index) < read_length) next
    matrix_reads <- matrix(as.character(Use_sequence[seq_index]), ncol=read_length, byrow=TRUE)
    seqs <- apply(matrix_reads, 1, paste, collapse="")
    for (j in seq_along(seqs)) {
      cat(sprintf("@NB551191:189:HTHCVBGX5:1:11101:%d:%d\n%s\n+\n%s\n",
                  Chr, nreads_total+j, seqs[j], arbitrary_qual),
          file=fq_conn)
    }
    nreads_total <- nreads_total + length(seqs)
  }
  close(pb)
  close(fq_conn)
  gzip(fq_file, destname=paste0(fq_file, ".gz"), overwrite=TRUE)
  file.remove(fq_file)
}

setwd(TIGER_folder)

# ---- 2. Import samtools output, generate alignability filter ----  
## line 119 in TIGER_generate_processing_files.m
%% Import alignment samtools output to Matlab, save coordinates to remove 
#  this codes generates the alignability filter- list of coordinates that are not uniquely alignable
# imports the output of generate_chromosome_mappability_mask.EDIT.sh

setwd(file.path(TIGER_folder, "Alignability_and_GC_filters", outdir))
load(sprintf("%s_sequence.RData", genome_build)) # loads chr_length
if (genome_build %in% c("hg19", "hg38")) chr_length <- chr_length[1:24]
if (genome_build == "mm10") chr_length <- chr_length[1:21]

Coordinates_to_remove <- vector("list", length(chr_length))
for (Chr in seq_along(chr_length)) {
  cat("Chromosome", Chr, "\n")
  filename <- sprintf("Chr%s_%dbp_filt.txt", chrnum(Chr, genome_build), read_length)
  if (!file.exists(filename)) stop(paste("Missing file:", filename))
  Coordinates <- scan(filename, what=integer(), quiet=TRUE)
  all_coords <- seq_len(chr_length[Chr])
  Coordinates_to_remove[[Chr]] <- setdiff(all_coords, Coordinates)
}
save(Coordinates_to_remove, file=sprintf("Coordinates_to_remove_%s_%dbp.RData", genome_build, read_length))
setwd(TIGER_folder)

# ---- 3. Define read number windows ----

setwd(file.path(TIGER_folder, "Alignability_and_GC_filters", outdir))
load(sprintf("Coordinates_to_remove_%s_%dbp.RData", genome_build, read_length))
load(sprintf("%s_sequence.RData", genome_build)) # chr_length
if (genome_build %in% c("hg19", "hg38")) chr_length <- chr_length[1:24]
if (genome_build == "mm10") chr_length <- chr_length[1:21]
load(gap_file) # Gap

Wins <- list()
for (Chr in seq_along(chr_length)) {
  cat("Chromosome", Chr, "\n")
  kept_coords <- setdiff(seq_len(chr_length[Chr]), Coordinates_to_remove[[Chr]])
  starts <- kept_coords[seq(1, length(kept_coords)-win_size+1, by=win_size)]
  ends   <- kept_coords[seq(win_size, length(kept_coords), by=win_size)]
  nwin <- min(length(starts), length(ends))
  wins_chr <- data.frame(chr=Chr, start=starts[1:nwin], end=ends[1:nwin],
                         center=round(rowMeans(cbind(starts[1:nwin],ends[1:nwin]))))
  # Remove windows that span gaps
  G <- Gap[Gap[,1]==Chr,2:3,drop=FALSE]
  if (nrow(G) > 0) {
    toremove <- integer(0)
    for (i in 1:nrow(G)) {
      in_gap <- which(wins_chr$start <= G[i,1] & wins_chr$end >= G[i,2])
      toremove <- c(toremove, in_gap)
    }
    if (length(toremove) > 0) wins_chr <- wins_chr[-unique(toremove),]
  }
  Wins[[Chr]] <- wins_chr
}
save(Wins, file=sprintf("%s_%dbp_wins_%d.RData", genome_build, read_length, win_size))
setwd(TIGER_folder)

# ---- 4. GC content normalization bins and Reads_expected_nominal ----

setwd(file.path(TIGER_folder, "Alignability_and_GC_filters", outdir))
load(sprintf("Coordinates_to_remove_%s_%dbp.RData", genome_build, read_length))
load(sprintf("%s_%dbp_wins_%d.RData", genome_build, read_length, win_size))
if (genome_build == "mm10") Wins <- Wins[1:20]
load(sprintf("%s_sequence.RData", genome_build)) # Sequence

bins <- round(seq(0, 1, by=1/400), 4)
Reads_expected_nominal <- list()
for (Chr in seq_along(Wins)) {
  cat("Chromosome", Chr, "\n")
  Use_sequence <- Sequence[[Chr]]
  nseq <- length(Use_sequence)
  GC_cont_win_chr <- data.frame(coord=201:(nseq-200), GC=NA)
  GC_vec <- as.numeric(Use_sequence == "G" | Use_sequence == "C" | Use_sequence == "g" | Use_sequence == "c")
  GC_cumsum <- cumsum(GC_vec)
  GC_window <- GC_cumsum[401:length(GC_cumsum)] - GC_cumsum[1:(length(GC_cumsum)-400)]
  NonN_vec <- as.numeric(Use_sequence != "N" & Use_sequence != "n")
  NonN_cumsum <- cumsum(NonN_vec)
  NonN_window <- NonN_cumsum[401:length(NonN_cumsum)] - NonN_cumsum[1:(length(NonN_cumsum)-400)]
  in_good <- which(NonN_window == 400)
  GC_cont_win_chr$GC[in_good] <- GC_window[in_good]/400
  
  # Filter for alignability
  coords_to_remove <- Coordinates_to_remove[[Chr]]
  GC_cont_win_chr <- GC_cont_win_chr[!(GC_cont_win_chr$coord %in% coords_to_remove),]
  GC_cont_win_chr <- GC_cont_win_chr[!is.na(GC_cont_win_chr$GC),]
  
  # Assign to GC bins
  GC_cont_win_chr_bins <- vector("list", length(bins))
  pb <- txtProgressBar(min=1, max=length(bins), style=3)
  for (i in seq_along(bins)) {
    setTxtProgressBar(pb, i)
    in_bin <- which(abs(GC_cont_win_chr$GC - bins[i]) < 1e-6)
    GC_cont_win_chr_bins[[i]] <- GC_cont_win_chr$coord[in_bin]
  }
  close(pb)
  save(GC_cont_win_chr_bins, file=sprintf("%s_%dbp_GC_cont_401_filtered_bins_chr%d.RData", genome_build, read_length, Chr))
  
  # Reads_expected_nominal
  win_coord <- c(Wins[[Chr]]$start, Wins[[Chr]]$end[nrow(Wins[[Chr]])])
  Reads_chr <- matrix(0, nrow=nrow(Wins[[Chr]]), ncol=401)
  for (GC_bin in 1:401) {
    Reads_chr[,GC_bin] <- hist(GC_cont_win_chr_bins[[GC_bin]], breaks=win_coord, plot=FALSE)$counts
  }
  Reads_chr <- cbind(Reads_chr, rowSums(Reads_chr))
  Reads_expected_nominal[[Chr]] <- Reads_chr
}
save(Reads_expected_nominal, file=sprintf("%s_%dbp_%dbp_wins_Reads_expected_nominal_401.RData",
                                          genome_build, read_length, win_size))
setwd(TIGER_folder)


# Explanation and notes:
# 
# - Each block translates the core MATLAB logic for that section.
# - Sequence should be a list of DNAString objects (from Biostrings); chr_length an integer vector.
# - Gap is expected to be a data.frame or matrix as in MATLAB.
# - This code does not include every possible error check (e.g., file existence, correct structure) but covers the main translation.
# - For the (optional) "additional window sizes" code, simply rerun the last block with different win_size and use existing GC bin files.
# - The alignment and filtering steps (samtools, BWA) remain as shell commands, not directly translatable to R.
# - If you need to read/write .mat files instead of .RData, use the R.matlab package.
