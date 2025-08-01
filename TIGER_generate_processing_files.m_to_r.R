# %# % Generate fastq with all possible sequence reads of a given length per a reference genome

# % need a reference genome in Matlab (e.g. use fastaread to import data). Must be the same genome as used to align sequence data to. 
# % save as "<reference_name>_sequence" in folder "Alignability_and_GC_filters"
# % also save in the same file a variable "chr_length" that contains only the length in bps of each chromosome by order

# % need a gap file in Matlab (e.g. convert a UCSC gap file for the same genome build into Matlab, save chromosome, start, end

# % code generates a fastq file with all possible sequences ("reads") of a chosen length, not including gaps. Generates arbitrary read metadata


### clear;clc

# Clear all variables from the environment
rm(list = ls())

# Clear the console
cat("\014")

######################################################


### genome_build = 'hg19'; # % hg19, hg38, mm10, dm6   # % other genomes can also be used- need genome sequence and gap file, and possibly code tweeks depending on chromosome naming convention
### read_length = 100; # % should be no longer than the lengths of the reads in the sequencing data, but no need to make it longer than 100
### TIGER_folder = '/TIGER';  # % ! change as appropriate

# Define genome and read parameters
genome_build <- 'hg19'  # Options: 'hg19', 'hg38', 'mm10', 'dm6'
### Other genomes can also be used — need genome sequence and gap file,
### and possibly code tweaks depending on chromosome naming convention

read_length <- 100  # Should be no longer than the lengths of the reads in the sequencing data,
### but no need to make it longer than 100

TIGER_folder <- '/TIGER'  # Change as appropriate


###  try
###     eval(['cd ' TIGER_folder '/Alignability_and_GC_filters/']) 
###  catch
###      eval(['mkdir ' TIGER_folder '/Alignability_and_GC_filters/']) 
###      eval(['cd ' TIGER_folder '/Alignability_and_GC_filters/']) 
###  end


# Define path
alignability_path <- file.path(TIGER_folder, 'Alignability_and_GC_filters')

# Try to set working directory; create the directory if it doesn't exist
if (!dir.exists(alignability_path)) {
  dir.create(alignability_path, recursive = TRUE)
}

setwd(alignability_path)


######################################################
    
###  eval(['mkdir ' genome_build '_' num2str(read_length) 'bp']) 
###  eval(['cd '  genome_build '_' num2str(read_length) 'bp']) 


# Create subdirectory name (e.g., "hg19_100bp")
subfolder_name <- paste0(genome_build, "_", read_length, "bp")

# Full path to the subfolder
subfolder_path <- file.path(alignability_path, subfolder_name)

# Create the subfolder if it doesn't exist
if (!dir.exists(subfolder_path)) {
  dir.create(subfolder_path)
}

######################################################


# % load Gap file; structure of file should be: chromosome, start, end
###  switch genome_build
###      case 'hg19'
###          load Gap_hg19 Gap
###      case 'hg38'
###          load Gap_hg38 Gap 
###      case 'mm10'
###          load Mouse_Gap Gap
###      case 'dm6'
###          load dm6_Gap Gap
###  end


# Decide file name based on genome build
gap_file <- switch(genome_build,
                   'hg19' = 'Gap_hg19.txt',
                   'hg38' = 'Gap_hg38.txt',
                   'mm10' = 'Mouse_Gap.txt',
                   'dm6'  = 'dm6_Gap.txt',
                   stop("Unsupported genome build: ", genome_build))

# Read the file (assumes tab-delimited format and no header)
Gap <- read.table(gap_file, header = FALSE, col.names = c("chromosome", "start", "end"), sep = "\t", stringsAsFactors = FALSE)


######################################################


###  line4{1} = repmat('A',1,read_length); # % arbitrary metadata
###  for Chr = 1:max(Gap(:,1))
###      disp(['Chromosome ' num2str(Chr)])


# Create line4 as a list and fill first element with 'A' repeated read_length times
line4 <- list()
line4[[1]] <- paste(rep('A', read_length), collapse = '')

# Assuming Gap$chromosome is numeric; if it's character like "chr1", you may need a different approach
# Convert chromosome column to numeric if possible
chromosomes_numeric <- as.numeric(as.character(Gap$chromosome))

# Loop from 1 to max chromosome number (ignoring NAs)
max_chr <- max(chromosomes_numeric, na.rm = TRUE)
for (Chr in 1:max_chr) {
  message(paste("Chromosome", Chr))
}

######################################################    


### % load genome sequence  (loading the entire genome every time and only keeping one chromosome takes more time but saves computer memory) (file needs to be in the path or same folder)
###   eval(['load ' genome_build '_sequence Sequence']) 
###   Use_sequence = Sequence{Chr};
###   clear Sequence
    
# Construct filename based on genome_build
sequence_file <- paste0(genome_build, "_sequence.rds")

# Load the entire sequence list/object
Sequence <- readRDS(sequence_file)

# Extract sequence for chromosome 'Chr' (assuming Chr is numeric index)
Use_sequence <- Sequence[[Chr]]

# Remove Sequence from memory
rm(Sequence)


######################################################


###    # % remove gaps
###    G = Gap(Gap(:,1)==Chr,2:3);
###    for i = size(G,1):-1:1
###        Use_sequence(G(i,1)+1:G(i,2)) = [];
###    end
###    clear G i
  
# Filter Gap rows where chromosome equals Chr
# Assuming Gap is a data.frame with columns: chromosome, start, end
G <- Gap[Gap$chromosome == Chr, c("start", "end")]

# Convert Use_sequence to a vector of single characters if it isn't already
# For example, if Use_sequence is a string, split it into chars:
Use_sequence_vec <- strsplit(Use_sequence, NULL)[[1]]


# Remove gaps in reverse order
for (i in nrow(G):1) {
  start_pos <- G$start[i] + 1  # +1 to match MATLAB indexing (if needed)
  end_pos <- G$end[i]
  
  # Remove the substring from start_pos to end_pos
  if (start_pos <= length(Use_sequence_vec)) {
    # Limit end_pos to length of sequence to avoid out-of-bounds
    end_pos <- min(end_pos, length(Use_sequence_vec))
    
    Use_sequence_vec <- Use_sequence_vec[-(start_pos:end_pos)]
  }
}



# Convert back to a single string if needed
Use_sequence <- paste0(Use_sequence_vec, collapse = "")

# Clean up
rm(G, i, start_pos, end_pos, Use_sequence_vec)


######################################################

###   # % write fastq file
###   filename = cell2mat(['Chr' chrnum(Chr,genome_build) '_' num2str(read_length) 'bp_R1.fastq']);
###   fid = fopen(filename, 'w');
###   line1txt = ['@NB551191:189:HTHCVBGX5:1:11101:' num2str(Chr) ':# %d']; # % arbitrary initial parameters, followed by chromosome and coordinate


# Construct the filename
filename <- paste0("Chr", chrnum(Chr, genome_build), "_", read_length, "bp_R1.fastq")

# Open the file for writing
fid <- file(filename, open = "w")

# Prepare the line1txt string (with a placeholder for coordinate)
# R doesn't have direct printf-style placeholders in strings; you can use sprintf when you want to insert numbers later.
line1txt <- sprintf("@NB551191:189:HTHCVBGX5:1:11101:%s:# %%d", Chr)
# Note: %%d escapes %d to be part of the string, as in MATLAB.

# Don't forget to close the connection when done writing:
close(fid)


######################################################

###      # % generate sequence reads [ordered by location]
###      lh = waitbar(0,['Chr' chrnum(Chr)]);  # % chrnum coverts chromosome index number to naming convention (e.g. 23 to X)
###      cnt = 0;
###      for i = 1:read_length
###          waitbar(i/read_length,lh)
###          clear Fast_I
###          seq_length_use = length(Use_sequence)-rem((length(Use_sequence)-i+1),read_length);
###          seqs_temp = ( reshape( Use_sequence(i:seq_length_use),read_length,floor(seq_length_use/read_length) ) )';
###          Seqs = cellstr(seqs_temp);
###  
###          seq_texts = compose(line1txt,(cnt+1:cnt+length(seqs_temp))');
###          Fast_I(1:4:length(Seqs)*4,1) = seq_texts;
###          Fast_I(2:4:length(Seqs)*4,1) = Seqs;
###          Fast_I(3:4:length(Seqs)*4,1) = {'+'}; # %line3; arbitrary
###          Fast_I(4:4:length(Seqs)*4,1) = line4;
###          
###          fprintf(fid,'# %s\n',Fast_I{:});
###          
###          cnt = cnt+length(seqs_temp);
###      end
###      close(lh)
###      eval(['gzip ' filename]) # % gzips the file - input for bwa is gzipped
###      eval(['delete ' filename]) # % delete the original
###      fclose('all'); 
###      clear fid Fast_I A line1 cnt seqs_temp seq_length_use seq_texts lh Seqs filename i
###      
###  end

###   Loop over offsets i from 1 to read_length.
###   For each offset, chop the sequence into reads of length read_length.
###   Compose FASTQ read headers (line1txt) with incrementing coordinates.
###   Assemble FASTQ entries (4 lines each).
###   Write to file.
###   Show a progress bar.
###   Finally gzip and delete the original file.


library(progress)  # For progress bar
library(R.utils)   # For gzip compression

# Assume:
# - Use_sequence is a single string representing the DNA sequence for Chr
# - line1txt is a format string with a placeholder for coordinate (e.g. "@...:%d")
# - line4 is a character vector (length = 1) representing the quality string, same length as read_length

# Convert Use_sequence string to vector of single characters
Use_sequence_vec <- strsplit(Use_sequence, NULL)[[1]]

seq_length <- length(Use_sequence_vec)

# Open connection to write FASTQ file
fid <- file(filename, open = "w")

# Progress bar
pb <- progress_bar$new(
  format = paste0("Chr", chrnum(Chr, genome_build), " [:bar] :percent eta: :eta"),
  total = read_length,
  clear = FALSE, width=60
)

cnt <- 0
for (i in 1:read_length) {
  pb$tick()

  # Calculate length to use, ensuring full reads
  seq_length_use <- seq_length - ((seq_length - i + 1) %% read_length)
  
  # Subset sequence vector for reads starting at offset i
  seq_indices <- i:seq_length_use
  n_reads <- length(seq_indices) / read_length
  
  # Reshape into matrix of reads (each row is one read)
  seqs_temp <- matrix(Use_sequence_vec[seq_indices], nrow = n_reads, ncol = read_length, byrow = TRUE)
  
  # Convert each row (read) to a single string
  Seqs <- apply(seqs_temp, 1, paste0, collapse = "")
  
  # Compose FASTQ read headers with coordinates
  coords <- cnt + 1:(length(Seqs))
  seq_texts <- sprintf(line1txt, coords)
  
  # Prepare FASTQ entries (4 lines per read)
  # line1: header
  # line2: sequence
  # line3: '+'
  # line4: quality (same for all reads)
  
  quality_string <- line4[[1]]  # assuming single string quality
  
  # Assemble all lines
  Fast_I <- character(length(Seqs) * 4)
  Fast_I[seq(1, length(Fast_I), by = 4)] <- seq_texts
  Fast_I[seq(2, length(Fast_I), by = 4)] <- Seqs
  Fast_I[seq(3, length(Fast_I), by = 4)] <- "+"
  Fast_I[seq(4, length(Fast_I), by = 4)] <- quality_string
  
  # Write to file
  writeLines(Fast_I, con = fid)
  
  cnt <- cnt + length(Seqs)
}

close(fid)

# gzip the file and remove the original
gzip(filename, remove = TRUE)


######################################################

###   eval(['cd ' TIGER_folder])

setwd(TIGER_folder)


# %# % Align fastq and extract alignment information
# % This part typically performed on a computer server


# % 1. copy fastq files to server

# % 2. use BWA-MEM to align to reference genome with duplicates marked. Should be the exact same reference used to generate the fastq files. 

# % 3. After alignment, extract the locations of reads that were uniquely aligned (change chromosomes and read length if needed):
# % for i in {1..22} X Y; do echo $i; samtools view -F 4 -F 16 -F 1024 -q 1  Chr${i}_100bp.bam ${i} | cut -f 4 > Chr${i}_100bp_filt.txt; done
# %  - or - (reads are labeled "chrXX", e.g. hg38, mm10)
# % for i in {1..22} X Y; do echo $i; samtools view -F 4 -F 16 -F 1024 -q 1  Chr${i}_100bp.bam chr${i} | cut -f 4 > Chr${i}_100bp_filt.txt; done
# %  - or - (drosophila)
# % for i in 2L 2R 3L 3R 4 X Y; do echo $i; samtools view -F 4 -F 16 -F 1024 -q 1  Chr${i}_100bp.bam chr${i} | cut -f 4 > Chr${i}_100bp_filt.txt; done

# % (-F 4: read unmapped; -F 16: read reverse strand (only the forward strand was written to fastq); -F 1024: PCR or optical duplicate; -q 1: only keep sequences with mapQ of 1 or more; -f 4: only extract the start position of the sequence)


# % * can then delete fastq files from computer and from server *



# %# % Import alignment samtools output to Matlab, save coordinates to remove 
# % this codes generates the alignability filter- list of coordinates that are not uniquely alignable


### clear;clc

# Clear all variables from the environment
rm(list = ls())

# Clear the console
cat("\014")


### genome_build = 'hg19'; # % hg19, hg38, mm10, dm6   # % other genomes can also be used- need genome sequence and gap file, and possibly code tweeks depending on chromosome naming convention
### read_length = 100; # % should be no longer than the lengths of the reads in the sequencing data, but no need to make it longer than 100
### TIGER_folder = '/TIGER';  # % ! change as appropriate

# Define genome and read parameters
genome_build <- 'hg19'  # Options: 'hg19', 'hg38', 'mm10', 'dm6'
### Other genomes can also be used — need genome sequence and gap file, and possibly code tweaks depending on chromosome naming convention

read_length <- 100  
# Should be no longer than the lengths of the reads in the sequencing data, but no need to make it longer than 100

TIGER_folder <- '/TIGER'  # Change as appropriate



###   eval(['cd ' TIGER_folder '/Alignability_and_GC_filters/' genome_build '_' num2str(read_length) 'bp']) 

# Construct the directory path
dir_path <- file.path(TIGER_folder, "Alignability_and_GC_filters", paste0(genome_build, "_", read_length, "bp"))

# Change working directory
setwd(dir_path)


######################################################

### # % load file with chromosome lengths
### eval(['load ' TIGER_folder '/Alignability_and_GC_filters/' genome_build '_sequence chr_length'])

# Construct file path
file_path <- file.path(TIGER_folder, "Alignability_and_GC_filters", paste0(genome_build, "_sequence.rds"))

# Load the object (assuming saved as .rds)
chr_length <- readRDS(file_path)

######################################################

###   # % remove "extra" chromosomes (e.g. unmapped contigs)
###   switch genome_build 
###       case 'hg19'
###           chr_length = chr_length(1:24);
###       case 'hg38'
###           chr_length = chr_length(1:24);
###       case 'mm10'
###           chr_length = chr_length(1:21);
###   end

chr_length <- switch(genome_build,
  'hg19' = chr_length[1:24],
  'hg38' = chr_length[1:24],
  'mm10' = chr_length[1:21],
  stop("Unsupported genome build: ", genome_build)
)

######################################################

# Define path
alignability_path <- file.path(TIGER_folder, 'Alignability_and_GC_filters')

# Try to set working directory; create the directory if it doesn't exist
if (!dir.exists(alignability_path)) {
  dir.create(alignability_path, recursive = TRUE)
}

setwd(alignability_path)


######################################################


###   for Chr = 1:length(chr_length)
###       disp(['Chromosome ' num2str(Chr)])
###       switch genome_build 
###           case 'mm10'
###               filename = cell2mat(['Chr' chrnum(Chr,'mm10') '_' num2str(read_length) 'bp_filt.txt']); % chrnum coverts chromosome index number to naming convention (e.g. 23 to X) 
###           case 'dm6'
###               filename = cell2mat(['Chr' chrnum(Chr,'dm6') '_' num2str(read_length) 'bp_filt.txt']);
###           otherwise
###               filename = cell2mat(['Chr' chrnum(Chr) '_' num2str(read_length) 'bp_filt.txt']);
###       end
###       fid = fopen(filename);
###       Coordinates = textscan(fid,'%f');  % read coordinate
###       fclose('all'); clear fid filename
###   
###       Coordinates_to_remove{Chr,1}(:,1) = setdiff(1:chr_length(Chr),Coordinates{1}); % all genomic coordinates that are not uniquely alignable
###       clear Coordinates
###   end


# Initialize list to store coordinates to remove
Coordinates_to_remove <- vector("list", length(chr_length))

for (Chr in seq_along(chr_length)) {
  message(paste("Chromosome", Chr))
  
  # Construct filename based on genome_build
  filename <- switch(genome_build,
    'mm10' = paste0("Chr", chrnum(Chr, "mm10"), "_", read_length, "bp_filt.txt"),
    'dm6'  = paste0("Chr", chrnum(Chr, "dm6"),  "_", read_length, "bp_filt.txt"),
    paste0("Chr", chrnum(Chr), "_", read_length, "bp_filt.txt")  # default case
  )
  
  # Read coordinates from the file
  # Assuming each line contains a coordinate number (no '#' prefix)
  Coordinates <- scan(filename, what = numeric(), quiet = TRUE)
  
  # Calculate coordinates to remove (those not uniquely alignable)
  Coordinates_to_remove[[Chr]] <- setdiff(seq_len(chr_length[Chr]), Coordinates)
  
  # Clean up
  rm(Coordinates)
}


######################################################

 
###   eval(['save -v7.3 Coordinates_to_remove_' genome_build '_' num2str(read_length) 'bp Coordinates_to_remove'])
###   eval(['cd ' TIGER_folder])

saveRDS(Coordinates_to_remove, file = paste0("Coordinates_to_remove_", genome_build, "_", read_length, "bp.rds"))
setwd(TIGER_folder)

######################################################

 
# %  * can then delete samtools files from server and from computer *

# %# % Define genomic DNA copy number windows from alignability filter coordinates to remove 
# % This code defines the windows in which read numbers will be counted ("read number windows").
# % smaller windows can be merged to larger ones (function "merge_windows") but not vice versa. Downside of using smaller windows is larger files, longer processing times and higher memory usage
 

###   clear;clc
###   
###   genome_build = 'hg19'; # % hg19, hg38, mm10, dm6
###   read_length = 100; 
###   win_size = 10000; # % this will be the size of the windows in uniquely-alignable bps
###   TIGER_folder = '/TIGER';  # % ! change as appropriate
###   
###   
###   eval(['cd ' TIGER_folder '/Alignability_and_GC_filters/' genome_build '_' num2str(read_length) 'bp']) 
###   eval(['load Coordinates_to_remove_' genome_build '_' num2str(read_length) 'bp']) # % load alignability filter

# Clear environment and console
rm(list = ls())
cat("\014")  # This sends Ctrl+L to clear console in RStudio

# Define variables
genome_build <- "hg19"  # hg19, hg38, mm10, dm6
read_length <- 100
win_size <- 10000       # size of windows in uniquely-alignable bps
TIGER_folder <- "/TIGER"  # change as appropriate

# Build the directory path and set working directory
dir_path <- file.path(TIGER_folder, "Alignability_and_GC_filters", paste0(genome_build, "_", read_length, "bp"))
setwd(dir_path)

# Load the RData file containing Coordinates_to_remove
load_file <- paste0("Coordinates_to_remove_", genome_build, "_", read_length, "bp.RData")
load(load_file)  # loads Coordinates_to_remove object into the environment


######################################################


###    # % load file with chromosome lengths
###    eval(['load ' TIGER_folder '/Alignability_and_GC_filters/' genome_build '_sequence chr_length'])
###    switch genome_build
###        case 'hg19'
###            chr_length = chr_length(1:24); 
###            load Gap_hg19 Gap
###        case 'hg38'
###            chr_length = chr_length(1:24); 
###            load Gap_hg38 Gap 
###        case 'mm10'
###            chr_length = chr_length(1:21);
###            load Mouse_Gap Gap
###        case 'dm6'
###            load dm6_Gap Gap        
###    end    


# Load chromosome length object
chr_file <- file.path(TIGER_folder, "Alignability_and_GC_filters", paste0(genome_build, "_sequence.RData"))
load(chr_file)  # loads `chr_length`

# Load gap file and subset chromosome lengths
gap_file <- switch(genome_build,
  "hg19" = {
    chr_length <- chr_length[1:24]
    "Gap_hg19.RData"
  },
  "hg38" = {
    chr_length <- chr_length[1:24]
    "Gap_hg38.RData"
  },
  "mm10" = {
    chr_length <- chr_length[1:21]
    "Mouse_Gap.RData"
  },
  "dm6" = "dm6_Gap.RData",
  stop("Unsupported genome build")
)

# Load the corresponding gap file (assumes it contains object named `Gap`)
load(file.path(TIGER_folder, "Alignability_and_GC_filters", gap_file))



######################################################


###   for Chr = 1:length(chr_length)
###       disp(['Chromosome ' num2str(Chr)])
###       kept_coordinates = setdiff(1:chr_length(Chr),Coordinates_to_remove{Chr});
###       Wins{Chr,1}(:,2) = kept_coordinates(1:win_size:end-win_size+1); # % start 
###       Wins{Chr,1}(:,3) = kept_coordinates(win_size+1:win_size:end)-1; # % end
###       Wins{Chr,1}(:,1) = Chr;
###       Wins{Chr,1}(:,4) = round(mean(Wins{Chr}(:,2:3)')'); # % window center
###       
###       # % remove windows that span gaps
###       G = Gap(Gap(:,1)==Chr,2:3);
###       in = [];
###       for i = 1:size(G,1)
###           in = [in;find(Wins{Chr,1}(:,2)<= G(i,1) & Wins{Chr,1}(:,3)>= G(i,2))];
###       end
###       Wins{Chr}(in,:) = [];
###       
###   end



Wins <- vector("list", length(chr_length))  # Initialize output

for (Chr in seq_along(chr_length)) {
  message(paste("Chromosome", Chr))
  
  # 1. Get alignable coordinates (remove unalignable)
  kept_coordinates <- setdiff(seq_len(chr_length[Chr]), Coordinates_to_remove[[Chr]])
  
  # 2. Ensure there’s enough data to make at least one window
  if (length(kept_coordinates) < win_size) next
  
  # 3. Define window start and end positions
  starts <- kept_coordinates[seq(1, length(kept_coordinates) - win_size + 1, by = win_size)]
  ends   <- kept_coordinates[seq(win_size + 1, length(kept_coordinates), by = win_size)] - 1
  
  # 4. Align start/end lengths
  n <- min(length(starts), length(ends))
  starts <- starts[1:n]
  ends   <- ends[1:n]
  
  # 5. Compute window centers
  centers <- round((starts + ends) / 2)
  
  # 6. Create window matrix: chr, start, end, center
  window_matrix <- cbind(chr = Chr, start = starts, end = ends, center = centers)
  
  # 7. Remove windows that span any gaps
  gaps_chr <- Gap[Gap[, 1] == Chr, 2:3, drop = FALSE]
  if (nrow(gaps_chr) > 0) {
    to_remove <- logical(nrow(window_matrix))
    for (i in seq_len(nrow(gaps_chr))) {
      gap_start <- gaps_chr[i, 1]
      gap_end   <- gaps_chr[i, 2]
      overlaps <- (window_matrix[, "start"] <= gap_start) & (window_matrix[, "end"] >= gap_end)
      to_remove <- to_remove | overlaps
    }
    window_matrix <- window_matrix[!to_remove, , drop = FALSE]
  }
  
  # 8. Save cleaned windows
  Wins[[Chr]] <- window_matrix
}



######################################################
 
###   eval(['save ' genome_build '_' num2str(read_length) 'bp_wins_'   num2str(win_size)  ' Wins'])
###   eval(['cd ' TIGER_folder])

# Define output file path
wins_filename <- file.path(
  TIGER_folder,
  paste0(genome_build, "_", read_length, "bp_wins_", win_size, ".RData") )

# Save the Wins object
save(Wins, file = wins_filename)

# (Optional) Set working directory to TIGER_folder
setwd(TIGER_folder)

######################################################

###   # %# % GC content normalization (files for calculating "expected" number of reads) 
###   # % Calculates GC content in 401 bp disregarding Ns
###   # % Saves files ("xbp_GC_cont_401_filtered_bins_chrx") with all the genomic coordinates (after alignability filter) that belong to each of the 401 GC# % bins
###   # % Saves "Reads_expected_nominal"- number of bps in each GC bin in each read number window
###   clear;clc
###   
###   genome_build = 'hg19'; # % hg19, hg38, mm10, dm6
###   read_length = 100;
###   win_size = 10000; # % only needed for the last part of the code
###   TIGER_folder = '/TIGER';  # % ! change as appropriate
###   
###   
###   eval(['cd ' TIGER_folder '/Alignability_and_GC_filters/' genome_build '_' num2str(read_length) 'bp']) 
###   eval(['load Coordinates_to_remove_' genome_build '_' num2str(read_length) 'bp']) # % load alignability filter


# Clear environment
rm(list = ls())

# Define parameters
genome_build <- "hg19"   # Options: hg19, hg38, mm10, dm6
read_length <- 100
win_size <- 10000        # Only needed later
TIGER_folder <- "/TIGER" # Change as appropriate

# Construct working directory path
work_dir <- file.path(TIGER_folder, "Alignability_and_GC_filters",
                      paste0(genome_build, "_", read_length, "bp"))

# Set working directory
setwd(work_dir)

# Load alignability filter file
coord_file <- paste0("Coordinates_to_remove_", genome_build, "_", read_length, "bp.RData")
load(coord_file)  # This loads Coordinates_to_remove into the environment



######################################################


##   # % load window coordinates
##   eval(['load ' genome_build '_' num2str(read_length) 'bp_wins_'   num2str(win_size)])
##   if strcmp(genome_build,'mm10')
##       Wins = Wins(1:20); # % no unique sequences on the Y chromosome
##   end

# Load saved Wins object
wins_file <- paste0(genome_build, "_", read_length, "bp_wins_", win_size, ".RData")
load(wins_file)  # Assumes this loads a list called 'Wins'

# If mm10, drop chromosome 21 (i.e., the Y chromosome)
if (genome_build == "mm10") {
  Wins <- Wins[1:20]  # Keep only Chr 1–20 (exclude Chr 21/Y)
}

######################################################

# % load genome sequence

###   eval(['load ' TIGER_folder '/Alignability_and_GC_filters/' genome_build '_sequence Sequence'])

seq_file <- file.path(TIGER_folder, "Alignability_and_GC_filters", paste0(genome_build, "_sequence.RData"))

# Load the Sequence object
load(seq_file)  # loads variable named 'Sequence'


######################################################

###   bins = 0:1/400:1;  
###   bins = round(bins*10000)/10000; # % prevents cases in which the value of bins isn't exact

bins <- seq(0, 1, by = 1/400)
bins <- round(bins * 10000) / 10000  # Prevent floating point precision issues


######################################################

###  for Chr = 1:size(Wins,1)
###      disp(['Chromosome ' num2str(Chr)])
###  

for (Chr in seq_along(Wins)) {
  message(paste("Chromosome", Chr))
  # loop body
}


######################################################


###    # %# %# %# %# %# %# %   calculate GC content for all bps in the genome [# % all alignable bp] # %# %# %# %# %
###    Use_sequence = Sequence{Chr};
###    GC_cont_win_chr(:,1) = 201:length(Use_sequence)-200;
### 
###    # % GC content for each bps
###    GC_sequence = zeros(length(Use_sequence),1);
###    GC_sequence(upper(Use_sequence)=='G' | upper(Use_sequence)=='C') = 1;
###    GC_sequence_sum = cumsum(GC_sequence);
###    GC_sequence_sum = [GC_sequence_sum(401); GC_sequence_sum(402:end)-GC_sequence_sum(1:end-401)];
###    GC_sequence_sum = GC_sequence_sum-GC_sequence(201:end-200); # % exclude the center position (the one actually tested)
###    
###    # % number of non-Ns in each GC window
###    Non_N = ones(length(Use_sequence),1);
###    Non_N(upper(Use_sequence)=='N') = 0;
###    Non_N_sum = cumsum(Non_N);
###    Non_N_sum = [Non_N_sum(401);Non_N_sum(402:end)-Non_N_sum(1:end-401)];
###    Non_N_sum = Non_N_sum-Non_N(201:end-200);
###    
###    GC_cont_win_chr(:,2) = NaN;
###    in = find(Non_N_sum==400); # % ignore GC windows that contain any Ns
###    GC_cont_win_chr(in,2) = GC_sequence_sum(in)./Non_N_sum(in);
###
###
###    # %# %# %# % filter for alignabilty  # %# %# %# %
###    [~, in] = intersect(GC_cont_win_chr(:,1),Coordinates_to_remove{Chr});
###    GC_cont_win_chr(in,:) = [];
###    in = find(isnan(GC_cont_win_chr(:,2)));
###    GC_cont_win_chr(in,:) = [];    
    

# Calculate GC content for all bps in the genome (all alignable bp)
Use_sequence <- Sequence[[Chr]]
GC_cont_win_chr <- matrix(NA, nrow = length(Use_sequence) - 400, ncol = 2)
GC_cont_win_chr[,1] <- 201:(length(Use_sequence) - 200)

# GC content for each bp
GC_sequence <- integer(length(Use_sequence))
Use_seq_upper <- toupper(Use_sequence)
GC_sequence[Use_seq_upper == 'G' | Use_seq_upper == 'C'] <- 1L
GC_sequence_sum <- cumsum(GC_sequence)

# Calculate sliding window sums of size 401 bp
GC_sequence_sum <- c(
  GC_sequence_sum[401],
  GC_sequence_sum[402:length(GC_sequence_sum)] - GC_sequence_sum[1:(length(GC_sequence_sum) - 401)]
)

# Exclude the center position in the window (position 201)
center_positions <- 201:(length(Use_sequence) - 200)
GC_sequence_sum <- GC_sequence_sum - GC_sequence[center_positions]

# Number of non-N bases in each GC window
Non_N <- integer(length(Use_sequence))
Non_N[Use_seq_upper != 'N'] <- 1L
Non_N_sum <- cumsum(Non_N)
Non_N_sum <- c(
  Non_N_sum[401],
  Non_N_sum[402:length(Non_N_sum)] - Non_N_sum[1:(length(Non_N_sum) - 401)]
)
Non_N_sum <- Non_N_sum - Non_N[center_positions]

# Initialize second column of GC_cont_win_chr to NA
GC_cont_win_chr[,2] <- NA_real_

# Indices where window has no Ns (exactly 400 non-N bases in window)
in <- which(Non_N_sum == 400)
GC_cont_win_chr[in, 2] <- GC_sequence_sum[in] / Non_N_sum[in]

# Filter for alignability: remove windows overlapping coordinates to remove
in_remove <- which(GC_cont_win_chr[,1] %in% Coordinates_to_remove[[Chr]])
if (length(in_remove) > 0) {
  GC_cont_win_chr <- GC_cont_win_chr[-in_remove, , drop = FALSE]
}

# Remove rows where GC content is NA
GC_cont_win_chr <- GC_cont_win_chr[!is.na(GC_cont_win_chr[,2]), , drop = FALSE]


######################################################

 
###    # % find the genomic coordinates that belong to each of the 401 GC# % bins
###    lh = waitbar(0,'wait');
###    for i = 1:length(bins)
###        waitbar(i/length(bins),lh)
###        in = find(GC_cont_win_chr(:,2)==bins(i));
###        GC_cont_win_chr_bins{i,1} = GC_cont_win_chr(in,1);
###    end
###    close(lh)
###   
###    eval(['save -v7.3 ' genome_build '_' num2str(read_length) 'bp_GC_cont_401_filtered_bins_chr' num2str(Chr) ' GC_cont_win_chr_bins'])
    


# Initialize progress bar
pb <- txtProgressBar(min = 0, max = length(bins), style = 3)

# Initialize list to store GC content bins
GC_cont_win_chr_bins <- vector("list", length(bins))

# Loop through each GC bin
for (i in 1:length(bins)) {
  # Update progress bar
  setTxtProgressBar(pb, i)
  
  # Find indices where GC content equals the current bin value
  in_bin <- which(GC_cont_win_chr[, 2] == bins[i])
  
  # Store the corresponding window centers in the list
  GC_cont_win_chr_bins[[i]] <- GC_cont_win_chr[in_bin, 1]
}

# Close progress bar
close(pb)

# Save the results
save(GC_cont_win_chr_bins, file = paste0(genome_build, "_", read_length, "bp_GC_cont_401_filtered_bins_chr", Chr, ".RData"))
  

######################################################


###       # %# %# %# %# %# % count number of bps in each GC bin in each read number window- creates variable "Reads_expected_nominal"  # %# %# %# %# %
###       win_coord = [Wins{Chr}(:,2);Wins{Chr}(end,3)];
###       for GC_bin = 1:401
###           Reads_expected_nominal{Chr,1}(:,GC_bin)  =  histcounts(GC_cont_win_chr_bins{GC_bin},win_coord); # % in each TIGER window, how many bp are found in each GC bin (in the genome)
###       end
###       Reads_expected_nominal{Chr}(:,402) = nansum(Reads_expected_nominal{Chr}'); # % this is for testing purposes only. Gets overwritten later. column 402 should be the same as the win_size for most or all windows
###       
###       clear GC_cont_win_chr GC_cont_win_chr_bins
###   end


# Define window coordinates
win_coord <- c(Wins[[Chr]][, 2], Wins[[Chr]][nrow(Wins[[Chr]]), 3])

# Initialize list to store expected reads
Reads_expected_nominal <- vector("list", length = length(chr_length))

# Loop through each GC bin
for (GC_bin in 1:401) {
  # Find indices of windows that fall into the current GC bin
  in_bin <- which(GC_cont_win_chr_bins[[GC_bin]] %in% win_coord)
  
  # Compute histogram counts for the current GC bin
  hist_counts <- hist(GC_cont_win_chr_bins[[GC_bin]][in_bin], breaks = win_coord, plot = FALSE)$counts
  
  # Store the counts in the list
  Reads_expected_nominal[[Chr]][, GC_bin] <- hist_counts
}

# Sum across rows to get total reads per window
Reads_expected_nominal[[Chr]][, 402] <- rowSums(Reads_expected_nominal[[Chr]], na.rm = TRUE)

# Clear unnecessary variables
rm(list = c("GC_cont_win_chr", "GC_cont_win_chr_bins"))


######################################################


###   eval(['save -v7.3 ' genome_build '_' num2str(read_length) 'bp_' num2str(win_size) 'bp_wins_Reads_expected_nominal_401 Reads_expected_nominal'])   
###   eval(['cd ' TIGER_folder])


# Define the file path
file_path <- file.path(TIGER_folder, paste0(genome_build, "_", read_length, "bp_", win_size, "bp_wins_Reads_expected_nominal_401.rds"))

# Save the object
saveRDS(Reads_expected_nominal, file = file_path)
 
# Change the working directory
setwd(TIGER_folder)


######################################################


# %# % Optional: calculate Reads_expected_nominal for additional windows sizes
# % For processing files with non-default read number windows, can run just this (GC# % coordinate files stay the same) 
###   clear;clc
###   
###   genome_build = 'hg19'; # % hg19, hg38, mm10, dm6
###   read_length = 100;
###   win_size = 10000; 
###   TIGER_folder = '/TIGER';  # % ! change as appropriate
###   
###   eval(['cd ' TIGER_folder '/Alignability_and_GC_filters/' genome_build '_' num2str(read_length) 'bp']) 


# Set parameters
genome_build <- "hg19"  # Options: hg19, hg38, mm10, dm6
read_length <- 100
win_size <- 10000
TIGER_folder <- "/TIGER"  # Adjust as necessary

# Define the directory path
dir_path <- file.path(TIGER_folder, "Alignability_and_GC_filters", paste0(genome_build, "_", read_length, "bp"))

# Change to the specified directory
setwd(dir_path)

# Load the required data
load(paste0("Coordinates_to_remove_", genome_build, "_", read_length, "bp.RData"))


######################################################


# % load window coordinates
##   eval(['load ' genome_build '_' num2str(read_length) 'bp_wins_'   num2str(win_size)])
##   if strcmp(genome_build,'mm10')
##       Wins = Wins(1:20); # % no unique sequences on the Y chromosome
##   end

# Load saved Wins object
wins_file <- paste0(genome_build, "_", read_length, "bp_wins_", win_size, ".RData")
load(wins_file)  # Assumes this loads a list called 'Wins'

# If mm10, drop chromosome 21 (i.e., the Y chromosome)
if (genome_build == "mm10") {
  Wins <- Wins[1:20]  # Keep only Chr 1–20 (exclude Chr 21/Y)
}
 
######################################################


###   for Chr = 1:size(Wins,1)
###       disp(['Chromosome ' num2str(Chr)])
###       eval(['load ' genome_build '_' num2str(read_length) 'bp_GC_cont_401_filtered_bins_chr' num2str(Chr)])
###       win_coord = [Wins{Chr}(:,2);Wins{Chr}(end,3)];
###       for GC_bin = 1:401
###           Reads_expected_nominal{Chr,1}(:,GC_bin)  =  histcounts(GC_cont_win_chr_bins{GC_bin},win_coord); # % in each TIGER window, how many bp are found in each GC bin (in the genome)
###       end
###       Reads_expected_nominal{Chr}(:,402) = nansum(Reads_expected_nominal{Chr}'); # % this is for testing purposes only. Gets overwritten later
###       
###       clear GC_cont_win_chr GC_cont_win_chr_bins
###   end


for (Chr in seq_len(nrow(Wins))) {
  cat("Chromosome", Chr, "\n")
  
  file_name <- paste0(genome_build, "_", read_length, "bp_GC_cont_401_filtered_bins_chr", Chr, ".RData")
  load(file_name)  # Assumes the data is saved in .RData format
  
  win_coord <- c(Wins[[Chr]][,2], Wins[[Chr]][nrow(Wins[[Chr]]),3])
  
  Reads_expected_nominal[[Chr]] <- matrix(0, nrow = length(win_coord) - 1, ncol = 402)
  
  for (GC_bin in 1:401) {
    Reads_expected_nominal[[Chr]][, GC_bin] <- hist(GC_cont_win_chr_bins[[GC_bin]], breaks = win_coord, plot = FALSE)$counts
  }
  
  Reads_expected_nominal[[Chr]][, 402] <- rowSums(Reads_expected_nominal[[Chr]], na.rm = TRUE)
  
  rm(GC_cont_win_chr, GC_cont_win_chr_bins)
}


######################################################

###   eval(['save -v7.3 ' genome_build '_' num2str(read_length) 'bp_' num2str(win_size) 'bp_wins_Reads_expected_nominal_401 Reads_expected_nominal'])   
###   eval(['cd ' TIGER_folder])

# Define the output file path
output_file <- file.path(TIGER_folder, paste0(genome_build, "_", read_length, "bp_", win_size, "bp_wins_Reads_expected_nominal_401.RData"))

# Save the Reads_expected_nominal object
save(Reads_expected_nominal, file = output_file)

# Set the working directory to the TIGER folder
setwd(TIGER_folder)
 
