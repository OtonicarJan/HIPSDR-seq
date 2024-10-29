library(optparse)
library(HMMcopy)
library(gtools)
library(parallel)
# Define options for the script
option_list <- list(
  make_option(c("-f", "--file"), type = "character", default = NULL, 
              help = ".txt file with paths to cells", metavar = "character"),
  make_option(c("-o", "--out"), type = "character", default = NULL, 
              help = "output directory", metavar = "character"),
  make_option(c("-r", "--reference"), type = "character", default = NULL, 
              help = "path to the reference files", metavar = "character"),
  make_option(c("-w", "--bin_kb"), type = "integer", default = 1000, 
              help = "window size in kb", metavar = "integer"),
  make_option(c("-c", "--min_count"), type = "integer", default = 0, 
              help = "minimal read count per cell", metavar = "integer"),
  make_option(c("-l", "--long"), action = "store_true", default = FALSE,
              help = "compute long segments"),
  make_option(c("-b", "--bulk"), action = "store_true", default = FALSE,
              help = "create bulk segments"),
  make_option(c("-s", "--keep_segs"), action = "store_true", default = FALSE,
              help = "keep segment files"),
  make_option(c("-m", "--mu"), default = NULL,type = "double",
              help = "provide lower bound for mu to define deletions")
)

# Parse command line options
opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

print("--------------------------------------")
cat("Input Parameters:\n")
cat("File: ", opt$file, "\n")
cat("Output directory: ", opt$out, "\n")
cat("Reference: ", opt$reference, "\n")
cat("Bin size (kb): ", opt$bin_kb, "\n")
cat("Minimum counts per cell: ", opt$min_count, "\n")
cat("Long segments: ", opt$long, "\n")
cat("Bulk segments: ", opt$bulk, "\n")
cat("Keep segments: ", opt$keep_segs, "\n")
cat("Mu lower bound: ", opt$mu, "\n")

print("--------------------------------------")

# Retrieve options from command line
bin_size_kb <- opt$bin_kb
input_files <- opt$file
output_dir <- opt$out
path_ref <- opt$reference
create_bulk <- opt$bulk
min_count <- opt$min_count
mu <- opt$mu

# Define a function to get segments
get_segments <- function(corrected_copy, long, output_dir,bin_size_kb, prefix = "", keep_seg = TRUE, mu = NULL) {
  if ((long) & (is.null(mu)) & (bin_size_kb > 100)) {
    print("Generating short segments")
    default_param <- HMMsegment(corrected_copy, getparam = TRUE)
    longseg_param <- default_param
    longseg_param$e <- 0.995
    longseg_segments <- HMMsegment(corrected_copy, longseg_param, verbose = FALSE)
  } else if ((long) & (is.null(mu)) & (bin_size_kb <= 100)) {
    print("Genrating long segments")
    default_param <- HMMsegment(corrected_copy, getparam = TRUE)
    longseg_param <- default_param
    longseg_param$e <- 0.999999999999999
    longseg_param$strength <- 1e30
    longseg_segments <- HMMsegment(corrected_copy, longseg_param, verbose = FALSE)
  } else if ((long) & !(is.null(mu))) {
    default_param <- HMMsegment(corrected_copy, getparam = TRUE)
    longseg_param <- default_param
    longseg_param$e <- 0.999999999999999
    longseg_param$strength <- 1e30
    mus = longseg_param$mu
    longseg_param$mu = c(mu,head(mus, -1))
    
    longseg_segments <- HMMsegment(corrected_copy, longseg_param, verbose = FALSE)
  }
  else if (!(long) & (is.null(mu))) {
    default_param <- HMMsegment(corrected_copy, getparam = TRUE)
    mus = default_param$mu
    default_param$mu = c(mu,head(mus, -1))
    
    longseg_segments <- HMMsegment(corrected_copy, default_param, verbose = FALSE)
  }
  else {
    longseg_segments <- HMMsegment(corrected_copy, verbose = FALSE)
  }
  
  all.segments <- longseg_segments$segs
  
  path_segments <- file.path(output_dir, paste0(prefix, "segments_long.bed"))
  path_windows <- file.path(output_dir, paste0(prefix, "windows.bed"))
  path_bins <-  file.path(output_dir, paste0(prefix, "per_bin_calls_long.bed"))
  bin = bin_size_kb*1000
  write.table(all.segments, file = path_segments, quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
  system(paste0("bedtools makewindows -b ", path_segments, " -w ", format(bin, scientific = F), " > ", path_windows))
  system(paste0("bedtools intersect -a ", path_segments, " -b ", path_windows, " > ", path_bins))
  system(paste0("rm ", path_windows))
  if (!keep_seg){
    system(paste0("rm ", path_segments))
  }
}

# Define file paths for reference files
path_gc_seg <- file.path(path_ref, paste0("hg38.", bin_size_kb, "kb.map.seg"))
path_map_seg <- file.path(path_ref, paste0("hg38.", bin_size_kb, "kb.map.seg"))

# Make output directory
dir.create(output_dir, showWarnings = FALSE)

# Read input files
files <- read.table(input_files, header = FALSE)$V1
# Define a variable to store all cells' read counts
all_cells_count <- NULL

process_cells <- function(cell){
  split <- strsplit(cell, "/")[[1]]
  cell_new <- strsplit(split[length(split)], ".seg")[[1]][1]
  an.error.occured <- FALSE
  tryCatch({
    uncorrected_reads <- wigsToRangedData(cell, path_gc_seg, path_map_seg)
    total_reads = sum(uncorrected_reads$reads)
    stopifnot(total_reads>min_count)
    
  }, error = function(e) {
    an.error.occured <<- TRUE
    print(e)
  })
  if (!an.error.occured) {
    an.error.occured <- FALSE
    tryCatch({
      corrected_copy <- correctReadcount(uncorrected_reads)
      print("pass correction")
    }, error = function(e) {
      an.error.occured <<- TRUE
    })
    if (!an.error.occured) {
      corrected_copy$chr <- factor(corrected_copy$chr, mixedsort(unique(corrected_copy$chr)))
      
      an.error.occured <- FALSE
      tryCatch({
        longseg_segments <- get_segments(corrected_copy, opt$long, output_dir, bin_size_kb,
                                         paste0(cell_new, "_"), keep_seg = opt$keep_segs)
        print("pass segmentation")
      }, error = function(e) {
        an.error.occured <<- TRUE
      })
    }
  }
}

if (!create_bulk) {
  
  mclapply(files,process_cells, mc.cores = 8)
} else {
  i = 1
  for (cell in files) {
    split <- strsplit(cell, "/")[[1]]
    cell_new <- strsplit(split[length(split)], ".seg")[[1]][1]
    an.error.occured <- FALSE
    tryCatch({
      uncorrected_reads <- wigsToRangedData(cell, path_gc_seg, path_map_seg)
      total_reads = sum(uncorrected_reads$reads)
      stopifnot(total_reads>min_count)
      
    }, error = function(e) {
      an.error.occured <<- TRUE
      print(e)
    })
    if (!an.error.occured) {
      if (i == 1) {
        all_cells_count = uncorrected_reads
        print(i)
        i = i + 1
      } else {
        all_cells_count$reads = all_cells_count$reads + uncorrected_reads$reads
        print(i)
        i = i + 1
      }
    }
  }
  
  
  corrected_copy <- correctReadcount(all_cells_count)
  split = strsplit(input_files, "/")[[1]]
  fname = split[length(split)]
  samp = paste0(strsplit(fname,".t")[[1]][1],"_")
  get_segments(corrected_copy, opt$long, output_dir,bin_size_kb, prefix = samp, keep_seg = opt$keep_segs)
  
}
