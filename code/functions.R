calculate_frip <- function(bam_files, peak_file){

  print("Reading Peak BED file")
  peaks <- rtracklayer::import(peak_file)
  
  print("Reading BAM files")
  reads <- as.list(bam_files) %>% 
    map(.f = readGAlignments)
  
  print("Extracting total reads")
  total_reads <- reads %>% 
    map_vec(.f = length)
  
  print("Calculating overlapping reads")
  overlap_reads <- reads %>% 
    map(.f = function(x) countOverlaps(subject = x, query = peaks)) %>% 
    map_vec(.f = sum, na.rm = TRUE)
  
  filenames <- basename(bam_files) 
  cell_id <- str_extract(string = filenames, pattern = "[0-9]{1,3}")
  
  print("Aggregating data")
  results <- tibble(cell = as.integer(cell_id), 
                    total_reads = total_reads, 
                    reads_in_peaks = overlap_reads) %>% 
    mutate(frip = reads_in_peaks/total_reads)
  
  return(results)
  
}


sample_genomic_intervals <- function(input_tibble, chrom_sizes){
  
  # chrom_sizes = tibble with 2 columns: chr and size
  # create tibble from chr, start, end
  
  regions = input_tibble %>% 
    dplyr::select(chr, start, end) %>% 
    mutate(width = abs(start - end)) %>% 
    left_join(chrom_sizes) %>% 
    mutate(sample_space = size - width)
  
  new_start <- apply(regions, 1, FUN = function(y) sample(x = 1:y[6], size = 1))
  new_end <- new_start + regions$width
  
  input_tibble <- input_tibble %>% 
    mutate(start = new_start, 
           end = new_end)
  
  return(input_tibble)
}



# File I/O ---------------------------------------------------------------------


## Single-cell MultiOmics output -----------------------------------------------

# primary SCMO output
parse_scmo_zmeth <- function(file, convert_meth = TRUE, includes_gene = FALSE, ...){
  
  if (includes_gene) {
    cols <- c("read_id", "chr", "start", "meth", "cutsite", "gene")
    classes <- "cciccc"
  } else {
    cols <- c("read_id", "chr", "start", "meth", "cutsite")
    classes <- "ccicc"
  }
  
  data <- data.table::fread(
    input = file,
    sep = "\t",
    col.names = cols, 
    ...)
  
  data <- as_tibble(data)
  
  if (convert_meth) {
    data <- data %>% 
      mutate(meth = ifelse(meth == "z", 0, 1))
  }
  
  return(data)
}


# read count stats from SCMO
parse_read_counts <- function(path){
  
  lib = str_extract(
    string = path,
    pattern = "(CG|JvB).*csv") %>% 
    str_replace(pattern = ".csv", replacement = "")
  
  data <- read_csv(path) %>% 
    janitor::clean_names() %>% 
    dplyr::rename("read" = x1) %>% 
    add_column(library = lib, .before = 1)
  
  return(data)
}

# convenience functions
extract_cell_id <- function(string){
  cell <- as.integer(str_replace_all(str_extract(string, pattern = "_[0-9]{1,3}"), "_", ""))
  return(cell)
}

extract_plate_id <- function(string){
  plate <- str_extract(string = string, pattern = "pl[0-9]{1,2}")
  return(plate)
}

extract_modification = function(string) {
  id = str_extract(string = string, pattern = "(k|K)[0-9]{1,2}(m|me)3")
  id = tolower(id)
  if(any(str_detect(id, pattern = "me"))) {
    id = id
  } else {
    id = str_replace_all(string = id, pattern = "m", replacement = "me") 
  }
  return(id)
}

extract_library = function(string) {
  lib = str_extract(string = string, pattern = "^[^_]*")
  return(lib)
}


## Miscellaneous ---------------------------------------------------------------


parse_mpileup <- function(path){
  
  mpileup_cols <- c("chrom", "pos", "ref_base", "coverage", "read_bases",
                    "base_quals", "read_coordinate", "cell_index", "umi")
  
  data <- data.table::fread(
    file = path, 
    sep = "\t", 
    header = FALSE, 
    col.names = mpileup_cols)
  
  data <- as_tibble(data)
  return(data)
}

  
  
parse_flagstat <- function(flagstat_file) {
  # Read the lines from the samtools flagstat output file
  lines <- readr::read_lines(flagstat_file)
  
  # Create a tibble, storing raw lines first
  data <- tibble(raw_line = lines) %>%
    # Extract the counts
    mutate(
      count1 = as.integer(str_extract(raw_line, "^\\d+")),
      count2 = as.integer(str_extract(raw_line, "(?<=\\+ )\\d+"))
    ) %>%
    # Remove the leading "<count1> + <count2> " portion
    mutate(
      remainder = str_remove(raw_line, "^\\d+ \\+ \\d+\\s+")
    ) %>%
    # Extract parenthetical text (if any)
    mutate(
      parentheses = str_extract(remainder, "\\(.*\\)$"),
      # Remove parentheses from the remainder
      parentheses = str_remove_all(parentheses, "^\\(|\\)$"),
      # Remove the parenthetical part from 'remainder'
      remainder = str_remove(remainder, "\\(.*\\)$")
    ) %>%
    # Rename 'remainder' to something more descriptive
    rename(remainder = "description") %>%
    # Attempt to parse a percentage if present inside the parentheses
    mutate(
      percentage = str_extract(parentheses, "[0-9.]+%")
    ) %>% 
    mutate(
      percentage = as.double(str_remove(percentage, "%"))
    ) %>% 
    mutate(
      # collapse description string
      description = str_replace_all(description, " ", "_")
    ) %>% 
    mutate(
      # remove trailing '_'
      description = str_replace(description, "_$", "")
    ) %>% 
    # Reorder or select columns as desired
    dplyr::select(count1, count2, percentage, description)
  
  return(data)
}

parse_rpe_frip <- function(path){
  filename <- basename(path)
  mod <- extract_modification(filename)
  plate <- extract_plate_id(filename)
  
  input <- read_tsv(file = path) %>% 
    mutate(
      cell = as.integer(str_extract(filename, pattern = "[0-9]{1,3}")), 
      mod = mod, 
      plate = plate
    ) %>% 
    dplyr::select(cell, mod, plate, total_reads, reads_in_peaks) %>% 
    mutate(
      frip = reads_in_peaks/total_reads
    )
  
  return(input)
  
}



## MACS output // narrow peaks file --------------------------------------------

read_narrow_peaks <- function(file, returnGranges = TRUE){
  data <- read.table(file)
  colnames(data) <- c('chr', 'start', 'end', 'name', 'score', 'strand', 
                      'signal', 'pval', 'qval', 'peak')
  if (returnGranges) {
    data = makeGRangesFromDataFrame(data, keep.extra.columns = TRUE)
    return(data)
  } else {
    return(data)
  }
}


## MACS output // broad peaks file ---------------------------------------------
parse_broad_peaks <- function(file, returnGranges = TRUE){
  data <- read.table(file)
  colnames(data) <- c('chr', 'start', 'end', 'name', 'score', 'strand', 
                      'signal', 'pval', 'qval')
  
  data <- data %>% 
    add_column(peak = NA)
  
  if (returnGranges) {
    data = makeGRangesFromDataFrame(data, keep.extra.columns = TRUE)
    return(data)
  } else {
    return(data)
  }
}



## ENCODE Methylbed files ------------------------------------------------------

parse_methyl_bed <- function(path){
  # Reads bed9 file used for methylation data from ENCODE
  # Silently drops columns used for visualization
  
  colnames <- c("chr", "start", "end", "name", "score", "strand", "display_start", 
                "display_end", "color", "cov", "meth", "ref", "sample", "quality")
  cols_keep <- c(1, 2, 3, 6, 10, 11)
  cols_drop <- c(4, 5, 7, 8, 9, 12, 13, 14)
  
  cat("Reading file... \n")
  data <- data.table::fread(file = path, header = FALSE, 
                            col.names = colnames[cols_keep], drop = cols_drop)
  cat("File read successfully, converting to tibble \n")
  data <- as_tibble(data)
  return(data)
}

convert_methyl_bed_to_rds <- function(path){
  cell_line <- str_split_1(path, pattern = "/")[5]
  file_name <- str_split_1(path, pattern = "/")[6] %>% 
    str_replace(pattern = ".bed.gz", replacement = "")
  cat(sprintf("Processing cell line %s, file %s \n", cell_line, file_name))
  
  data <- read_methyl_bed(path)
  cat("Done ! \n")
  out_path <- sprintf("./temp/encode_%s_wgbs_%s.rds", cell_line, file_name)
  cat("Saving to RDS")
  saveRDS(object = data, file = out_path)
}


## ENCODE Chip Data ------------------------------------------------------------

read_chip_bed <- function(x) read_table(x, col_names = c("chr", "start", "end", "count"))



## Zeller et al. 2022 data -----------------------------------------------------


# helper function to read files
read_zeller_k562 <- function(path){
  data <- read_csv(file = path, col_types = paste0(strrep("c", 1), strrep("n", 1154)))
  colnames(data)[1:3] <- c("chr", "start", "end")
  data <- data[-1, ]
  return(data)
}


# Miscellaneous ----------------------------------------------------------------

# overlap in bp between two granges objects
overlap_granges <- function(gr1, gr2) {
  overlap_bp <- GenomicRanges::intersect(gr1, gr2, ignore.strand = TRUE) %>% 
    as.data.frame() %>% 
    pull(width) %>% 
    sum()
  return(overlap_bp)
}


# calculate concordance between CpG sites (concordance = 1 fo all Os or all 1s)
concordance_pairwise <- function(x) {
  # make sure that beta values are in the range 0 .. 1 and that there is only 1
  # measurement per cell
  p <- mean(x)
  return(p^2 + (1 - p)^2)
}

# this concordance measure estimates what % of cells shares the same value
concordance_majority <- function(x) {
  # x is a numeric or integer vector of 0s and 1s
  p <- mean(x)               # proportion of 1s
  return(max(p, 1 - p))      # proportion of the majority class
}

scale_beta <- function(tibble, scale_to_100 = TRUE){
  # expects tibble with column 'meth'
  # if methylation is expressed as 0 to 1, convert to 0 to 100
  
  # input check
  if(!"meth" %in% colnames(tibble)) {
    stop("Input data must contain column 'meth'")
  }
  
  # get max beta value
  max_meth <- max(tibble$meth)
  
  # logic
  if (max_meth == 100) {
    if (scale_to_100) {
      cat("Data in range 0..100, target = 100, returning data unscaled \n")
      return(tibble)
    } else {
      cat("Data in range 0..100, target = 1, scaling down by 100... \n")
      tibble <- tibble %>% 
        mutate(meth = meth/100)
      return(tibble)
    }
  } else if (max_meth == 1) {
    if (scale_to_100) {
      cat("Data in range 0..1, target = 100, scaling up by 100... \n")
      tibble <- tibble %>% 
        mutate(meth = meth*100)
      return(tibble)
    } else {
      cat("Data in range 0..1, target = 1, returning data unscaled \n")
    }
  }
}


simulate_taps_detection <- function(x, alpha, beta) {
  x <- sapply(x, function(y){
    if (y == 1) {
      # True = 1
      rbinom(n = 1, size = 1, prob = alpha)
    } else {
      # True = 0
      rbinom(n = 1, size = 1, prob = beta)
    }
  })
}

prefix_chromosome_names <- function(tibble){
  if(!("chr" %in% colnames(tibble))){
    stop("Tibble should contain the column 'chr'")
  }
  
  if(any(str_detect(string = tibble$chr, pattern = "chr"))){
    cat("Chromosomes already prefixed with 'chr', returning data unchanged")
    return(tibble)
  } else {
    tibble <- tibble %>% 
      mutate(chr = paste0("chr", chr))
    return(tibble)
  }
}

prettify_mod <- function(x) return(paste0("H3", str_to_title(x)))


normalize_zscore <- function(vector){
  sd <- sd(vector, na.rm = TRUE)
  mean <- mean(vector, na.rm = TRUE)
  norm <- (vector-mean)/sd
  return(norm)
}
