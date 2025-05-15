library(tidyverse)
library(GenomicRanges)
library(biomaRt)
library(GenomicFeatures)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(org.Hs.eg.db)
library(Signac)
library(Seurat)
library(EnsDb.Hsapiens.v86)  # For gene annotation (assuming human hg38)
library(ggplot2)
library(patchwork)  # For combining plots if desired
library(BSgenome.Hsapiens.UCSC.hg38)

source("./scripts/functions.R")

genome_size = 3000000000

color_scheme <- c("#D674BA","#86BF88","#7A416A")
color_scheme_2 <- c("#0D3B66", "#FAF0CA", "#F95738", "#6C969D", "#2E5339")

setwd(dir = "/mnt/ssd/private/cgeisenberger/projects/sc-epi-seq/")


# General notes
#
# 1: chromosome names should be prefixed with 'chr'
# 2: beta values are in the range (0, 100)



# AVO // Stats across libraries ------------------------------------------------


## Mpileup files // sequencing errors -----------------------------------------

pileup_files <- list.files(
  path = "./input/avo/mpileup/",
  pattern = "*.pileup",
  full.names = TRUE)

pileup_data <- pileup_files %>% 
  map(.f = function(x){
    sprintf("Reading %s \n", x)
    lib = str_extract(string = x, pattern = "CG.*pl[0-9]{1,2}")
    data = parse_mpileup(x) %>% 
      mutate(library = lib) %>% 
      mutate(
        matches = str_count(read_bases, pattern = "[,.]"), 
        mismatches = coverage - matches, 
        mismatch_rate = mismatches/coverage, 
        context = paste0(lag(ref_base, 1), ref_base, lead(ref_base, 1))
      ) %>% 
      filter(str_detect(context, "N", negate = TRUE)) %>% 
      group_by(library, context) %>% 
      summarise(taps = sum(mismatches)/sum(coverage))
  })

pileup_data <- Reduce(f = bind_rows, x = pileup_data)

## Lambda spike-in conversion efficiencies (all libraries) ---------------------

conversion_files <- list.files(path = "./input/avo/conversion_efficiencies/", pattern = ".txt", recursive = TRUE, full.names = TRUE)

conversion_efficiency <- as.list(conversion_files) %>% 
  map(function(x){
    name = str_extract(x, pattern = "[^/]*$") %>% 
      str_replace(".txt", "")
    data = read_table(x) %>% 
      mutate(library = name)
    return(data)
  })

conversion_efficiency <- Reduce(f = bind_rows, x = conversion_efficiency)

#saveRDS(object = conversion_efficiency, file = "./temp/avo_all_conversion_efficiency.rds")


## Mapping statistics (scMO librarystats) --------------------------------------

lib_stats <- list.files(path = "./input/avo/library_stats/", 
                        pattern = "csv", 
                        full.names = TRUE) %>% 
  map(.f = parse_read_counts) %>% 
  Reduce(f = bind_rows)

lib_stats <- lib_stats %>% 
  filter(read != "R?") %>% 
  mutate(demux_rate = demultiplexed_reads/raw_reads*100, 
         mapping_rate = mapped_reads/demultiplexed_reads*100, 
         dedup_rate = deduplicated_reads/assigned_site_reads*100) 
saveRDS(object = lib_stats, file = "./temp/avo_all_scmo_library_stats.rds")



# AVO // K562 cells ------------------------------------------------------------


## Prepare Signac object for single-gene plots ---------------------------------

# Parse metadata 

# generate GRanges object with 50 kb genomic bins
hg38_seqinfo <- seqinfo(BSgenome.Hsapiens.UCSC.hg38)[paste0("chr", 1:22)]
signac_bins <- tileGenome(hg38_seqinfo, tilewidth = 50000, cut.last.tile.in.chrom = TRUE)

# load adapter sequences (used as colnames in count matrices)
chic_adapters <- readxl::read_xlsx("./meta/sc-chic-adapters.xlsx") %>% 
  dplyr::select(id, well, barcode)

# detect Signac input files
signac_files <- list.files(path = "./input/avo/tagged_bam/", 
                           pattern = "fragments.sort.bed.gz$", 
                           recursive = TRUE, 
                           full.names = TRUE)
signac_data_sets <- c("H3K27me3", "H3K36me3", "H3K9me3", rep("H3K27me3", 4))
signac_libs <- str_extract(signac_files, "(CG|JvB)[^/]*")

# Create fragment objects and count matrices for each histone modification dataset
signac_frags <- signac_files %>% 
  map(.f = CreateFragmentObject, cells = chic_adapters$barcode)

signac_counts <- signac_frags %>% 
  map(.f = FeatureMatrix, features = signac_bins, cells = NULL)

chrom_assay <- map2(.x = signac_counts, 
                    .y = signac_frags, 
                    .f = function(x, y){ 
                      CreateChromatinAssay(counts = x, 
                                           min.cells = 10, 
                                           min.features = 1000,
                                           genome = hg38_seqinfo,
                                           assay = "peaks", 
                                           fragments = y,
                                           ranges = signac_bins)
                    })

seurat <- map(chrom_assay, CreateSeuratObject, assay = "peaks")

for (i in 1:length(seurat)){
  seurat[[i]]$dataset <- signac_data_sets[i]
}

seurat <- merge(
  x = seurat[[1]],
  y = seurat[2:7],
  add.cell.ids = signac_libs
)

seurat <- RunTFIDF(seurat)
seurat <- FindTopFeatures(seurat, min.cutoff = 20)
seurat <- RunSVD(seurat)
seurat <- RunUMAP(seurat, dims = 2:50, reduction = 'lsi')

DimPlot(seurat, group.by = 'dataset', pt.size = 0.5)


# add gene annotations 
gene_annot <- GetGRangesFromEnsDb(EnsDb.Hsapiens.v86) %>% 
  as_tibble()

gene_annot <- gene_annot %>% 
  mutate(seqnames = paste0("chr", seqnames)) %>% 
  dplyr::filter(seqnames %in% paste0("chr", 1:22)) %>% 
  makeGRangesFromDataFrame(keep.extra.columns = TRUE)

Annotation(seurat) <- gene_annot

# save object to disk
saveRDS(object = seurat, file = "./temp/avo_k562_seurat_all_datasets.rds")



## Epi-Seq data ----------------------------------------------------------------

# detect files, add library annotation
z_meth_anno <- tibble(
  path = list.files(path = "./input/avo/z_methylation/", full.names = TRUE, pattern = "tsv.gz"), 
  file = list.files(path = "./input/avo/z_methylation/", pattern = "tsv.gz"), 
  library = str_replace(file, ".tsv.gz", ""))

z_meth_anno <- left_join(
  z_meth_anno,
  library_anno)

# filter annotation
z_meth_anno_k562 <- z_meth_anno %>% 
  filter(cell == "k562")

# read data
z_meth_k562 <- z_meth_anno_k562$path %>% 
  map(.f = parse_scmo_meth) %>% 
  map(.f = function(x) {
    y = x %>% 
      mutate(
        library = extract_library(read_id), 
        cell = extract_cell_id(read_id))
    return(y)
  }) %>% 
  Reduce(f = bind_rows)

# add information about mod and plate no. 
z_meth_k562 <- left_join(z_meth_k562,
                         dplyr::select(.data = z_meth_anno_k562, c("library", "mod", "plate")))

z_meth_k562 <- z_meth_k562 %>% 
  dplyr::select(read_id, mod, plate, cell, chr, start, meth, cutsite)

# scale beta values to (0,100) and prefix chromosome names
z_meth_k562 <- z_meth_k562 %>% 
  scale_beta %>% 
  prefix_chromosome_names

# include end coordinate
z_meth_k562 <- z_meth_k562 %>% 
  add_column(end = (z_meth_k562$start + 1), .after = "start")

z_meth_k562 = readRDS("./temp/avo_k562_full_data.rds")



### QC stats per cell ----------------------------------------------------------

z_meth_k562_qc <- z_meth_k562 %>% 
  group_by(mod, plate, cell) %>% 
  dplyr::summarise(n_cpg = n(),
                   avg_beta = sum(meth)/n(),
                   unique_cuts = length(unique(read_id)),
                   ta_frac = sum(cutsite == "TA", na.rm = TRUE)/n(),
                   mod = sample(mod, size = 1)) %>%
  ungroup

# add metadata for QC cutoffs
qc_thresholds_k562 <- readxl::read_xlsx("./meta/qc_cutoffs_k562.xlsx")
z_meth_k562_qc <- left_join(z_meth_k562_qc, qc_thresholds_k562)

# determine which cells passed cutoffs
z_meth_k562_qc <- z_meth_k562_qc %>% 
  mutate(qc = ifelse(log10(unique_cuts) > n_min_log10 &
                       avg_beta > f_min &
                       avg_beta < f_max,
                     "passed",
                     "failed"))

# create QC filtered version of the dataset
z_meth_k562_filtered <- left_join(z_meth_k562,
                                  z_meth_k562_qc,
                                  by = c("mod", "plate", "cell"))

z_meth_k562_filtered <- z_meth_k562_filtered %>% 
  filter(qc == "passed") %>% 
  dplyr::select(read_id, mod, plate, cell, chr, start, meth, cutsite)

# save data to disk
#saveRDS(object = z_meth_k562_qc, file = "./temp/avo_k562_qc_stats.rds")
#saveRDS(object = z_meth_k562_filtered, file = "./temp/avo_k562_qc_filtered.rds")
#saveRDS(object = z_meth_k562, file = "./temp/avo_k562_full_data.rds")



## Fraction Reads in Peaks (FRiP) per cell -------------------------------------

# gather MACS3 peaks for each modification
frip_ref <- list(
  k27me3 = "./input/avo/macs3_peaks/k562/k27me3/macs3_peaks.broadPeak", 
  k36me3 = "./input/avo/macs3_peaks/k562/k36me3/macs3_peaks.broadPeak",
  k9me3 = "./input/avo/macs3_peaks/k562/k9me3/macs3_peaks.broadPeak", 
  k27me3 = "./input/avo/macs3_peaks/k562/k27me3/macs3_peaks.broadPeak",
  k27me3 = "./input/avo/macs3_peaks/k562/k27me3/macs3_peaks.broadPeak",
  k27me3 = "./input/avo/macs3_peaks/k562/k27me3/macs3_peaks.broadPeak",
  k27me3 = "./input/avo/macs3_peaks/k562/k27me3/macs3_peaks.broadPeak")

# list input directories - each contains a 'by_cell' subfolder
frip_bam_folders <- list.dirs(path = "input/avo/tagged_bam/", recursive = FALSE) %>% 
  paste0("/by_cell")

# extract library from folder name
frip_libraries <- str_extract(frip_bam_folders, pattern = "(CG|JvB)[^/]*")

# list cell-level BAMs for each folder
frip_bam_files <- frip_bam_folders %>% 
  map(function(x) list.files(path = x, pattern = "*.bam", full.names = TRUE))

# loop over list with BAM files and peaks, extract frip
frip_data <- map2(
  .x = frip_bam_files, 
  .y = frip_ref, 
  .f = calculate_frip)

# add library information to each element of the list
for (i in 1:length(frip_data)) {
  frip_data[[i]] <- frip_data[[i]] %>% 
    mutate(
      library = frip_libraries[i], 
      cell = as.integer(cell))
}

# aggregate into 1 tibble
frip_data <- frip_data %>% 
  Reduce(f = bind_rows)


frip_data <- left_join(frip_data, library_anno, by = "library")

frip_data <- frip_data %>% 
  dplyr::select(cell, total_reads, reads_in_peaks, frip, mod, plate)

saveRDS(object = frip_data, file = "./temp/avo_k562_frip_statistics.rds")

rm(frip_ref, frip_bam_folders, frip_libraries, frip_bam_files, frip_data)



## MACS output: CHiC peaks -----------------------------------------------------

avo_k562_peaks_k27 <- parse_broad_peaks(
  file = "./input/avo/macs3_peaks/k562/k27me3/macs3_peaks.broadPeak",
  returnGranges = FALSE)

avo_k562_peaks_k36 <- parse_broad_peaks(
  file = "./input/avo/macs3_peaks/k562/k36me3/macs3_peaks.broadPeak",
  returnGranges = FALSE)

avo_k562_peaks_k9 <- parse_broad_peaks(
  file = "./input/avo/macs3_peaks/k562/k9me3/macs3_peaks.broadPeak",
  returnGranges = FALSE)

avo_k562_peaks_k27 <- avo_k562_peaks_k27 %>% 
  mutate(mod = "k27me3")

avo_k562_peaks_k36 <- avo_k562_peaks_k36 %>% 
  mutate(mod = "k36me3")

avo_k562_peaks_k9 <- avo_k562_peaks_k9 %>% 
  mutate(mod = "k9me3")

avo_k562_peaks <- bind_rows(avo_k562_peaks_k27, avo_k562_peaks_k36, avo_k562_peaks_k9)

avo_k562_peaks <- avo_k562_peaks %>% 
  prefix_chromosome_names()

#saveRDS(object = avo_k562_peaks, file = "./temp/avo_k562_broad_peaks_combined.rds")
rm(avo_k562_peaks_k27, avo_k562_peaks_k36, avo_k562_peaks_k9)
rm(avo_k562_peaks)



# ENCODE // K562 cells ---------------------------------------------------------


## ChIP controls (for normalization) -------------------------------------------

# grep files
chip_controls <- list.files(path = "./input/encode/chip/control/", 
                            pattern = "counts_50k.txt", full.names = TRUE) 

# read data
chip_controls <- chip_controls %>% 
  map(.f = read_chip_bed) %>% 
  Reduce(f = function(x, y) full_join(x, y, b = c("chr", "start", "end")), x = .)

# aggregate data across controls
chip_controls <- chip_controls %>% 
  mutate(control = apply(chip_controls[, 4:6], 1, sum)) %>% 
  dplyr::select(chr, start, end, control)

# save to disk
#saveRDS(object = chip_controls, file = "./temp/encode_k562_chip_controls_aggregated.rds")



## Chip Data -------------------------------------------------------------------


# list all chip files
chip_files <- list.files(path = "./input/encode/chip/", 
                         pattern = "counts_50k.txt", full.names = TRUE, recursive = TRUE)

# remove controls
chip_files <- chip_files[str_detect(chip_files, pattern = "control", negate = TRUE)]

# read data
chip_data <- chip_files %>% 
  map(.f = read_chip_bed)

chip_data <- chip_data %>% 
  Reduce(f = function(x, y) full_join(x, y, by = c("chr", "start", "end")))

# replace colnames with histone mod
colnames(chip_data)[4:9] <- paste0(str_extract(chip_files, "h3k[0-9]{1,2}me3"), 
                                   paste0("_rep", rep(c(1, 2), 3)))

#saveRDS(object = chip_data, file = "./temp/encode_k562_chip_data_aggregated.rds")



## Whole Genome Bisulfite Sequencing (WGBS) ------------------------------------

# find all files
encode_wgbs_files <- list.files(path = "./input/encode/wgbs", pattern = "bed.gz", 
                                recursive = TRUE, full.names = TRUE)

# process
encode_wgbs_files %>% 
  as.list %>% 
  map(.f = convert_methyl_bed_to_rds)




# AvO // RPE1 cells ------------------------------------------------------------

## Prepare Signac object for single-gene plots ---------------------------------

# Parse metadata 

# generate GRanges object with 50 kb genomic bins
hg38_seqinfo <- seqinfo(BSgenome.Hsapiens.UCSC.hg38)[paste0("chr", 1:22)]
signac_bins <- tileGenome(hg38_seqinfo, tilewidth = 50000, cut.last.tile.in.chrom = TRUE)

# load adapter sequences (used as colnames in count matrices)
chic_adapters <- readxl::read_xlsx("./meta/sc-chic-adapters.xlsx") %>% 
  dplyr::select(id, well, barcode) %>% 
  rename(id = "cell")

# detect Signac input files
signac_files <- list.files(path = "./input/avo/rpe1_fragment_files/", 
                           pattern = "fragments.sorted.bed.gz$", 
                           recursive = TRUE, 
                           full.names = TRUE)
signac_data_sets <- rep(c("H3K27me3", "H3K36me3", "H3K9me3"), each = 4)
signac_libs <- str_extract(signac_files, "(CG|JvB)[^\\.]*")

# Create fragment objects and count matrices for each histone modification dataset
signac_frags <- signac_files %>% 
  map(.f = CreateFragmentObject, cells = chic_adapters$barcode)

signac_counts <- signac_frags %>% 
  map(.f = FeatureMatrix, features = signac_bins, cells = NULL)

chrom_assay <- map2(.x = signac_counts, 
                    .y = signac_frags, 
                    .f = function(x, y){ 
                      CreateChromatinAssay(counts = x, 
                                           min.cells = 10, 
                                           min.features = 1000,
                                           genome = hg38_seqinfo,
                                           assay = "peaks", 
                                           fragments = y,
                                           ranges = signac_bins)
                    })

seurat <- map(chrom_assay, CreateSeuratObject, assay = "peaks")

for (i in 1:length(seurat)){
  seurat[[i]]$dataset <- signac_data_sets[i]
}

seurat <- merge(
  x = seurat[[1]],
  y = seurat[2:length(seurat)],
  add.cell.ids = signac_libs
)

seurat <- RunTFIDF(seurat)
seurat <- FindTopFeatures(seurat, min.cutoff = 20)
seurat <- RunSVD(seurat)
seurat <- RunUMAP(seurat, dims = 2:50, reduction = 'lsi')

DimPlot(seurat, group.by = 'dataset', pt.size = 0.5)


# add gene annotations 
gene_annot <- GetGRangesFromEnsDb(EnsDb.Hsapiens.v86) %>% 
  as_tibble()

gene_annot <- gene_annot %>% 
  mutate(seqnames = paste0("chr", seqnames)) %>% 
  dplyr::filter(seqnames %in% paste0("chr", 1:22)) %>% 
  makeGRangesFromDataFrame(keep.extra.columns = TRUE)

Annotation(seurat) <- gene_annot

# add RPE1 annotation
rpe1_qc <- readRDS(file = "./temp/avo_rpe_qc_stats.rds") %>% 
  left_join(chic_adapters) %>% 
  mutate(cell_id = paste0("CG-ChIC-TAPS-RPE-exp1-", mod, "-", plate, "_", barcode))

rpe1_qc <- rpe1_qc %>% 
  mutate(unique_cuts_log = log10(unique_cuts + 1), 
         n_cpg_log = log10(n_cpg + 1))

rpe1_qc <- left_join(rpe1_qc, readRDS(file = "./temp/avo_rpe_frip.rds"))

rpe_seurat_metadata_old <- seurat@meta.data

rpe_seurat_metadata_new <- rpe_seurat_metadata_old %>% 
  as_tibble(rownames = "cell_id") %>% 
  left_join(rpe1_qc)

seurat@meta.data <- as.data.frame(rpe_seurat_metadata_new)
rownames(seurat@meta.data) <- seurat@meta.data$cell_id

# save object to disk
saveRDS(object = seurat, file = "./temp/avo_rpe1_seurat_all_datasets.rds")



## index  data -------------------------------------------------------

avo_rpe_index_files <- list.files(path = "./meta/index_data/rpe1/", 
                                  pattern = "*.csv", 
                                  full.names = TRUE)

# extract modification and plate
avo_rpe_index_mod <- str_extract(avo_rpe_index_files, "k[0-9]{1,2}me3")
avo_rpe_index_plate <- str_extract(avo_rpe_index_files, "pl[0-9]{1,2}")

avo_rpe_index_data <- avo_rpe_index_files %>% 
  map(.f = read_csv) %>% 
  map(janitor::clean_names)

avo_rpe_index_data <- avo_rpe_index_data %>% 
  map2(.y = avo_rpe_index_mod, .f = function(x, y) {
    x  <- x %>% 
      mutate(mod = y)
  })

avo_rpe_index_data <- avo_rpe_index_data %>% 
  map2(.y = avo_rpe_index_plate, .f = function(x, y) {
    x  <- x %>% 
      mutate(plate = y)
  })

avo_rpe_index_data <- Reduce(f = bind_rows, avo_rpe_index_data)

avo_rpe_index_data <- avo_rpe_index_data %>% 
  mutate(cell = ((match(str_extract(well, pattern = "[A-P]"), LETTERS) -1) * 24) +
           as.integer(str_extract(well, pattern = "[0-9]{1,2}")))
#saveRDS(object = avo_rpe_index_data, file = "./temp/avo_rpe1_index_data_full.rds")

avo_rpe_index_data <- avo_rpe_index_data %>% 
  dplyr::select(mod, plate, cell, fsc, ssc, x488_530_40_2, x561_585_29, x405_460_50)
#saveRDS(object = avo_rpe_index_data, file = "./temp/avo_rpe1_index_data_compressed.rds")



## whole-genome bisulfite sequencing data --------------------------------------

avo_rpe_wgbs <- read_tsv(
  file = "./input/avo/rpe1_bisulfite/context_files_merged/CpG_context_JVL-HMW-RPE-gDNA-BSseq-baseQC-15-02-22_merged_trimmed_bismark_bt2.deduplicated.txt.gz.bismark.cov.gz",
  col_names = c("chr", "start", "end", "meth", "m", "u"))

avo_rpe_wgbs <- avo_rpe_wgbs %>% 
  mutate(cov = m + u) %>% 
  dplyr::select(chr, start, end, meth, cov)

saveRDS(
  object = avo_rpe_wgbs, 
  file = "./temp/avo_rpe_wgbs.rds")



## scEpi-seq data --------------------------------------------------------------

# list input files
avo_rpe_files <- list.files(
  path = "./input/avo/z_methylation/rpe1/", 
  pattern = "*.tsv", 
  full.names = TRUE)

# load raw data
avo_rpe_data <- avo_rpe_files %>% 
  as.list %>% 
  map(.f = parse_scmo_zmeth)


# convert all chromosome columns to character and remove 'gene' column
avo_rpe_data <- avo_rpe_data %>% 
  map(.f = prefix_chromosome_names) %>% 
  map(.f = scale_beta)

# extract cell IDs and the modification
avo_rpe_data <- avo_rpe_data %>%
  map(.f = function(x){
    x <- x %>% 
      mutate(plate = extract_plate_id(read_id), 
             cell = extract_cell_id(read_id), 
             mod = extract_modification(read_id))
    return(x)
  })

# combine into single tibble
avo_rpe_data <- Reduce(f = bind_rows, x = avo_rpe_data)
#saveRDS(object = avo_rpe_data, file = "./temp/avo_rpe_full_data.rds")

# extract QC stats
avo_rpe_qc <- avo_rpe_data %>% 
  group_by(plate, cell, mod) %>% 
  dplyr::summarise(n_cpg = n(),
                   avg_beta = mean(meth),
                   unique_cuts = length(unique(read_id)),
                   ta_frac = sum(cutsite == "TA", na.rm = TRUE)/n(),
                   mod = sample(mod, size = 1)) %>%
  ungroup


# add metadata for QC cutoffs
qc_thresholds_rpe <- readxl::read_xlsx("./meta/qc_cutoffs_rpe.xlsx") %>% 
  mutate(f_min = f_min * 100, 
         f_max = f_max * 100)

avo_rpe_qc <- left_join(avo_rpe_qc, qc_thresholds_rpe)

avo_rpe_qc <- avo_rpe_qc %>% 
  mutate(passed_qc = ifelse(
    log10(unique_cuts) > n_min_log10 & avg_beta > f_min & avg_beta < f_max,
    "passed",
    "failed"))

#saveRDS(object = avo_rpe_qc, file = "./temp/avo_rpe_qc_stats.rds")

filterlim <- tibble(mod = c("k27me3","k36me3", "k9me3"),
                    n_min_log10 = c(3.8, 3.8, 3.8),
                    n_max_log10 = rep(log10(max(avo_rpe_qc$unique_cuts)), 3),
                    f1 = c(0, 0.65, 0),
                    f2 = c(0.45, max(avo_rpe_qc$avg_beta), 0.45), 
                    f = c(0.45, 0.65, 0.45))

avo_rpe_qc %>%
  ggplot(aes(log10(unique_cuts), avg_beta, col = mod)) +
  geom_point(size = 0.8, alpha = ifelse(avo_rpe_qc$passed_qc == "passed", 1, 0.1)) +
  facet_grid(cols = vars(mod)) +
  scale_color_manual(values = color_scheme) +
  theme_bw(base_size = 14)

avo_rpe_qc_passed <- avo_rpe_qc %>% 
  dplyr::select(cell, mod, plate, passed_qc)

# create filtered version of the QC filtered version of the dataset
avo_rpe_data <- left_join(avo_rpe_data, avo_rpe_qc_passed)

avo_rpe_data_qc_filtered <- avo_rpe_data %>% 
  filter(passed_qc == "passed") %>% 
  dplyr::select(read_id, cell, chr, start, meth, cutsite, mod, plate)

avo_rpe_data_qc_filtered <- avo_rpe_data_qc_filtered %>% 
  add_column(end = avo_rpe_data_qc_filtered$start + 1, .after = 4)

#saveRDS(object = avo_rpe_data_qc_filtered, file = "./temp/avo_rpe_qc_filtered.rds")



## FRIP statistics -------------------------------------------------------------

rpe_frip_files <- list.files(path = "./input/avo/rpe1_frip/", 
                             pattern = ".txt", 
                             full.names = TRUE)

rpe_frip_data <- rpe_frip_files %>% 
  map(.f = parse_rpe_frip) %>% 
  Reduce(f = bind_rows)

saveRDS(object = rpe_frip_data, file = "./temp/avo_rpe_frip.rds")

avo_rpe_qc <- left_join(avo_rpe_qc, rpe_frip_data)


## MACS output: CHiC peaks -----------------------------------------------------

avo_rpe_peaks_k27 <- parse_broad_peaks(
  file = "./input/avo/macs3_peaks/rpe1/k27me3_broad/macs3_peaks.broadPeak",
  returnGranges = FALSE)

avo_rpe_peaks_k36 <- parse_broad_peaks(
  file = "./input/avo/macs3_peaks/rpe1/k36me3_broad/macs3_peaks.broadPeak",
  returnGranges = FALSE)

avo_rpe_peaks_k9 <- parse_broad_peaks(
  file = "./input/avo/macs3_peaks/rpe1/k9me3_broad/macs3_peaks.broadPeak",
  returnGranges = FALSE)

avo_rpe_peaks_k27 <- avo_rpe_peaks_k27 %>% 
  mutate(mod = "k27me3")

avo_rpe_peaks_k36 <- avo_rpe_peaks_k36 %>% 
  mutate(mod = "k36me3")

avo_rpe_peaks_k9 <- avo_rpe_peaks_k9 %>% 
  mutate(mod = "k9me3")

avo_rpe_peaks <- bind_rows(avo_rpe_peaks_k27, avo_rpe_peaks_k36, avo_rpe_peaks_k9)

avo_rpe_peaks <- avo_rpe_peaks %>% 
  prefix_chromosome_names()

#saveRDS(object = avo_rpe_peaks, file = "./temp/avo_rpe_broad_peaks_combined.rds")
rm(avo_rpe_peaks_k27, avo_rpe_peaks_k36, avo_rpe_peaks_k9)
rm(avo_rpe_peaks)




# K562 // Zeller et al. 2022 (ChIC only) ---------------------------------------

# Files from 2022 Nature Genetics publication
zeller_k562_single_cell_files <- c("./input/zeller/zeller2022/GSM5018608_K562-EtOH-H3K9me3.G1sorted.merged.sorted.tagged.countTable.csv.gz", 
                                   "./input/zeller/zeller2022/GSM5018610_K562-EtOH-H3K27me3.G1sorted.merged.sorted.tagged.countTable.csv.gz")

# In addition, there is bulk data from a H3K36me3 titration experiment
zeller_k562_titration_file <- "./input/zeller/PZ-K562-sap-mAB-H3K36me3-1in1600_10kb.bw"

# create emtpy list
zeller_chic_data <- list()

# first, import bulk titration bigwig file and wrangle data
zeller_chic_data[["k36me3"]] <- rtracklayer::import.bw(con = zeller_k562_titration_file) %>% 
  as.data.frame() %>% 
  as_tibble() %>% 
  dplyr::select(seqnames, start, end, score) %>% 
  rename(seqnames = "chr") %>% 
  prefix_chromosome_names

# aggregate the data into 50 kb bins (resolution of other data)
zeller_chic_data[["k36me3"]] <- zeller_chic_data[["k36me3"]] %>% 
  group_by(chr) %>% 
  mutate(start = floor(start/50000)*50000, 
         end = start + 50000) %>% 
  group_by(chr, start, end) %>% 
  summarise(count = sum(score)) %>% 
  ungroup


# read data
zeller_chic_data[["k9me3"]] <- read_zeller_k562(zeller_k562_single_cell_files[1])
zeller_chic_data[["k27me3"]] <- read_zeller_k562(zeller_k562_single_cell_files[2])

# aggregate data
zeller_chic_data[["k9me3"]] <- zeller_chic_data[["k9me3"]] %>% 
  prefix_chromosome_names() %>% 
  mutate(count = apply(zeller_chic_data[["k9me3"]][, 4:1155], 1, sum, na.rm = TRUE)) %>% 
  dplyr::select(chr, start, end, count)

zeller_chic_data[["k27me3"]] <- zeller_chic_data[["k27me3"]] %>% 
  prefix_chromosome_names() %>% 
  mutate(count = apply(zeller_chic_data[["k27me3"]][, 4:1155], 1, sum, na.rm = TRUE)) %>% 
  dplyr::select(chr, start, end, count)

#saveRDS(object = zeller_chic_data, file = "./temp/zeller_k562_chic_counts_binned.rds")





