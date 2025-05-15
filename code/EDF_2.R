library(tidyverse)
#library(encodeChIPqc)
library(GenomicRanges)
library(biomaRt)

library(GenomicAlignments)
library(GenomicFeatures)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(BSgenome.Hsapiens.UCSC.hg38)

source("./scripts/functions.R")

color_scheme <- c("#D674BA","#86BF88","#7A416A")
color_scheme_2 <- c("#0D3B66", "#FAF0CA", "#F95738", "#6C969D", "#2E5339")


library_anno <- readxl::read_xlsx(path = "./meta/library_annotation.xlsx")



# AvO K562: Quality control (QC) -----------------------------------------------


## Insert size distributions ---------------------------------------------------

insert_size_data <- readRDS(file = "./temp/avo_k562_insert_sizes.rds")

insert_size_data_median <- insert_size_data %>% 
  mutate(insert_size = abs(insert_size)) %>% 
  group_by(mod) %>% 
  summarise(median_insert = median(insert_size)) %>% 
  mutate(mod = paste0("H3", str_to_title(mod)))

insert_size_data %>% 
  mutate(mod = paste0("H3", str_to_title(mod))) %>% 
  ggplot(aes(abs(insert_size), fill = mod)) +
  geom_histogram(bins = 50, alpha = 0.5) +
  geom_vline(aes(xintercept = median_insert, col = mod), 
             data = insert_size_data_median, 
             lty = 2, lwd = 1) +
  facet_grid(facets = vars(mod), scales = "free_y") +
  theme_bw(base_size = 14) +
  scale_fill_manual(values = color_scheme) +
  scale_color_manual(values = color_scheme) +
  theme(legend.position = "none") +
  labs(
    title = "K562 sequencing library insert size distributions", 
    x = "Insert size of 100,000 sampled reads (bp)", 
    y = "Count")



## import data -------

avo_k562_qc <- readRDS(file = "./temp/avo_k562_qc_stats.rds")
frip_k562 <- readRDS(file = "./temp/avo_k562_frip_statistics.rds")

avo_k562_qc <- left_join(avo_k562_qc, frip_k562)

#k562_layout <- read_csv(file = "./meta/sort_layouts/20200303_k562.csv")
avo_k562_peaks <- readRDS(file = "./temp/avo_k562_broad_peaks_combined.rds")


# add information about empty controls (indices 357::360 and 381::384)
avo_k562_qc <- avo_k562_qc %>%
  mutate(ctrl = ifelse(
    plate == "pl1",
    ifelse(cell %in% c(357:360, 381:384), "empty", "cell"), 
    ifelse(cell %in% c(358:360, 382:384), "empty", "cell")))


## Scatter plot // QC selection // ChIC counts vs. average 5mC -----------------

avo_k562_qc %>% 
  mutate(mod = paste0("H3", str_to_title(mod))) %>%
  ggplot(aes(log10(unique_cuts + 1), avg_beta, col = mod)) +
  facet_grid(cols = vars(mod)) +
  geom_point(size = 1) +
  geom_hline(aes(yintercept = f_min), lty = 2) +
  geom_vline(aes(xintercept = n_min_log10), lty = 2) +
  theme_bw(base_size = 14) +
  scale_color_manual(values = color_scheme) +
  theme(legend.position = "none") +
  labs(x = "Unique ChIC counts (log10)", y = "Mean methylation (per cell)")
#ggsave(filename = "./plots/avo_k562_qc_log_chic_vs_avg_meth.pdf", width = 10, height = 4)


## Histogram of unique CpGs for cells which pass QC vs. those which don't ------
avo_k562_qc %>% 
  mutate(mod = paste0("H3", str_to_title(mod))) %>% 
  mutate(passed_qc = paste0("QC ", qc)) %>% 
  ggplot(aes(log10(n_cpg + 1), fill = mod)) +
  geom_density(alpha = 0.5) +
  facet_grid(rows = vars(passed_qc), cols = vars(mod)) +
  scale_fill_manual(values = color_scheme) +
  theme_bw(base_size = 14) +
  theme(legend.position = "none") +
  labs(x = "Detected CpGs (log10)", y = "Density")
#ggsave(filename = "./plots/avo_k562_qc_log_numberCpgs_passed_vs_failed.pdf", width = 10, height = 4)

# unique reads per cell stratified for QC status
avo_k562_qc %>% 
  mutate(mod = paste0("H3", str_to_title(mod))) %>% 
  mutate(passed_qc = paste0("QC ", qc)) %>% 
  ggplot(aes(log10(unique_cuts + 1), fill = mod)) +
  geom_density(alpha = 0.5) +
  facet_grid(rows = vars(passed_qc), cols = vars(mod)) +
  scale_fill_manual(values = color_scheme) +
  theme_bw(base_size = 14) +
  theme(legend.position = "none") +
  labs(x = "Unique cut sites (log10)", y = "Density")

avo_k562_qc %>% 
  #filter(qc == "passed") %>% 
  group_by(qc, mod) %>% 
  dplyr::summarise(
    mean_cpg = mean(n_cpg),
    median_cpg = median(n_cpg),
    mean_cuts = mean(unique_cuts),
    median_cuts = median(unique_cuts))
# 
# qc     mod    mean_cpg median_cpg mean_cuts median_cuts
# <chr>  <chr>     <dbl>      <int>     <dbl>       <int>
#   1 failed k27me3    8297.        514     3105.         229
# 2 failed k36me3   12613.       3493     4113.        1212
# 3 failed k9me3     9284.       1743     4085.         796
# 4 passed k27me3   73737.      54602    26754.       19793
# 5 passed k36me3   73179.      65436    23936.       21215
# 6 passed k9me3    87516.      82392    39642.       37734

# violin plot of FRIP across histone mods
avo_k562_qc %>% 
  #filter(qc == "passed") %>% 
  mutate(mod = paste0("H3", str_to_title(mod))) %>% 
  ggplot(aes(mod, frip)) +
  geom_violin(aes(fill = mod), alpha = 0.7) +
  geom_boxplot(width = 0.1, outliers = FALSE, ) +
  theme_bw(base_size = 18) +
  scale_fill_manual(values = color_scheme) +
  labs(x = NULL, y = "Fraction reads in peaks (FRIP)") +
  theme(legend.position = "none")
#ggsave(filename = "./plots/avo_k562_frip.pdf")


# no. of cells passing qc
avo_k562_qc %>% 
  group_by(mod) %>% 
  summarise(n_passed = sum(qc == "passed"), 
            n_total = n()) %>% 
  mutate(perc_passed = n_passed/n_total*100)


# mod    n_passed n_total perc_passed
# <chr>     <int>   <int>       <dbl>
#   1 k27me3     1495    1920        77.9
# 2 k36me3      255     384        66.4
# 3 k9me3       231     384        60.2


mean(avo_k562_qc$frip[avo_k562_qc$passed_qc == "passed"]) # 0.6638748
mean(avo_k562_qc$frip) # 0.6638748

avo_k562_qc %>% 
  group_by(qc == "passed") %>% 
  group_by(mod) %>% 
  summarise(mean_frip = mean(frip), 
            median_frip = median(frip), 
            lower_bound = min(frip), 
            upper_bound = max(frip))

avo_k562_qc %>% 
  filter(qc == "passed") %>% 
  group_by(mod) %>% 
  summarise(mean_frip = mean(frip), 
            median_frip = median(frip), 
            lower_bound = min(frip), 
            upper_bound = max(frip))


# Average No. of CpGs per read per modification
avo_k562_qc %>% 
  filter(qc == "passed") %>% 
  mutate(mod = paste0("H3", str_to_title(mod))) %>% 
  mutate(cpg_per_cut = n_cpg/unique_cuts) %>% 
  filter(cpg_per_cut < 3.5) %>% 
  ggplot(aes(cpg_per_cut, fill = plate)) +
  geom_histogram(bins = 50, alpha = 0.7) +
  theme_bw(base_size = 14) +
  labs(x = "CpGs per Read", y = "Count") +
  facet_grid(rows = vars(mod))
#scale_fill_manual(values = color_scheme) +
#theme(legend.position = "none")


## empty wells / negative controls (R3.6) -----

with(avo_k562_qc, median(unique_cuts[label == 'cell' & passed_qc == "passed"]))/
  with(avo_k562_qc, median(unique_cuts[label == 'empty']))

with(avo_k562_qc, median(unique_cuts[label == 'cell']))/
  with(avo_k562_qc, median(unique_cuts[label == 'empty']))


avo_k562_qc %>% 
  mutate(mod = paste0("H3", str_to_title(mod))) %>% 
  ggplot(aes(ctrl, unique_cuts, fill = mod)) +
  geom_boxplot(width = 0.5) +
  scale_y_log10() +
  facet_grid(rows = vars(mod)) +
  theme_bw(base_size = 14) +
  theme(legend.position = "none") +
  scale_fill_manual(values = color_scheme) +
  labs(x = NULL, y = "Unique ChiC cuts (log10)")

with(avo_k562_qc, table(ctrl, qc))



# Genome coverage as measured by ChIC peaks // % detected CpGs ----------------

# read chromosome annotation from reference genome ('samtools faidx')
genome_baseline <- read_table(file = "./meta/genomes/Homo_sapiens.GRCh38.dna.primary_assembly.fa.fai", 
                              col_names = c("chr", "size", "x", "y", "z")) %>% 
  dplyr::select(chr, size) %>% 
  rename(size = "chr_size") %>% 
  prefix_chromosome_names

# calculate the total size of chic peaks identified by MACS
genome_covered_chic <- avo_k562_peaks %>% 
  dplyr::select(chr, start, end) %>% 
  mutate(ws = abs(start - end)) %>% 
  group_by(chr) %>% 
  summarise(total_ws = sum(ws)) %>% 
  ungroup %>% 
  dplyr::select(chr, total_ws)
rm(avo_k562_peaks)

# calculate the no. of CpGs in the human genome
cpg_positions <- lapply(names(Hsapiens)[1:24], function(x) start(matchPattern("CG", Hsapiens[[x]])))

# need to multiply by two to account for both strands
cpg_anno <- tibble(chr = names(Hsapiens)[1:24], 
                   n_cpg_strand_agnostic = sapply(cpg_positions, length), 
                   n_cpg = n_cpg_strand_agnostic*2)
sum(test$cov==0)

# calculate no. of CpGs with non-zero coverage in ENCODE WGBS data
k562_wgbs <- readRDS("./temp/encode_K562_wgbs_ENCFF721JMB.rds")

k562_wgbs_coverage <- k562_wgbs %>% 
  filter(cov != 0) %>% 
  group_by(chr) %>% 
  summarise(cov_wgbs = n())
rm(k562_wgbs)

# calculate no of CpGs with coverage of at least one read in Episeq data
k562_episeq <- readRDS("./temp/avo_k562_qc_filtered.rds")

k562_episeq_cov <- k562_episeq %>% 
  group_by(chr) %>% 
  summarise(cov_cpg = length(unique(start))) %>% 
  ungroup

k562_episeq_cov_strand_agnostic <- k562_episeq %>% 
  mutate(strand = str_sub(read_id, -1)) %>% 
  mutate(start = ifelse(strand == "-", start, start - 1), 
         end = ifelse(strand == "-", end, end - 1)) %>% 
  dplyr::select(chr, start) %>% 
  group_by(chr) %>% 
  summarise(cov_cpg_strand_agnostic = length(unique(start))) %>% 
  ungroup
rm(k562_episeq)

# which chromosomes to look at 
chrs <- c(paste0("chr", 1:22), "chrX")

# aggregate data
genome_coverage <- left_join(genome_baseline, cpg_anno) %>% 
  left_join(k562_wgbs_coverage) %>% 
  left_join(genome_covered_chic) %>% 
  left_join(k562_episeq_cov) %>% 
  left_join(k562_episeq_cov_strand_agnostic)

genome_coverage <- genome_coverage %>% 
  mutate(cov_rel_peaks = total_ws/chr_size, 
         cov_rel_cpg = cov_cpg/n_cpg, 
         cov_rel_cpg_strand_agnostic = cov_cpg_strand_agnostic/n_cpg_strand_agnostic) %>% 
  filter(chr %in% chrs)

genome_coverage %>% 
  pivot_longer(cols = starts_with("cov_rel_cpg")) %>% 
  ggplot(aes(factor(chr, levels = chrs), value*100)) +
  geom_col(fill = "chocolate", alpha = 0.6) +
  facet_grid(rows = vars(name)) +
  theme(axis.text.x = element_text(angle = 90)) +
  labs(x = "Chromosome", y = "Coverage (%)")

genome_coverage %>% 
  summarise(all_cpgs = sum(n_cpg), 
            all_cpgs_agnostic = sum(n_cpg_strand_agnostic), 
            n_cpgs = sum(cov_cpg), 
            n_cpgs_agnostic = sum(cov_cpg_strand_agnostic)) %>% 
  mutate(total_cov = n_cpgs/all_cpgs*100, 
         total_cov_agnostic = n_cpgs_agnostic/all_cpgs_agnostic*100)

sum(genome_coverage$cov_cpg)/sum(genome_coverage$n_cpg)



## Coverage per CpG per single cell (R3.2) -------------------------------------


avo_k562_data_qc_filtered_cov <- 
  avo_k562_data_qc_filtered %>% 
  group_by(mod, cell, chr, start) %>% 
  dplyr::summarise(cov = n(), 
                   beta = mean(meth)) %>% 
  ungroup

sum(avo_k562_data_qc_filtered_cov$beta %in% c(0, 50, 100))/nrow(avo_k562_data_qc_filtered_cov)

avo_k562_data_qc_filtered_cov %>% 
  #filter(cov > 1) %>% 
  ggplot(aes(beta)) +
  geom_histogram() +
  scale_y_log10()

avo_k562_data_qc_filtered %>% 
  filter(cell == 4) %>% 
  filter(chr == "chr9") %>% 
  filter(start == 130937266)

sum(avo_k562_data_qc_filtered_cov$beta %in% c(0, 100))/nrow(avo_k562_data_qc_filtered_cov)
table(avo_k562_data_qc_filtered_cov$beta)



# Correlation between ChiC-TAPS and ENCODE WGBS --------------------------------

# import K562 Epi-seq data (QC-filtered)
k562_episeq <- readRDS("./temp/avo_k562_qc_filtered.rds")

# bin data (10 kb)
k562_episeq_binned <- k562_episeq %>% 
  scale_beta %>% 
  mutate(cov = 1) %>% 
  bin_methylation(bin_size = 10000)

# collapse CpGs across modifications and cells
k562_episeq <- k562_episeq %>% 
  mutate(cov = 1) %>% 
  group_by(chr, start, end) %>% 
  summarise(cov = sum(cov), 
            meth = mean(meth))
#saveRDS(k562_episeq, file = "./temp/avo_k562_qc_filtered_summarised_cpg.rds")  

# Read ENCODE K562 data (single-CpG resolution)
encode_wgbs_k562 <- as.list(c("./temp/wgbs_encode_K562_ENCFF721JMB.rds",
                              "./temp/wgbs_encode_K562_ENCFF867JRG.rds")) %>% 
  map(.f = readRDS)

# account for 0-based encoding
encode_wgbs_k562 <- encode_wgbs_k562 %>% 
  map(.f = function(x){
    x <- x %>% 
      mutate(start = start + 1, 
             end = end + 1)
  })

# single-CpG resolution

# simulate low-coverage data from ENCODE
alpha = 0.875 # average conversion rate across plates
beta = 0.0023 # false positive

simulated_episeq <- encode_wgbs_k562[[1]] %>% 
  filter(cov > 9) %>% 
  sample_n(size = 1000000) %>% 
  scale_beta(scale_to_100 = FALSE)

simulated_episeq <- simulated_episeq %>% 
  rowwise() %>% 
  mutate(cov_sim = sample(k562_episeq$cov, size = 1), 
         meth_sim = rbinom(n = sample(cov_sim, size = 1), size = 1, prob = meth) %>% 
           simulate_taps_detection(beta = beta, alpha = alpha) %>% 
           mean)

simulated_episeq %>% 
  filter(cov_sim > 0) %>% 
  with(., cor(meth, meth_sim, method = "pearson"))


# combine single-CpG Epi-seq and WGBS data
k562_episeq <- left_join(k562_episeq, encode_wgbs_k562[[1]], by = c("chr", "start", "end")) %>% 
  filter(cov.y > 9)

k562_episeq %>% 
  sample_n(size = 1000000) %>% 
  filter(cov.x > 0) %>% 
  with(., cor(meth.x, meth.y, method = "pearson"))


# for different cutoffs w.r.t. to Episeq coverage, determine actual and simulated correlation
cutoff <- 0:5
cor_simulated <- vector(mode = "numeric", length = 5L)
cor_actual <- vector(mode = "numeric", length = 5L)

for (i in 1:length(cutoffs)){
  print(i)
  cor_simulated[i] <- simulated_episeq %>% 
    filter(cov_sim > cutoff[i]) %>% 
    with(., cor(meth, meth_sim, method = "pearson", use = "pairwise"))
  cor_actual[i] <- k562_episeq %>% 
    sample_n(size = 1000000) %>% 
    filter(cov.x > cutoff[i]) %>% 
    with(., cor(meth.x, meth.y, method = "pearson", use = "pairwise"))
}

# collect data
cor_episeq_wgbs_single_cpg <- tibble(cutoff, cor_simulated, cor_actual)

cor_episeq_wgbs_single_cpg %>% 
  pivot_longer(cols = 2:3) %>% 
  mutate(name = str_replace(name, pattern = "cor_", "")) %>% 
  ggplot(aes(cutoff, value, col = name)) +
  geom_point(size = 2) +
  geom_line(lty = 2) +
  lims(y = c(0, 1)) +
  theme_bw(base_size = 14) +
  theme(legend.position = "inside", legend.position.inside = c(0.15, 0.1), legend.title = element_blank()) +
  labs(x = "Cutoff coverage (n reads)", y = "Correlation (Pearson's r)")


# Read ENCODE data across 4 different cell types (10 kb bins)
encode_wgbs_binned <- list.files(path = "./temp", pattern = "wgbs_encode", full.names = TRUE) %>% 
  as.list %>% 
  map(.f = function(x){
    cat(sprintf("Reading %s \n", x))
    binned_meth <- readRDS(x) %>% 
      bin_methylation(bin_size = 10000) %>% 
      mutate(cell_line = str_extract(string = x, pattern = "(GM12878|H1|HepG2|K562)"), 
             file = str_extract(string = x, pattern = "ENC[0-9A-Z]*"))
  })

# compare binned methylation values

encode_wgbs_binned <- encode_wgbs_binned %>% 
  Reduce(f = bind_rows)

comparison_episeq_wgbs_binned <- left_join(encode_wgbs_binned, k562_episeq_binned, by = c("chr", "bin")) %>% 
  drop_na() %>% 
  filter(chr %in% paste0("chr", c(1:22, "X")))

# convert file names (of BED) to "replicate"
comparison_episeq_wgbs_binned <- comparison_episeq_wgbs_binned %>% 
  mutate(replicate = ifelse(file %in% c("ENCFF570TIL", "ENCFF434CNG", "ENCFF453UDK", "ENCFF721JMB"), "replicate1", "replicate2"))

comparison_episeq_wgbs_binned %>%
  ggplot()+
  geom_bin2d(aes(x = meth.x, y = meth.y), bins=100)+
  geom_smooth(aes(x = meth.x, y = meth.y), method = 'lm', linetype = 2, col = "green1") +
  coord_fixed() +
  theme_bw() +
  scale_fill_viridis_c(option = "C", trans='log10') +
  ylab('TET-assisted Pyrdine Borane sequencing [TAPS] ') +
  xlab('Whole Genome Bisulfite Sequencing [WGBS]') +
  facet_grid(rows = vars(cell_line), cols = vars(replicate)) +
  lims(x = c(0, 100), y = c(0, 100))
#ggsave(filename = "./plots/comparison_betas_10kb_episeq_wgbs_k562.pdf", height = 10, width = 9)

comparison_episeq_wgbs_binned %>% 
  group_by(cell_line, replicate) %>% 
  summarise(cor = cor(meth.x, meth.y, use = "pairwise", method = "pearson"))



# Concordance of single CpGs between cells for Epi-seq (R3.8) ------------------

k562_episeq <- readRDS("./temp/avo_k562_qc_filtered.rds")

k562_episeq_concordance <- k562_episeq %>% 
  filter(meth %in% c(0, 100)) %>% 
  mutate(meth = meth/100) %>% 
  group_split(mod)

k562_episeq_concordance <- k562_episeq_concordance %>% 
  map(.f = function(x){
    y = x %>%
      group_by(chr, start) %>% 
      mutate(n_obs = n(), 
             n_cells = length(unique(cell))) %>% 
      filter(n_obs > 1)
    return(y)})

k562_episeq_concordance <- k562_episeq_concordance %>% 
  map(.f = function(x){
    y = x %>% 
      summarise(avg_meth = mean(meth), 
                mod = sample(mod, size = 1), 
                n_obs = sample(n_obs, size = 1), 
                n_cells = sample(n_cells, size = 1), 
                pairwise = concordance_pairwise(meth), 
                majority = concordance_majority(meth)) %>% 
      ungroup
    return(y)
  })

# combine across modifications
k562_episeq_concordance <- Reduce(x = k562_episeq_concordance, f = bind_rows)
#saveRDS(object = k562_episeq_concordance, file = "./temp/avo_k562_concordance_per_cpg.rds")


# Distribution by histone marks
k562_episeq_concordance %>%
  group_by(mod) %>% 
  filter(pairwise < 0.99) %>% 
  slice_sample(n = 1000000) %>% 
  mutate(mod = paste0("H3", str_to_title(mod))) %>% 
  ggplot(aes(pairwise, fill = mod)) +
  geom_histogram() +
  theme_bw(base_size = 12) +
  theme(legend.position = "none") +
  scale_fill_manual(values = color_scheme) +
  labs(x = NULL, y = "Concordance (%)") +
  facet_grid(rows = vars(mod), scales = "free_y")


# plot concordance per histone mark
k562_episeq_concordance %>%
  group_by(mod) %>% 
  mutate(mod = paste0("H3", str_to_title(mod))) %>% 
  summarise(fully_homogeneous = sum(pairwise == 1)/n()*100) %>% 
  ggplot(aes(mod, fully_homogeneous, fill = mod)) +
  geom_col() +
  theme_bw(base_size = 12) +
  theme(legend.position = "none") +
  scale_fill_manual(values = color_scheme) +
  labs(x = NULL, y = "Concordance (%)") +
  lims(y = c(0, 100))
#ggsave(filename = "./plots/cpg_concordance_k562_episeq.pdf", width = 5, height = 7)


# now, compare to intra-CpG concordance in ENCODE data overlapping ChiC peaks

# read AvO K562 ChIC peaks
avo_k562_peaks <- readRDS(file = "./temp/avo_k562_broad_peaks_combined.rds")
avo_k562_peaks_granges <- makeGRangesFromDataFrame(avo_k562_peaks)

# read ENCODE WGBS data for K562 (2 replicates)
encode_wgbs <- c("./temp/k562-encode-wgbs-1-ENCFF721JMB.rds", 
                 "./temp/k562-encode-wgbs-2-ENCFF867JRG.rds") %>% 
  map(readRDS)

# find overlaps between data
encode_wgbs_overlaps <- encode_wgbs %>% 
  map(.f = makeGRangesFromDataFrame) %>% 
  map(.f = findOverlaps,  subject = avo_k562_peaks_granges)

concordance_encode_wgbs <- list()

for(i in 1:length(encode_wgbs)){
  concordance_encode_wgbs[[i]] <- bind_cols(avo_k562_peaks[encode_wgbs_overlaps[[i]]@to, ], 
                                            encode_wgbs[[i]][encode_wgbs_overlaps[[i]]@from, ])
}

# clean up
rm(avo_k562_peaks, avo_k562_peaks_granges, encode_wgbs, encode_wgbs_overlaps)

# update colnames
new_colnames <- c("chr_peak", "start_peak", "end_peak", "strand_peak", "chr_wgbs", "start_wgbs", "end_wgbs", "strand_wgbs")
colnames(concordance_encode_wgbs[[1]])[c(1, 2, 3, 6, 12, 13, 14, 15)] <- new_colnames
colnames(concordance_encode_wgbs[[2]])[c(1, 2, 3, 6, 12, 13, 14, 15)] <- new_colnames

# convert to tibble
concordance_encode_wgbs[[1]] <- as_tibble(concordance_encode_wgbs[[1]])
concordance_encode_wgbs[[2]] <- as_tibble(concordance_encode_wgbs[[2]])

# calculate concordance - here we only use pairwise concordance b/c it supports vector input
concordance_encode_wgbs <- concordance_encode_wgbs %>% 
  map(.f = function(x){
    y = x %>% 
      dplyr::select(chr_wgbs, start_wgbs, mod, name, cov, beta) %>% 
      mutate(beta = beta/100) %>% 
      mutate(pairwise = beta^2 + (1 - beta)^2)
  })

# add information about replicate / file
concordance_encode_wgbs[[1]] <- concordance_encode_wgbs[[1]] %>% 
  mutate(file = "ENCFF721JMB")
concordance_encode_wgbs[[2]] <- concordance_encode_wgbs[[2]] %>% 
  mutate(file = "ENCFF867JRG")

# merge
concordance_encode_wgbs <- Reduce(bind_rows, x = concordance_encode_wgbs)

# plot
concordance_encode_wgbs %>% 
  group_by(file, mod) %>% 
  mutate(mod = paste0("H3", str_to_title(mod))) %>% 
  summarise(fully_homogeneous = sum(pairwise == 1)/n()*100) %>% 
  ggplot(aes(mod, fully_homogeneous, fill = mod)) +
  geom_col() +
  theme_bw(base_size = 12) +
  theme(legend.position = "none") +
  scale_fill_manual(values = color_scheme) +
  labs(x = NULL, y = "Concordance (%)") +
  lims(y = c(0, 100)) +
  facet_grid(cols = vars(file))
#ggsave(filename = "./plots/cpg_concordance_k562_encode.pdf", width = 9, height = 7)

# clean up
rm(concordance_encode_wgbs, k562_episeq_concordance, k562_episeq)



# ENCODE DNA methylation w.r.t. to ChIC-seq peaks (R3.7) -----------------------

# read AvO K562 ChIC peaks
avo_k562_peaks <- readRDS(file = "./temp/avo_k562_broad_peaks_combined.rds")

# read ENCODE WGBS data for K562 (2 replicates)
encode_wgbs <- list()
encode_wgbs[[1]] <- readRDS(file = "./temp/encode_K562_wgbs_ENCFF721JMB.rds")
encode_wgbs[[2]] <- readRDS(file = "./temp/encode_K562_wgbs_ENCFF867JRG.rds")

# find overlaps with ChiC peaks
encode_wgbs_overlaps <- encode_wgbs %>% 
  map(.f = makeGRangesFromDataFrame)

encode_wgbs_overlaps <- encode_wgbs_overlaps %>% 
  map(.f = findOverlaps,  subject = makeGRangesFromDataFrame(avo_k562_peaks))

encode_wgbs_in_peaks <- list()
encode_wgbs_in_peaks[[1]] <- bind_cols(avo_k562_peaks[encode_wgbs_overlaps[[1]]@to, ], 
                                       encode_wgbs[[1]][encode_wgbs_overlaps[[1]]@from, ])
encode_wgbs_in_peaks[[2]] <- bind_cols(avo_k562_peaks[encode_wgbs_overlaps[[2]]@to, ], 
                                       encode_wgbs[[2]][encode_wgbs_overlaps[[2]]@from, ])

new_colnames <- c("chr_peak", "start_peak", "end_peak", "strand_peak", "chr_wgbs", "start_wgbs", "end_wgbs", "strand_wgbs")
colnames(encode_wgbs_in_peaks[[1]])[c(1, 2, 3, 6, 12, 13, 14, 15)] <- new_colnames
colnames(encode_wgbs_in_peaks[[2]])[c(1, 2, 3, 6, 12, 13, 14, 15)] <- new_colnames
rm(new_colnames)

# avg. beta per bin
encode_wgbs_in_peaks <- encode_wgbs_in_peaks %>% 
  map(.f = function(x){
    y = x %>% 
      group_by(name, mod) %>% 
      dplyr::summarise(avg_beta = mean(meth), 
                       no_cpgs = n(), 
                       total_cov = sum(cov), 
                       cov_per_cpg = total_cov/no_cpgs) %>% 
      ungroup()
    return(y)
  })

# add file name, condense into 1 tibble
encode_wgbs_in_peaks[[1]] <- encode_wgbs_in_peaks[[1]] %>% 
  mutate(file = "ENCFF721JMB")

encode_wgbs_in_peaks[[2]] <- encode_wgbs_in_peaks[[2]] %>% 
  mutate(file = "ENCFF867JRG")

encode_wgbs_in_peaks <- Reduce(f = bind_rows, x = encode_wgbs_in_peaks)

# plot
encode_wgbs_in_peaks %>% 
  filter(total_cov > 50) %>% 
  filter(no_cpgs > 10) %>% 
  mutate(mod = paste0("H3", str_to_title(mod))) %>% 
  ggplot(aes(avg_beta, fill = mod)) +
  geom_histogram(bins = 100) +
  facet_grid(rows = vars(mod), cols = vars(file), scales = "free_y") +
  theme_bw(base_size = 12) +
  theme(legend.position = "none") +
  labs(x = "DNA methylation (%)", y = "Count") +
  scale_fill_manual(values = color_scheme)
ggsave(filename = "./plots/encode_k562_avg_beta_in_chic_peaks.pdf", width = 6, height = 7)


# get statistics (mean methylation per mark)
encode_wgbs_in_peaks %>% 
  filter(total_cov > 50) %>% 
  filter(no_cpgs > 10) %>% 
  group_by(file, mod) %>% 
  summarise(mean_beta = mean(avg_beta))

# 1 ENCFF721JMB k27me3     14.2 
# 2 ENCFF721JMB k36me3     47.4 
# 3 ENCFF721JMB k9me3       8.98
# 4 ENCFF867JRG k27me3     13.8 
# 5 ENCFF867JRG k36me3     47.2 
# 6 ENCFF867JRG k9me3       8.59

# save temporary file
saveRDS(object = encode_wgbs_in_peaks, file = "./temp/encode_k562_wgbs_in_peaks.rds")

# clean up
rm(avo_k562_peaks, encode_wgbs, encode_wgbs_overlaps)
rm(encode_wgbs_in_peaks)



# Correlation with Zeller et al. ChiC and ENCODE Chip-seq ----------------------

# load data 
chip_data <- readRDS(file = "./temp/encode_k562_chip_data_aggregated.rds")
chip_controls <- readRDS(file = "./temp/encode_k562_chip_controls_aggregated.rds")
zeller_chic_data <- readRDS(file = "./temp/zeller_k562_chic_counts_binned.rds")
k562_episeq <- readRDS(file = "./temp/avo_k562_qc_filtered.rds")


## Correlation with Zeller et al. 2022 -----------------------------------------

# aggregate data for Zeller et al. into 1 tibble
zeller_chic_data <- Reduce(f = function(x, y) full_join(x, y, by = c("chr", "start", "end")), x = zeller_chic_data)
colnames(zeller_chic_data)[4:6] <- c("k36me3_zeller", "k9me3_zeller", "k27me3_zeller")

k562_episeq_binned_cuts <- k562_episeq %>% 
  dplyr::select(mod, chr, start, end, read_id) %>% 
  group_by(mod, chr) %>% 
  mutate(start = floor(start/50000)*50000, 
         end = start + 50000) %>% 
  group_by(mod, chr, start, end) %>% 
  summarise(counts = length(unique(read_id)))

# split by modification
k562_episeq_binned_cuts <- k562_episeq_binned_cuts %>% 
  ungroup() %>% 
  pivot_wider(names_from = mod, values_from = counts)

# combine datasets
comparison_k562_chic <- full_join(k562_episeq_binned_cuts, zeller_chic_data)

# calculate correlations
cor(comparison_k562_chic[, 4:9], use = "pairwise") %>% 
  pheatmap::pheatmap()

cor(comparison_k562_chic[, 4:9], use = "pairwise", method = "spearman")

# k27me3 0.9380113
# k36me3 0.81718392
# k9me3 0.94918115

chic_k27 <- comparison_k562_chic %>% 
  ggplot(aes(k27me3, k27me3_zeller)) +
  geom_bin2d(bins=100)+
  geom_smooth(method = 'lm', linetype = 2, col = "green1") +
  theme_bw() +
  scale_fill_viridis_c(option = "C", trans='log10') +
  scale_x_log10() +
  scale_y_log10() +
  ggtitle(label = "Zeller et al. H3K27me3") +
  xlab('ChiC counts (this study)') +
  ylab('ChiC counts (Zeller et al.)')+
  theme(legend.position = "none")

chic_k9 <- comparison_k562_chic %>% 
  ggplot(aes(k9me3, k9me3_zeller)) +
  geom_bin2d(bins=100)+
  geom_smooth(method = 'lm', linetype = 2, col = "green1") +
  theme_bw() +
  scale_fill_viridis_c(option = "C", trans='log10') +
  scale_x_log10() +
  scale_y_log10() +
  ggtitle(label = "Zeller et al. H3K9me3") +
  xlab('ChiC counts (this study)') +
  ylab('ChiC counts (Zeller et al.)')+
  theme(legend.position = "none")

chic_k36 <- comparison_k562_chic %>% 
  ggplot(aes(k36me3, k36me3_zeller)) +
  geom_bin2d(bins=100)+
  geom_smooth(method = 'lm', linetype = 2, col = "green1") +
  theme_bw() +
  scale_fill_viridis_c(option = "C", trans='log10') +
  scale_x_log10() +
  scale_y_log10() +
  ggtitle(label = "Zeller et al. H3K36me3") +
  xlab('ChiC counts (this study)') +
  ylab('ChiC counts (Zeller et al.)') +
  theme(legend.position = "none")


## Correlation with ENOCDE Chip-Seq --------------------------------------------

chip_data <- chip_data %>% 
  pivot_longer(cols = starts_with("h3"))

chip_data <- chip_data %>% 
  left_join(chip_controls)

chip_data <- chip_data %>% 
  mutate(norm = log2((value+1)/(control+1)))

chip_data_norm <- chip_data %>% 
  dplyr::select(chr, start, end, name, norm) %>% 
  pivot_wider(names_from = name, values_from = norm)

chip_data_counts <- chip_data %>% 
  dplyr::select(chr, start, end, name, value) %>% 
  pivot_wider(names_from = name, values_from = value)


test1 <- left_join(comparison_k562_chic, chip_data_norm)
test2 <- left_join(comparison_k562_chic, chip_data_counts)
test3 <- bind_cols(chip_data_counts[, 1:7], chip_data_norm[, 8:9])
test4 <- left_join(comparison_k562_chic, test3)

test4[, 4:15] %>% 
  as.matrix %>% 
  cor(use = "pairwise", method = "spearman")

# h3k36me3_rep2 0.6322502
# h3k9me3_rep1 0.7446769
# h3k27me3_rep1 0.9506541


chip_k36 <- test4 %>% 
  ggplot(aes(k36me3, h3k36me3_rep2)) +
  geom_bin2d(bins=100)+
  geom_smooth(method = 'lm', linetype = 2, col = "green1") +
  theme_bw() +
  scale_fill_viridis_c(option = "C", trans='log10') +
  scale_x_log10() +
  scale_y_log10() +
  ggtitle(label = "ENCODE H3K36me3") +
  xlab('ChiC counts (this study)') +
  ylab('ChiP signal (ENCODE)') +
  theme(legend.position = "none")


chip_k9 <- test4 %>% 
  filter(k9me3>1) %>% 
  ggplot(aes(k9me3, h3k9me3_rep1)) +
  geom_bin2d(bins=100)+
  geom_smooth(method = 'lm', linetype = 2, col = "green1") +
  theme_bw() +
  scale_fill_viridis_c(option = "C", trans='log10') +
  scale_x_log10() +
  #scale_y_log10() +
  ggtitle(label = "ENCODE H3K9me3") +
  xlab('ChiC counts (this study)') +
  ylab('Normalized ChiP signal (ENCODE)') +
  theme(legend.position = "none")

chip_k27 <- test4 %>% 
  ggplot(aes(k27me3, h3k27me3_rep1)) +
  geom_bin2d(bins=100)+
  geom_smooth(method = 'lm', linetype = 2, col = "green1") +
  theme_bw() +
  scale_fill_viridis_c(option = "C", trans='log10') +
  scale_x_log10() +
  scale_y_log10() +
  ggtitle(label = "ENCODE H3K27me3") +
  xlab('ChiC counts (this study)') +
  ylab('ChiP signal (ENCODE)') +
  theme(legend.position = "none")

(chic_k27 | chic_k9 | chic_k36) / (chip_k27 | chip_k9 | chip_k36)


# Assocation of 5mC / ChiC with repressed genes (R1.2) -------------------------

k562_episeq <- readRDS(file = "./temp/avo_k562_qc_filtered.rds")
k562_rnaseq <- readRDS(file = "./temp/encode_k562_bulk_rnaseq.rds")

# gene annotation
anno_rnaseq_promoters <- readRDS(file = "./temp/encode_annotation_promoters.rds")
anno_rnaseq_genebodies <- readRDS(file = "./temp/encode_annotation_gene_bodies.rds")
anno_rnaseq_biotype <- readRDS(file = "./temp/encode_annotation_biotype.rds") %>% 
  dplyr::select(entrez_id, symbol, biotype)

anno_rnaseq_genebodies <- anno_rnaseq_genebodies %>% 
  left_join(anno_rnaseq_biotype)



## Gene promoter methylation ---------------------------------------------------

# map RNA-seq data to promoters
k562_rnaseq_promoters <- left_join(
  k562_rnaseq, 
  anno_rnaseq_promoters) %>% 
  drop_na()

# find overlapping data points
k562_overlap_prom_episeq <- findOverlaps(
  query = makeGRangesFromDataFrame(k562_episeq[, c("chr", "start", "end")]),
  subject = makeGRangesFromDataFrame(k562_rnaseq_promoters[, c("seqnames", "start", "end")]))

# merge data sets
promoter_meth <- bind_cols(
  k562_rnaseq_promoters[k562_overlap_prom_episeq@to, ],
  k562_episeq[k562_overlap_prom_episeq@from, c("mod", "plate", "cell", "start", "meth")])

colnames(promoter_meth)[c(10, 17)] <- c("start", "cpg_pos")

# compress data set (calculate average)
promoter_meth_compressed <- promoter_meth %>% 
  drop_na() %>% 
  mutate(mod = paste0("H3", str_to_title(mod))) %>% 
  mutate(dist = abs(start - cpg_pos)-2000) %>%
  mutate(dist = ifelse(strand == "-", -dist, dist)) %>% 
  group_by(mod, exprs_quantile, dist) %>% 
  summarise(meth = mean(meth), 
            n = n())

# plot
promoter_meth_compressed %>%
  ggplot(aes(dist/1000, meth, col = mod)) +
  geom_line(alpha = 0.4) +
  geom_smooth() +
  facet_grid(rows = vars(exprs_quantile), 
             cols = vars(mod), scales = "free_y") +
  theme_bw(base_size = 14) +
  theme(legend.position = "none") +
  scale_colour_manual(values = color_scheme) +
  labs(x = "Distance to TSS (x 1,000 bp)", y = "Beta")


## Gene body methylation -------------------------------------------------------

# map RNA-seq data to gene bodies
k562_rnaseq_gb <- left_join(
  k562_rnaseq, 
  anno_rnaseq_genebodies) %>% 
  drop_na()

# find overlapping data points
k562_overlap_gb_episeq <- findOverlaps(
  query = makeGRangesFromDataFrame(k562_episeq[, c("chr", "start", "end")]),
  subject = makeGRangesFromDataFrame(k562_rnaseq_gb[, c("seqnames", "start", "end")]))

# merge data sets
gb_meth <- bind_cols(
  k562_rnaseq_gb[k562_overlap_gb_episeq@to, ],
  k562_episeq[k562_overlap_gb_episeq@from, c("mod", "plate", "cell", "start", "meth")])

colnames(gb_meth)[c(10, 18)] <- c("start", "cpg_pos")

# this is a large dataset - calculate avg. methylation as 
gb_meth_compressed <- gb_meth %>% 
  drop_na() %>% 
  mutate(mod = paste0("H3", str_to_title(mod))) %>% 
  mutate(dist = ifelse(strand == "-", abs(end - cpg_pos), abs(start - cpg_pos))) %>%
  mutate(dist_rel = dist/width)

gb_meth_compressed <- gb_meth_compressed %>% 
  mutate(dist_rel = round(dist_rel, digits = 3)) %>% 
  group_by(mod, exprs_quantile, dist_rel) %>% 
  summarise(meth = mean(meth), 
            n = n()) %>% 
  ungroup()

gb_meth_compressed %>% 
  ggplot(aes(dist_rel, meth, col = mod)) +
  geom_line(alpha = 0.4) +
  geom_smooth() +
  facet_grid(rows = vars(exprs_quantile), 
             cols = vars(mod)) +
  theme_bw(base_size = 14)  +
  theme(legend.position = "none") +
  scale_colour_manual(values = color_scheme) +
  labs(x = "Position in gene body (normalized)", y = "Beta")


# extract gene body methylation per gene
gb_meth_per_gene <- gb_meth %>% 
  group_by(ensembl_id, mod) %>% 
  summarise(symbol = symbol[1],
            length = length[1], 
            tpm = TPM[1],
            fpkm = FPKM[1],
            exprs_quantile = exprs_quantile[1], 
            width = width[1], 
            strand = strand[1], 
            meth = mean(meth), 
            n = n()) %>% 
  ungroup

# extract no. of detected CpGs per mark to normalize
total_reads_per_mark <- k562_episeq %>% 
  group_by(mod) %>% 
  summarise(n_total = n())

# normalize no. of detected CpGs by (1) total no. of CpGs per mark and (2) gene size
# also, z-score normalize TPM
gb_meth_per_gene <- gb_meth_per_gene %>% 
  left_join(total_reads_per_mark) %>% 
  mutate(n_norm = n/n_total/width) %>% 
  group_by(mod) %>% 
  mutate(n_norm_mean = mean(n_norm), 
         n_norm_sd = sd(n_norm)) %>% 
  ungroup %>% 
  mutate(n_norm_zscore = (n_norm - n_norm_mean)/n_norm_sd, 
         tpm_norm = log10(tpm))


gb_meth_per_gene %>% 
  mutate(znf = str_detect(pattern = "ZNF", symbol)) %>% 
  mutate(mod = paste0("H3", str_to_title(mod))) %>% 
  filter(n_norm_zscore < 5, 
         width > 1000) %>% 
  dplyr::select(ensembl_id, znf, mod, n_norm_zscore) %>% 
  pivot_wider(names_from = mod, values_from = n_norm_zscore) %>% 
  ggplot(aes(H3K9me3, H3K36me3)) +
  geom_point(alpha = 0.7) +
  geom_bin2d(bins = 200) +
  scale_fill_viridis_c(option = "C", trans='log10') +
  theme_bw(base_size = 14)

gb_meth_per_gene %>% 
  mutate(znf = ifelse(str_detect(pattern = "ZNF", symbol), "ZNF", "other")) %>% 
  mutate(mod = paste0("H3", str_to_title(mod))) %>% 
  filter(n_norm_zscore < 5, 
         width > 1000) %>% 
  dplyr::select(ensembl_id, znf, mod, n_norm_zscore) %>% 
  pivot_wider(names_from = mod, values_from = n_norm_zscore) %>% 
  ggplot(aes(H3K9me3, H3K36me3, col = znf)) +
  geom_point(alpha = 0.5) +
  theme_bw(base_size = 14) +
  scale_colour_manual(values = c("grey85", color_scheme_2[1]))

gb_meth_per_gene %>% 
  mutate(znf = ifelse(str_detect(pattern = "ZNF", symbol), "ZNF", "other")) %>% 
  mutate(mod = paste0("H3", str_to_title(mod))) %>% 
  filter(n_norm_zscore < 5, 
         width > 1000) %>% 
  dplyr::select(ensembl_id, znf, mod, n_norm_zscore) %>% 
  pivot_wider(names_from = mod, values_from = n_norm_zscore) %>% 
  mutate(dual = ifelse(H3K36me3 > 1 & H3K9me3 > 1, "yes", "no")) %>% 
  with(., table(znf, dual))
# pvalue 9.94135e-17

codecorated_znf_k36 <- gb_meth_per_gene %>% 
  filter(mod == "k36me3") %>% 
  filter(n_norm_zscore > 0.5) %>% 
  filter(str_detect(symbol, patter = "^ZNF")) %>% 
  pull(symbol)

codecorated_znf_k9 <- gb_meth_per_gene %>% 
  filter(mod == "k9me3") %>% 
  filter(n_norm_zscore > 0.5) %>% 
  filter(str_detect(symbol, patter = "^ZNF")) %>% 
  pull(symbol)

codecorated_znf <- intersect(codecorated_znf_k9, codecorated_znf_k36)



# gene body methylation for ZNF genes ------------------------------------------

gb_meth_znf <- gb_meth %>% 
  mutate(znf = ifelse(str_detect(symbol, pattern = "^ZNF"), "znf", "non_znf"), 
         codecorated = ifelse(symbol %in% codecorated_znf, "codecorated", "other"))

gb_meth_znf <- gb_meth_znf %>% 
  drop_na() %>% 
  mutate(mod = paste0("H3", str_to_title(mod))) %>% 
  mutate(dist = ifelse(strand == "-", abs(end - cpg_pos), abs(start - cpg_pos))) %>%
  mutate(dist_rel = dist/width) %>% 
  mutate(dist_rel = round(dist_rel, digits = 3)) %>% 
  group_by(mod, codecorated, exprs_quantile, dist_rel) %>% 
  summarise(meth = mean(meth), 
            n = n()) %>% 
  ungroup()


gb_meth_znf %>% 
  filter(mod %in% c("H3K9me3", "H3K36me3")) %>% 
  filter(exprs_quantile == 4) %>% 
  ggplot(aes(dist_rel, meth, col = mod, lty = codecorated)) +
  geom_line(alpha = 0.4) +
  geom_smooth() +
  facet_grid(cols = vars(mod)) +
  theme_bw(base_size = 14)  +
  scale_colour_manual(values = color_scheme[c(2,3)]) +
  labs(x = "Position in gene body (normalized)", y = "Beta")


# ChromHMM enrichments ---------------------------------------------------------



test <- read_tsv(file = "./input/encode/mnase-seq/overlap_mnase_chromHmm.bed", col_names = c("chr", "start", "end", "state_id", "count"))

test2 <- read_tsv(file = "./input/avo-k562/bam/k27me3/overlap_k27me3_chromHmm.bed", col_names = c("chr", "start", "end", "state_id", "count"))
test3 <- read_tsv(file = "./input/avo-k562/bam/k36me3/overlap_k36me3_chromHmm.bed", col_names = c("chr", "start", "end", "state_id", "count"))
test3b <- read_tsv(file = "./input/avo-k562/bam/k9me3/overlap_k9me3_chromHmm.bed", col_names = c("chr", "start", "end", "state_id", "count"))

test4 <- tibble(chr = test$chr, start = test$start, end = test$end, state_id = test$state_id, 
                background = test$count, k27me3 = test2$count, k36me3 = test3$count, k9me3 = test3b$count)

test5 <- test4 %>% 
  group_by(state_id) %>% 
  summarise(background = sum(background), 
            k27me3_norm = sum(k27me3)/sum(background), 
            k36me3_norm = sum(k36me3)/sum(background), 
            k9me3_norm = sum(k9me3)/sum(background))

test5 <- test5 %>% 
  mutate(k27me3_norm = log2(k27me3_norm/mean(k27me3_norm)), 
         k36me3_norm = log2(k36me3_norm/mean(k36me3_norm)), 
         k9me3_norm = log2(k9me3_norm/mean(k9me3_norm)))


test5 <- left_join(test5, chromhmm_meta)


reorder_state_id <- function(x){
  nr <- str_extract(x, pattern = "^[0-9]{1,3}") %>% 
    as.integer()
  x <- x[order(nr)]
  return(x)
}

test5 %>% 
  pivot_longer(cols = c(k27me3_norm, k36me3_norm, k9me3_norm)) %>% 
  mutate(state_id = fct_relevel(state_id, reorder_state_id)) %>% 
  mutate(name = str_extract(name, pattern = "k[0-9]{1,2}me3")) %>% 
  filter(!str_detect(state_id, "Gap")) %>% 
  ggplot(aes(state_id, value, fill = group)) +
  geom_col() +
  facet_grid(rows = vars(name)) +
  theme_classic(base_size = 14) +
  theme(legend.position = "top", legend.title = element_blank()) +
  theme(axis.text.x = element_text(vjust = 0.5, angle = 90, size = 8)) +
  labs(x = NULL, y = "Enrichment (log2)")

test5 %>% 
  mutate(state_id = fct_relevel(state_id, reorder_state_id)) %>% 
  pivot_longer(cols = c(k27me3_norm, k36me3_norm, k9me3_norm)) %>% 
  
  ggplot(aes(name, state_id, fill = value)) +
  geom_tile(colour="white", size=0.25) 


test5_mat <- test5 %>% 
  mutate(state_id = fct_relevel(state_id, reorder_state_id)) %>% 
  dplyr::select(k27me3_norm, k36me3_norm, k9me3_norm) %>% 
  as.matrix

rownames(test5_mat) <- test5$state_id
test5_mat <- test5_mat[reorder_state_id(rownames(test5_mat)), ]

pheatmap::pheatmap(test5_mat, cluster_rows = TRUE, border_color = "white", scale = "row", fontsize_row = 6)




# Seurat single-gene plots -----------------------------------------------------

library(Seurat)
library(Signac)

seurat <- readRDS(file = "./temp/avo_k562_seurat_all_datasets.rds")

roi = "chr19-19600000-19700000"
roi <- "chr22-22400000-22700000" # PRAME
roi <- "chr11-2100000-2400000" # IGF2 locus

cov_plot <- CoveragePlot(
  object = seurat,
  group.by = "dataset", 
  region = roi,
  annotation = FALSE,
  peaks = FALSE
) + scale_fill_manual(values = color_scheme)

gene_plot <- AnnotationPlot(
  object = seurat,
  region = roi
)

CombineTracks(
  plotlist = list(cov_plot,gene_plot),
  heights = c(10, 2),
  widths = c(10, 1)
)

ggsave(filename = paste0("./plots/gene_plots/avo_k562_", roi, ".pdf"), width = 10, height = 6)

