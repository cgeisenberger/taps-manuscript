library(tidyverse)
setwd(dir = "/mnt/ssd/private/cgeisenberger/projects/sc-epi-seq/")

source("./scripts/functions.R")

color_scheme <- c("#D674BA","#86BF88","#7A416A")
color_scheme_2 <- c("#0D3B66", "#FAF0CA", "#F95738", "#6C969D", "#2E5339")



# library stats (singleCellMultiOmics) -----------------------------------------

lib_stats <- readRDS(file = "./temp/avo_all_scmo_library_stats.rds")

# libraries which contain "qcChiC" were not TAPS converted
lib_stats <- lib_stats %>% 
  mutate(taps = ifelse(str_detect(library, "qcChIC"), "Unconverted", "TAPS"))

lib_stats %>% 
  dplyr::filter(taps == "TAPS") %>% 
  pull(mapping_rate) %>% 
  range

# save data as table for supplementary files
write_csv(
  x = lib_stats, 
  file = "./output/mapping_statistics_all_libraries.csv"
)

lib_stats %>% 
  dplyr::filter(taps == "TAPS") %>% 
  dplyr::select(library, read, ends_with("rate")) %>% 
  pivot_longer(cols = ends_with("rate")) %>% 
  mutate(name = str_replace(name, "_", " ")) %>% 
  mutate(name = str_to_title(name)) %>% 
  group_by(library, name) %>% 
  summarise(mean_r1_r2 = mean(value)) %>% 
  ggplot(aes(name, mean_r1_r2, fill = name)) +
  geom_violin(alpha = 0.5) +
  geom_boxplot(fill = "white", width = 0.1) +
  theme_bw(base_size = 12) +
  theme(legend.position = "none") +
  coord_flip() +
  labs(x = NULL, y = "Percent (%)")
ggsave(filename = "./plots/avo_all_library_stats.pdf", width = 5, height = 2.5)

lib_stats %>% 
  dplyr::select(library, read, taps, ends_with("rate")) %>% 
  pivot_longer(cols = ends_with("rate")) %>% 
  group_by(taps, name) %>% 
  summarise(mean_r1_r2 = mean(value))

lib_stats %>% 
  dplyr::filter(str_detect(library, "RPE")) %>% 
  ggplot(aes(taps, mapping_rate, fill = taps)) +
  geom_violin() +
  geom_boxplot(fill = "white", width = 0.1) +
  scale_fill_manual(values = color_scheme_2) +
  theme_bw(base_size = 16) +
  theme(legend.position = "none") +
  labs(x = "Treatment", y = "Mapping Rate (%)")
ggsave(filename = "./plots/avo_all_mapping_rate_unconverted_vs_taps.pdf", width = 4, height = 4)

lib_stats %>% 
  dplyr::filter(str_detect(library, "RPE")) %>% 
  ggplot(aes(taps, dedup_rate, fill = taps)) +
  geom_violin() +
  geom_boxplot(fill = "white", width = 0.1) +
  scale_fill_manual(values = color_scheme_2) +
  theme_bw(base_size = 16) +
  theme(legend.position = "none") +
  labs(x = "Treatment", y = "Deduplication Rate (%)")
ggsave(filename = "./plots/avo_all_deduplication_rate_unconverted_vs_taps.pdf", width = 4, height = 4)



# Lambda phage spike-in conversion rates ---------------------------------------

conversion_efficiency <- readRDS(file = "./temp/avo_all_conversion_efficiency.rds")

# summarise data for each library
conversion_efficiency <- conversion_efficiency %>% 
  group_by(library) %>% 
  summarise(m = sum(num_converted), 
            u = sum(num_unconverted), 
            beta = m/(m+u)*100)

# extract statistics
conversion_efficiency %>% 
  summarise(
    mean_efficiency = mean(beta), 
    min_efficiency = min(beta), 
    max_efficiency = max(beta)
  )

# plot
conversion_efficiency %>% 
  mutate(experiment = str_extract(library, "(RPE|K562|mouse)")) %>% 
  mutate(experiment = str_replace(experiment, "mouse", "Mouse \n Intestine")) %>% 
  ggplot(aes(experiment, beta)) +
  geom_boxplot(alpha = 0.1, outliers = FALSE, col = "lightblue") +
  geom_jitter(width = 0.3, col = "steelblue") +
  theme_bw(base_size = 12) +
  theme(legend.position = "none") +
  ylim(0, 100) +
  labs(x = NULL, y = "Spike-in conversion (%)")
ggsave(filename = "./plots/avo_all_conversion_efficiencies.pdf", width = 3, height = 5)

conversion_efficiency %>% 
  mutate(experiment = str_extract(library, "(RPE|K562|mouse)")) %>% 
  mutate(experiment = str_replace(experiment, "mouse", "Mouse \n Intestine")) %>% 
  ggplot(aes(experiment, beta)) +
  geom_boxplot(alpha = 0.1, outliers = FALSE, col = "lightblue") +
  geom_jitter(width = 0.3, col = "steelblue") +
  theme_bw(base_size = 16) +
  theme(legend.position = "none") +
  ylim(85, 100) +
  labs(x = NULL, y = NULL)
ggsave(filename = "./plots/avo_all_conversion_efficiencies_zooomed_in.pdf", width = 3, height = 5)



# Mismatch rate: samtools mpileup data for chr18 -------------------------------

pileup_data <- readRDS(file = "./temp/avo_all_mpileup_mismatch_rate.rds")

pileup_data <- pileup_data %>% 
  mutate(
    cpg = ifelse(str_detect(context, "CG"), "CpG", "CpH"),
    treatment = ifelse(str_detect(library, "qcChIC"), "Unconverted", "TAPS"),
    mod = str_extract(library, pattern = "k[0-9]{1,2}me3"))

pileup_data %>% 
  group_by(context, treatment) %>% 
  summarise(cpg = cpg[1], 
            taps = mean(taps)) %>% 
  ungroup %>% 
  mutate(taps = taps*100) %>% 
  pivot_wider(names_from = treatment, values_from = taps) %>% 
  ggplot(aes(Unconverted, TAPS, col = cpg)) +
  geom_point(size = 2) +
  geom_abline(slope = 1, intercept = 0, col = "grey", lty = 2) +
  scale_color_manual(values = color_scheme_2[c(3,4)]) +
  theme_bw(base_size = 16) +
  theme(legend.position = "inside", legend.position.inside = c(0.2, 0.8),
        legend.title = element_blank()) +
  labs(x = "Unconverted (%)", 
       y = "TAPS (%)")
ggsave(filename = "./plots/avo_all_mismatch_per_context_scatter.pdf", width = 4, height = 4)

pileup_data %>% 
  filter(treatment == "TAPS") %>% 
  group_by(context) %>% 
  summarise(cpg = cpg[1], 
            taps = mean(taps)) %>% 
  ungroup %>% 
  mutate(taps = taps*100) %>%  
  ggplot(aes(context, taps, fill = cpg)) +
  geom_col() +
  scale_fill_manual(values = color_scheme_2[c(3,4)]) +
  theme_bw(base_size = 14) +
  theme(legend.position = "none", axis.text.y = element_text(size = 6, vjust = 0.5)) +
  labs(x = "Context (3 bp)", y = "Mismatch Rate (%)") +
  coord_flip()
ggsave(filename = "./plots/avo_all_mismatch_per_context_bar.pdf", width = 4, height = 6)

pileup_data %>% 
  filter(cpg == "CpH") %>% 
  mutate(treatment = str_to_title(treatment)) %>% 
  ggplot(aes(treatment, taps * 100, fill = treatment)) +
  geom_violin(width = 0.5) +
  geom_boxplot(width = 0.1, fill = "white") +
  theme_bw(base_size = 16) +
  theme(legend.position = "none") +
  labs(x = "Treatment", y = "Mismatch rate (%)") +
  scale_fill_manual(values = color_scheme_2)
ggsave(filename = "./plots/avo_all_mismatch_rate_unconverted_vs_taps.pdf", width = 4, height = 4)


# clean up --------------------

rm(pileup_data, conversion_efficiency, insert_sizes, lib_stats)
