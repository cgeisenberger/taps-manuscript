library(data.table)
library(tidyverse)


setwd("~/Desktop/post-doc/Experiments/exp176-scEpi2-seq-mouse1-3-H3K27me3-prox-mid-distal/analysis/CpG/scTAPS/")
xa <- fread("~/Desktop/post-doc/Experiments/exp176-scEpi2-seq-mouse1-3-H3K27me3-prox-mid-distal/analysis/CpG/JvB-176-SI-mouse1-H3K27me3-scEpi2-seq-pl2/Z_methylation_combined.tsv.gz")

xb <- xa[, cell := sub(":.*", "", V1)][V5 %in% c("TA", "TT")][, 
                                 `:=`(
                                   methylated = as.integer(V4 == "Z"), 
                                   unmethylated = as.integer(V4 == "z"), 
                                   percentage = (V4 == "Z") / ((V4 == "Z") + (V4 == "z"))*100, 
                                   end = V3
                                 )][, .(cell, V2, V3, end, percentage, methylated, unmethylated)]


split_list <-  split(xb[, -1, with = FALSE], xb[[1]])


setwd("~/Desktop/post-doc/Experiments/exp176-scEpi2-seq-mouse1-3-H3K27me3-prox-mid-distal/analysis/CpG/scTAPS/")

for (name in names(split_list)) {
  fwrite(split_list[[name]], file = paste0(name, ".cov"), sep = "\t", col.names = FALSE)
}



# Define the file pattern (adjust the filename as needed)
filename <- "Z_methylation_combined.tsv.gz"  # Change this to the actual filename you're looking for

# Get all file paths matching the filename in different folders
file_paths <- list.files(path = "~/Desktop/post-doc/Experiments/exp176-scEpi2-seq-mouse1-3-H3K27me3-prox-mid-distal/analysis/CpG/", pattern = filename, full.names = TRUE, recursive = TRUE)

# Read all files and concatenate into one data.table
dt_combined <- rbindlist(lapply(file_paths, fread), use.names = TRUE, fill = TRUE)

# View the result
print(dt_combined)

gc()

xb <- dt_combined[, cell := sub(":.*", "", V1)][V5 %in% c("TA", "TT")][, 
                                                              `:=`(
                                                                methylated = as.integer(V4 == "Z"), 
                                                                unmethylated = as.integer(V4 == "z"), 
                                                                percentage = (V4 == "Z") / ((V4 == "Z") + (V4 == "z"))*100, 
                                                                end = V3
                                                              )][, .(cell, V2, V3, end, percentage, methylated, unmethylated)]

fwrite(xb, file = "allcells_methscan_input.tsv.gz", sep = "\t", compress = "gzip")

xb <- fread("~/Desktop/post-doc/Experiments/exp176-scEpi2-seq-mouse1-3-H3K27me3-prox-mid-distal/analysis/CpG/scTAPS/allcells_methscan_input.tsv.gz")
xb <- xb %>% dplyr::filter(V2 %in% c(1:19))
split_list <-  split(xb[, -1, with = FALSE], xb[[1]])

head(xb)
#JvB-176-SI-mouse1-H3K27me3-scEpi2-seq-pl2_130

  
setwd("~/Desktop/post-doc/Experiments/exp176-scEpi2-seq-mouse1-3-H3K27me3-prox-mid-distal/analysis/CpG/scTAPS/")

for (name in names(split_list)) {
  fwrite(split_list[[name]], file = paste0(name, ".cov"), sep = "\t", col.names = FALSE)
}



MD2 %>% 
  rownames_to_column() %>%
  separate(rowname, into = c("mouse_pl", "CN"), sep = "[_]") %>%
  separate(mouse_pl, into = c("mouse", "plate"), sep = "[-]") %>%
  mutate(cell = paste0("JvB-176-SI-", mouse, "-H3K27me3-scEpi2-seq-", plate, "_", CN )) %>%
  dplyr::select(cell, supercluster, subcluster)


DimPlot(chictaps.combined, reduction = "umap",group.by = "AP_axis",shuffle = T)+scale_color_manual(values= c("magenta", "blue", "cyan"))



