library(tidyverse)
library(data.table)
library(irlba)
library(patchwork)
library(ggdendro)


prcomp_iterative <- function(x, n=10, n_iter=50, min_gain=0.001, ...) {
  mse <- rep(NA, n_iter)
  na_loc <- is.na(x)
  x[na_loc] = 0  # zero is our first guess
  
  for (i in 1:n_iter) {
    prev_imp <- x[na_loc]  # what we imputed in the previous round
    # PCA on the imputed matrix
    pr <- prcomp_irlba(x, center = F, scale. = F, n = n, ...)
    # impute missing values with PCA
    new_imp <- (pr$x %*% t(pr$rotation))[na_loc]
    x[na_loc] <- new_imp
    # compare our new imputed values to the ones from the previous round
    mse[i] = mean((prev_imp - new_imp) ^ 2)
    # if the values didn't change a lot, terminate the iteration
    gain <- mse[i] / max(mse, na.rm = T)
    if (gain < min_gain) {
      message(paste(c("\n\nTerminated after ", i, " iterations.")))
      break
    }
  }
  pr$mse_iter <- mse[1:i]
  pr
}

fread("~/Desktop/post-doc/Experiments/exp176-scEpi2-seq-mouse1-3-H3K27me3-prox-mid-distal/analysis/CpG/scTAPS/compact_data/cell_stats.csv")%>% 
  mutate(QC = if_else((global_meth_frac >0.375 & global_meth_frac < 0.80 & n_obs >4000), true ="Pass", "Fail")) %>%
  ggplot(aes(x = global_meth_frac * 100, y = n_obs, color=QC)) +
  geom_point(alpha=0.5) +
  labs(x = "global DNA methylation %", y = "# of observed CpG sites")+
  scale_y_log10()+theme_bw()


meth_mtx <- read.csv("~/Desktop/post-doc/Experiments/exp176-scEpi2-seq-mouse1-3-H3K27me3-prox-mid-distal/analysis/CpG/scTAPS/VMR_matrix_int_sec/mean_shrunken_residuals.csv.gz", row.names=1) %>%
  as.matrix()

pca <- meth_mtx %>%
  scale(center = T, scale = F) %>%
  prcomp_iterative(n = 15)  # increase this value to e.g. 15 for real data sets

pca_tbl <- as_tibble(pca$x) %>% 
  add_column(cell=rownames(meth_mtx))

pca_tbl %>% 
  ggplot(aes(x = PC1, y = PC2)) +
  geom_point() +
  coord_fixed() +
  labs(title="PCA based on VMR methylation")


library(uwot)  # R package for UMAP

umap_obj <- uwot::umap(pca$x, min_dist=0.5, n_neighbors=5, seed=2, ret_nn=T)
umap_tbl <- umap_obj$embedding %>%
  magrittr::set_colnames(c("UMAP1", "UMAP2")) %>% 
  as_tibble() %>% 
  add_column(cell=rownames(meth_mtx))

umap_tbl %>% 
  ggplot(aes(x = UMAP1, y = UMAP2)) +
  geom_point() +
  coord_fixed() +
  labs(title="UMAP based on VMR methylation")+theme_bw()

library(igraph)  # R package for graph manipulation, also implements the Leiden algorithm

# get the edges of the neighbor graph from the UMAP object
neighbor_graph_edges <- 
  tibble(from = rep(1:nrow(umap_obj$nn$euclidean$idx), times=ncol(umap_obj$nn$euclidean$idx)),
         to = as.vector(umap_obj$nn$euclidean$idx),
         weight = as.vector(umap_obj$nn$euclidean$dist)) %>%
  filter(from != to) %>%
  mutate(from = rownames(meth_mtx)[from],
         to = rownames(meth_mtx)[to])

# run Leiden clustering
clust_obj <- neighbor_graph_edges %>%
  igraph::graph_from_data_frame(directed=F) %>% 
  igraph::cluster_leiden(resolution = 0.005)  # adjust the resolution parameter to your needs

# put the clustering results into a data frame (tibble) for plotting
clust_tbl <- tibble(
  leiden_cluster = as.character(clust_obj$membership),
  cell = clust_obj$names
) %>% 
  full_join(umap_tbl, by="cell")

clust_tbl %>% 
  ggplot(aes(x = UMAP1, y = UMAP2, color = leiden_cluster)) +
  geom_point() +
  coord_fixed()+theme_bw()+ggtitle("MethSCAN VMR UMAP")

clust_tbl %>%
  mutate(cell_group = case_when(
    leiden_cluster == "1" ~ "group1",
    leiden_cluster == "2" ~ "group2")) %>% 
  dplyr::select(cell, cell_group) %>%
  write_csv("~/Desktop/post-doc/Experiments/exp176-scEpi2-seq-mouse1-3-H3K27me3-prox-mid-distal/analysis/CpG/scTAPS/cell_groups.csv", col_names=F)

fread("~/Desktop/post-doc/Experiments/exp176-scEpi2-seq-mouse1-3-H3K27me3-prox-mid-distal/analysis/CpG/scTAPS/DMRs.bed") %>%
  group_by(V1) %>%summarize(n=n()) %>%
  ggplot()+
  geom_col(aes(x=V1,y=n))+theme_bw()




bb <- clust_tbl %>%
  mutate(cell_group = case_when(
    leiden_cluster == "1" ~ "group1",
    leiden_cluster == "2" ~ "group2")) %>% 
  dplyr::select(cell, cell_group) %>%
  separate(cell, into = c("name", "exp", "tissue", "mouse", "mark", "methods", "seq", "plate", "cell" ), sep= "[-]") %>%
  separate(plate, into = c("plate", "cell"), sep = "[_]") %>%
  group_by(mouse, plate, cell_group) %>%
  summarize(n=n()) %>%
  group_by(plate) %>%
  mutate(norm = n/ sum(n))%>%
  ggplot()+geom_col(aes(x=plate, y=norm, group=cell_group, fill=cell_group))+theme_bw()
  

cc <- clust_tbl %>%
  mutate(cell_group = case_when(
    leiden_cluster == "1" ~ "group1",
    leiden_cluster == "2" ~ "group2")) %>% 
  dplyr::select(cell, cell_group) %>%
  separate(cell, into = c("name", "exp", "tissue", "mouse", "mark", "methods", "seq", "plate", "cell" ), sep= "[-]") %>%
  separate(plate, into = c("plate", "cell"), sep = "[_]") %>%
  group_by(mouse, plate, cell_group) %>%
  summarize(n=n()) %>%
  group_by(mouse) %>%
  mutate(norm = n/ sum(n))%>%
  ggplot()+geom_col(aes(x=mouse, y=norm, fill=cell_group))+theme_bw()


  aa/bb/cc

chictaps.combined <- readRDS("~/Desktop/post-doc/Experiments/exp176-scEpi2-seq-mouse1-3-H3K27me3-prox-mid-distal/analysis/merge/cchictaps.chr1.19.integrated.rds")

MD2 <- chictaps.combined@meta.data %>% 
  rownames_to_column() %>%
  separate(rowname, into = c("mouse_pl", "CN"), sep = "[_]") %>%
  separate(mouse_pl, into = c("mouse", "plate"), sep = "[-]") %>%
  mutate(cell = paste0("JvB-176-SI-", mouse, "-H3K27me3-scEpi2-seq-", plate, "_", CN )) %>%
  dplyr::select(cell, supercluster, subcluster) %>%as.data.table()

clusterMS <- clust_tbl %>%
  mutate(cell_group = case_when(
    leiden_cluster == "1" ~ "group1",
    leiden_cluster == "2" ~ "group2")) %>% as.data.table()

merge.data.table(x = MD2, y=clusterMS, by.y = "cell", by.x = "cell") %>%
  ggplot(aes(x = UMAP1, y = UMAP2, color = supercluster)) +
  geom_point() +
  coord_fixed()+theme_bw()+ggtitle("MethSCAN VMR UMAP")+scale_color_manual(values = c("red2", "orange"))+


merge.data.table(x = MD2, y=clusterMS, by.y = "cell", by.x = "cell") %>%
  ggplot(aes(x = UMAP1, y = UMAP2, color = subcluster)) +
  geom_point() +
  coord_fixed()+theme_bw()+ggtitle("MethSCAN VMR UMAP")+scale_color_manual(values= c("red", "brown","orange1", "yellow2","yellow4",  "red4", "maroon" ))
  
  
fread("~/Desktop/post-doc/Experiments/exp176-scEpi2-seq-mouse1-3-H3K27me3-prox-mid-distal/analysis/CpG/scTAPS/DMRs.bed") %>%
  group_by(V1, V10) %>%summarize(n=n()) %>%
  mutate(CT = case_when(
    V10 == "group1" ~ "Intestinal+Secretory",
    V10 == "group2" ~ "Immune")) %>%
  ggplot()+
  geom_col(aes(x=V1 ,y=n, fill=CT))+theme_bw()+
  facet_grid(cols=vars(CT))+
  xlab("Chromosome")+
  ylab("Number of Differentially Methylated Regions")

head(xb)


fread("~/Desktop/post-doc/Experiments/exp176-scEpi2-seq-mouse1-3-H3K27me3-prox-mid-distal/analysis/CpG/scTAPS/DMRs.bed") %>%
  dplyr::filter()
  mutate(binsize =  V3-V2) %>%
  mutate(CT = case_when(
    V10 == "group1" ~ "Intestinal+Secretory",
    V10 == "group2" ~ "Immune")) %>%
  ggplot()+#geom_boxplot(aes(x=V10, y=binsize, fill=V10), alpha=0.2,outlier.shape = NA)+
  geom_boxplot(aes(x=CT, y=V4, fill=CT), outlier.shape = NA)+
  ggbeeswarm::geom_quasirandom(aes(x=CT, y=V4, col=CT), alpha =0.014)+
  theme_bw()+
  #facet_grid(cols=vars(CT))+
  ylab("Methylation difference in pecentage per DMR")
  
  
  
  fread("~/Desktop/post-doc/Experiments/exp176-scEpi2-seq-mouse1-3-H3K27me3-prox-mid-distal/analysis/CpG/scTAPS/DMRs.bed") %>%
dplyr::filter((V4 < - 5 | V4 > 5 ) & V6 > 100 & V7 > 100) %>%
    select(V1, V2, V3) %>%
    mutate(V1 = paste0("chr", V1))%>%
    setNames(c("chromosome", "chromStart", "chromEnd")) %>%
    write_tsv('~/Desktop/post-doc/Experiments/exp176-scEpi2-seq-mouse1-3-H3K27me3-prox-mid-distal/analysis/CpG/scTAPS/DMR_small.bed',col_names = F)
    

  clusterMS1 <- clust_tbl %>%
    mutate(cell_group = case_when(
      leiden_cluster == "1" ~ "Intestinal+Secretory",
      leiden_cluster == "2" ~ "Immune")) %>% as.data.table()
  
xb19 <- xb %>%
  dplyr::filter(V2 == "19" & V3 > 43553761 & V3 <  43656761) %>% as.data.table()

aa <- merge.data.table(x = xb19, y = clusterMS1, by.y = "cell", by.x = "cell") %>%
  #group_by(cell_group) %>%
  mutate(bin = round(V3/100)*100) %>%
  group_by(cell_group, bin) %>%
  summarize(meth_mean= mean(percentage)) %>%
  ggplot()+geom_col(aes(x= bin, y= meth_mean, fill= cell_group), width=100)+
  geom_vline(xintercept = 43593761, linetype=2, color= "blue")+
  geom_vline(xintercept = 43606761, linetype=2, color= "blue")+
  scale_fill_manual(values = c("red3", "yellow3"))+
  facet_grid(rows= vars(cell_group))+xlab('')+ylab('')+
  theme_bw()+theme(legend.position = "none")

aaa <- merge.data.table(x = xb19, y = clusterMS1, by.y = "cell", by.x = "cell") %>%
  group_by(cell_group) %>% mutate(bin = round(V3/5000)*5000) %>%
  group_by(cell_group, bin) %>%
  summarize(cov= n()) %>%
  ggplot()+geom_line(aes(x= bin, y= cov, color= cell_group), linewidth = 1)+
  geom_vline(xintercept = 43593761, linetype=2, color= "blue")+
  geom_vline(xintercept = 43606761, linetype=2, color= "blue")+
  scale_color_manual(values = c("red3", "yellow3"))+
  #facet_grid(rows= vars(cell_group))+
  xlab("")+ylab('CpG coverage (log10)')+scale_y_log10()+
  #ylab("mean CpG methylation 100bp-1 bin")+
  theme_bw()+theme(legend.position = "none")+ggtitle("Chr19:43553761-43656761")


  

xb7 <- xb %>%
  dplyr::filter(V2 == "7" & V3 > 70285227 & V3 <  70406227) %>% as.data.table()

bb <- merge.data.table(x = xb7, y = clusterMS1, by.y = "cell", by.x = "cell") %>%
  #group_by(cell_group) %>%
  mutate(bin = round(V3/100)*100) %>%
  group_by(cell_group, bin) %>%
  summarize(meth_mean= mean(percentage)) %>%
  ggplot()+geom_col(aes(x= bin, y= meth_mean, fill= cell_group), width=100)+
  geom_vline(xintercept = 70335227, linetype=2, color= "blue")+
  geom_vline(xintercept = 70356227, linetype=2, color= "blue")+
  scale_fill_manual(values = c("red3", "yellow3"))+
  facet_grid(rows= vars(cell_group))+ylab('Mean CpG methylation (100bp-1)')+xlab('')+
  theme_bw()+theme(legend.position = "none")

bbb <- merge.data.table(x = xb7, y = clusterMS1, by.y = "cell", by.x = "cell") %>%
  group_by(cell_group) %>% mutate(bin = round(V3/5000)*5000) %>%
  group_by(cell_group, bin) %>%
  summarize(cov= n()) %>%
  ggplot()+geom_line(aes(x= bin, y= cov, color= cell_group), linewidth = 1)+
  geom_vline(xintercept = 70335227, linetype=2, color= "blue")+
  geom_vline(xintercept = 70356227, linetype=2, color= "blue")+
  scale_color_manual(values = c("red3", "yellow3"))+
  #facet_grid(rows= vars(cell_group))+
  xlab("")+ylab('CpG Coverage (log10)')+scale_y_log10()+
  #ylab("mean CpG methylation 100bp-1 bin")+
  theme_bw()+theme(legend.position = "none")+ggtitle("Chr7:70285227-70406227")



xb3 <- xb %>%
  dplyr::filter(V2 == "3" & V3 > 129143175 & V3 <  129263175) %>% as.data.table() 


cc <- merge.data.table(x = xb3, y = clusterMS1, by.y = "cell", by.x = "cell") %>%
  #group_by(cell_group) %>%
  mutate(bin = round(V3/100)*100) %>%
  group_by(cell_group, bin) %>%
  summarize(meth_mean= mean(percentage)) %>%
  ggplot()+geom_col(aes(x= bin, y= meth_mean, fill= cell_group), width=100)+
  geom_vline(xintercept = 129193175, linetype=2, color= "blue")+
  geom_vline(xintercept = 129213175, linetype=2, color= "blue")+
  scale_fill_manual(values = c("red3", "yellow3"))+
  facet_grid(rows= vars(cell_group))+
  xlab("position (bp)")+ylab('')+
  #ylab("mean CpG methylation 100bp-1 bin")+
  theme_bw()+theme(legend.position = "none")

ccc <- merge.data.table(x = xb3, y = clusterMS1, by.y = "cell", by.x = "cell") %>%
  mutate(bin = round(V3/5000)*5000) %>%
  group_by(cell_group, bin) %>%
  summarize(cov= n()) %>%
  ggplot()+geom_line(aes(x= bin, y= cov, color= cell_group), linewidth = 1)+
  geom_vline(xintercept = 129193175, linetype=2, color= "blue")+
  geom_vline(xintercept = 129213175, linetype=2, color= "blue")+
  scale_color_manual(values = c("red3", "yellow3"))+
  #facet_grid(rows= vars(cell_group))+
  xlab("")+ylab('CpG Coverage(log10)')+scale_y_log10()+
  #ylab("mean CpG methylation 100bp-1 bin")+
  theme_bw()+theme(legend.position = "none")+ggtitle("Chr3:129193175-129213175")
  
  
acombi <- (aaa/aa)+plot_layout(heights = c(1,4))
bcombi <- (bbb/bb)+plot_layout(heights = c(1,4))
ccombi  <- (ccc/cc)+plot_layout(heights = c(1,4))

acombi/bcombi/ccombi
#chr6:52,118,157-52,326,456

(ccc/cc)+plot_layout(heights = c(2,6))




xb_hoxa <- xb %>%
  dplyr::filter(V2 == "6" & V3 > 52118157 & V3 <  52326456) %>% as.data.table() 


hoxa <- merge.data.table(x = xb_hoxa, y = clusterMS1, by.y = "cell", by.x = "cell") %>%
  #group_by(cell_group) %>%
  mutate(bin = round(V3/100)*100) %>%
  group_by(cell_group, bin) %>%
  summarize(meth_mean= mean(percentage)) %>%
  ggplot()+geom_col(aes(x= bin, y= meth_mean, fill= cell_group), width=100)+
  geom_vline(xintercept = 52231197, linetype=2, color= "blue")+ 
  geom_vline(xintercept = 52234939, linetype=2, color= "blue")+
  scale_fill_manual(values = c("red3", "yellow3"))+
  facet_grid(rows= vars(cell_group))+
  xlab("position (bp)")+ylab('')+
  #ylab("mean CpG methylation 100bp-1 bin")+
  theme_bw()+theme(legend.position = "none")+ggtitle("Hoxa cluster")

hoxaa <- merge.data.table(x = xb_hoxa, y = clusterMS1, by.y = "cell", by.x = "cell") %>%
  mutate(bin = round(V3/5000)*5000) %>%
  group_by(cell_group, bin) %>%
  summarize(cov= n()) %>%
  ggplot()+geom_line(aes(x= bin, y= cov, color= cell_group), linewidth = 1)+
  geom_vline(xintercept = 52231197, linetype=2, color= "blue")+ 
  geom_vline(xintercept = 52234939, linetype=2, color= "blue")+
  scale_color_manual(values = c("red3", "yellow3"))+
  #facet_grid(rows= vars(cell_group))+
  xlab("")+ylab('Coverage')+
  #ylab("mean CpG methylation 100bp-1 bin")+
  theme_bw()+theme(legend.position = "none")+ggtitle("Hoxa cluster")

(hoxaa/hoxa)+plot_layout(heights = c(1.5,6))

xb_hoxb <- xb %>%
  dplyr::filter(V2 == "11" & V3 > 96249484 & V3 <  96382560) %>% as.data.table() 


hoxb <- merge.data.table(x = xb_hoxb, y = clusterMS1, by.y = "cell", by.x = "cell") %>%
  #group_by(cell_group) %>%
  mutate(bin = round(V3/100)*100) %>%
  group_by(cell_group, bin) %>%
  summarize(meth_mean= mean(percentage)) %>%
  ggplot()+geom_col(aes(x= bin, y= meth_mean, fill= cell_group), width=100)+
  #geom_vline(xintercept = 52231197, linetype=2, color= "blue")+ 
  #geom_vline(xintercept = 52234939, linetype=2, color= "blue")+
  scale_fill_manual(values = c("red3", "yellow3"))+
  facet_grid(rows= vars(cell_group))+
  xlab("position (bp)")+ylab('')+
  #ylab("mean CpG methylation 100bp-1 bin")+
  theme_bw()+theme(legend.position = "none")+ggtitle("Hoxb cluster")

hoxbb <- merge.data.table(x = xb_hoxb, y = clusterMS1, by.y = "cell", by.x = "cell") %>%
  mutate(bin = round(V3/5000)*5000) %>%
  group_by(cell_group, bin) %>%
  summarize(cov= n()) %>%
  ggplot()+geom_line(aes(x= bin, y= cov, color= cell_group), linewidth = 1)+
  #geom_vline(xintercept = 52231197, linetype=2, color= "blue")+ 
  #geom_vline(xintercept = 52234939, linetype=2, color= "blue")+
  scale_color_manual(values = c("red3", "yellow3"))+
  #facet_grid(rows= vars(cell_group))+
  xlab("")+ylab('Coverage')+
  #ylab("mean CpG methylation 100bp-1 bin")+
  theme_bw()+theme(legend.position = "none")+ggtitle("Hoxb cluster")

(hoxbb/hoxb)+plot_layout(heights = c(1.5,6))


xb_hoxc <- xb %>%
  dplyr::filter(V2 == "15" & V3 > 102919917 & V3 <  103042185) %>% as.data.table() 


hoxc <- merge.data.table(x = xb_hoxc, y = clusterMS1, by.y = "cell", by.x = "cell") %>%
  #group_by(cell_group) %>%
  mutate(bin = round(V3/100)*100) %>%
  group_by(cell_group, bin) %>%
  summarize(meth_mean= mean(percentage)) %>%
  ggplot()+geom_col(aes(x= bin, y= meth_mean, fill= cell_group), width=100)+
  #geom_vline(xintercept = 52231197, linetype=2, color= "blue")+ 
  #geom_vline(xintercept = 52234939, linetype=2, color= "blue")+
  scale_fill_manual(values = c("red3", "yellow3"))+
  facet_grid(rows= vars(cell_group))+
  xlab("position (bp)")+ylab('')+
  #ylab("mean CpG methylation 100bp-1 bin")+
  theme_bw()+theme(legend.position = "none")+ggtitle("Hoxc cluster")

hoxcc <- merge.data.table(x = xb_hoxc, y = clusterMS1, by.y = "cell", by.x = "cell") %>%
  mutate(bin = round(V3/5000)*5000) %>%
  group_by(cell_group, bin) %>%
  summarize(cov= n()) %>%
  ggplot()+geom_line(aes(x= bin, y= cov, color= cell_group), linewidth = 1)+
  #geom_vline(xintercept = 52231197, linetype=2, color= "blue")+ 
  #geom_vline(xintercept = 52234939, linetype=2, color= "blue")+
  scale_color_manual(values = c("red3", "yellow3"))+
  #facet_grid(rows= vars(cell_group))+
  xlab("")+ylab('Coverage')+
  #ylab("mean CpG methylation 100bp-1 bin")+
  theme_bw()+theme(legend.position = "none")+ggtitle("Hoxc cluster")

(hoxcc/hoxc)+plot_layout(heights = c(1.5,6))


xb_hoxd <- xb %>%
  dplyr::filter(V2 == "2" & V3 > 74663995 & V3 <  74784006) %>% as.data.table() 


hoxd <- merge.data.table(x = xb_hoxd, y = clusterMS1, by.y = "cell", by.x = "cell") %>%
  #group_by(cell_group) %>%
  mutate(bin = round(V3/100)*100) %>%
  group_by(cell_group, bin) %>%
  summarize(meth_mean= mean(percentage)) %>%
  ggplot()+geom_col(aes(x= bin, y= meth_mean, fill= cell_group), width=100)+
  geom_vline(xintercept = 74721978, linetype=2, color= "blue")+ 
  geom_vline(xintercept = 74729160, linetype=2, color= "blue")+
  scale_fill_manual(values = c("red3", "yellow3"))+
  facet_grid(rows= vars(cell_group))+
  xlab("position (bp)")+ylab('')+
  #ylab("mean CpG methylation 100bp-1 bin")+
  theme_bw()+theme(legend.position = "none")+ggtitle("Hoxd cluster")

hoxdd <- merge.data.table(x = xb_hoxd, y = clusterMS1, by.y = "cell", by.x = "cell") %>%
  mutate(bin = round(V3/5000)*5000) %>%
  group_by(cell_group, bin) %>%
  summarize(cov= n()) %>%
  ggplot()+geom_line(aes(x= bin, y= cov, color= cell_group), linewidth = 1)+
  geom_vline(xintercept = 74721978, linetype=2, color= "blue")+ 
  geom_vline(xintercept = 74729160, linetype=2, color= "blue")+
  scale_color_manual(values = c("red3", "yellow3"))+
  #facet_grid(rows= vars(cell_group))+
  xlab("")+ylab('Coverage')+
  #ylab("mean CpG methylation 100bp-1 bin")+
  theme_bw()+theme(legend.position = "none")+ggtitle("Hoxd cluster")

(hoxdd/hoxd)+plot_layout(heights = c(1.5,6))

methQC <- fread("~/Desktop/post-doc/Experiments/exp176-scEpi2-seq-mouse1-3-H3K27me3-prox-mid-distal/analysis/CpG/scTAPS/compact_data/cell_stats.csv") %>%
dplyr::filter(global_meth_frac >0.375 & global_meth_frac < 0.80 & n_obs >4000) %>%
  as.data.table()

merge.data.table(x= MD2, y= methQC, by.y = "cell_name", by.x = "cell" ) %>%
  ggplot()+
  geom_boxplot(aes(x=supercluster, y =  global_meth_frac, fill=supercluster), outlier.shape = NA,alpha=0.3)+
  geom_jitter(aes(x=supercluster, y =  global_meth_frac, col=supercluster), alpha =0.3)+
  theme_bw()+
  scale_color_manual(values = c("red2","yellow2", "orange"))+
  scale_fill_manual(values = c("red2","yellow2", "orange"))+
  ylab("Mean CpG Methylation per cell")+
  xlab("")
  

merge.data.table(x= MD2, y= methQC, by.y = "cell_name", by.x = "cell" ) %>%
  ggplot()+
  geom_boxplot(aes(x=subcluster, y =  global_meth_frac, fill=subcluster), outlier.shape = NA,alpha=0.3)+
  geom_jitter(aes(x=subcluster, y =  global_meth_frac, col=subcluster), alpha =0.3)+
  theme_bw()+
  scale_color_manual(values= c("red", "brown","orange1", "yellow2","yellow4",  "red4", "maroon" ))+
  scale_fill_manual(values= c("red", "brown","orange1", "yellow2","yellow4",  "red4", "maroon" ))+
  ylab("Mean CpG Methylation per cell")+
  xlab("")+
  facet_grid(cols=vars(fct_rev(AP_axis)))

merge.data.table(x= MD2, y= methQC, by.y = "cell_name", by.x = "cell" ) %>%
  dplyr::filter(sex == "male") %>%
  ggplot()+
  geom_boxplot(aes(x=supercluster, y =  global_meth_frac, fill=supercluster), outlier.shape = NA,alpha=0.3)+
  geom_jitter(aes(x=supercluster, y =  global_meth_frac, col=supercluster), alpha =0.3)+
  theme_bw()+
  scale_color_manual(values = c("red2","yellow2", "orange"))+
  scale_fill_manual(values = c("red2","yellow2", "orange"))+
  ylab("Mean CpG Methylation per cell")+
  xlab("")+
  facet_grid(cols=vars(fct_rev(AP_axis)))


MD2 %>% dplyr::filter(supercluster %in% c("Immune")) %>%
  dplyr::select(cell) %>% 
  write_tsv("~/Desktop/post-doc/Experiments/exp176-scEpi2-seq-mouse1-3-H3K27me3-prox-mid-distal/analysis/CpG/scTAPS/keep_immune.txt", col_names = F)

header_cells <- fread("~/Desktop/post-doc/Experiments/exp176-scEpi2-seq-mouse1-3-H3K27me3-prox-mid-distal/analysis/CpG/scTAPS/compact_data/cell_stats.csv")%>% 
  dplyr::filter(global_meth_frac >0.375 & global_meth_frac < 0.80 & n_obs >4000) %>%
  select(cell_name) %>%
  pull()

MD2 %>%
  dplyr::filter(cell %in% header_cells) %>%
  mutate(cell_group = case_when(
    subcluster == "EECs" ~ "EECs",
    subcluster == "Goblet" ~ "Goblet")) %>%
  dplyr::select(cell, cell_group) %>%
  mutate(across(everything(), ~replace_na(.x, "-"))) %>% 
  write_csv("~/Desktop/post-doc/Experiments/exp176-scEpi2-seq-mouse1-3-H3K27me3-prox-mid-distal/analysis/CpG/scTAPS/EECs1_Goblet2_groups.csv", col_names=F)

MD2 %>%
  dplyr::filter(cell %in% header_cells) %>%
  mutate(cell_group = case_when(
    subcluster == "Intestinal1" ~ "Intestinal1",
    subcluster == "Goblet" ~ "Goblet")) %>%
  dplyr::select(cell, cell_group) %>%
  mutate(across(everything(), ~replace_na(.x, "-"))) %>% 
  write_csv("~/Desktop/post-doc/Experiments/exp176-scEpi2-seq-mouse1-3-H3K27me3-prox-mid-distal/analysis/CpG/scTAPS/Intestinal1_1_Goblet2_groups.csv", col_names=F)


MD2 %>%
  dplyr::filter(cell %in% header_cells) %>%
  mutate(cell_group = case_when(
    subcluster == "Intestinal2" ~ "Intestinal2",
    subcluster == "Goblet" ~ "Goblet")) %>%
  dplyr::select(cell, cell_group) %>%
  mutate(across(everything(), ~replace_na(.x, "-"))) %>% 
  write_csv("~/Desktop/post-doc/Experiments/exp176-scEpi2-seq-mouse1-3-H3K27me3-prox-mid-distal/analysis/CpG/scTAPS/Intestinal2_1_Goblet2_groups.csv", col_names=F)

MD2 %>%
  dplyr::filter(cell %in% header_cells) %>%
  mutate(cell_group = case_when(
    subcluster == "EECs" ~ "EECs",
    subcluster == "Intestinal2" ~ "Intestinal2")) %>%
  dplyr::select(cell, cell_group) %>%
  mutate(across(everything(), ~replace_na(.x, "-"))) %>% 
  write_csv("~/Desktop/post-doc/Experiments/exp176-scEpi2-seq-mouse1-3-H3K27me3-prox-mid-distal/analysis/CpG/scTAPS/EECs1_Intestinal2_2_groups.csv", col_names=F)

MD2 %>%
  dplyr::filter(cell %in% header_cells) %>%
  mutate(cell_group = case_when(
    subcluster == "EECs" ~ "EECs",
    subcluster == "Intestinal1" ~ "Intestinal1")) %>%
  dplyr::select(cell, cell_group) %>%
  mutate(across(everything(), ~replace_na(.x, "-"))) %>% 
  write_csv("~/Desktop/post-doc/Experiments/exp176-scEpi2-seq-mouse1-3-H3K27me3-prox-mid-distal/analysis/CpG/scTAPS/EECs1_Intestinal1_2_groups.csv", col_names=F)

MD2 %>%
  dplyr::filter(cell %in% header_cells) %>%
  mutate(cell_group = case_when(
    subcluster == "Intestinal1" ~ "Intestinal1",
    subcluster == "Intestinal2" ~ "Intestinal2")) %>%
  dplyr::select(cell, cell_group) %>%
  mutate(across(everything(), ~replace_na(.x, "-"))) %>% 
  write_csv("~/Desktop/post-doc/Experiments/exp176-scEpi2-seq-mouse1-3-H3K27me3-prox-mid-distal/analysis/CpG/scTAPS/Intestinal1_1_Intestinal2_2_groups.csv", col_names=F)




MD2 %>%
   dplyr::filter(cell %in% header_cells) %>%
  mutate(cell_group = case_when(
    subcluster == "B-cell" ~ "B-cell",
    subcluster == "T-cell" ~ "T-cell")) %>%
  dplyr::select(cell, cell_group) %>%
  mutate(across(everything(), ~replace_na(.x, "-"))) %>% 
  write_csv("~/Desktop/post-doc/Experiments/exp176-scEpi2-seq-mouse1-3-H3K27me3-prox-mid-distal/analysis/CpG/scTAPS/BCell1_Tcell2_groups.csv", col_names=F)



MD2 %>%
  dplyr::filter(cell %in% header_cells) %>%
  mutate(cell_group = case_when(
    subcluster == "B-cell" ~ "B-cell",
    subcluster == "Monocyte" ~ "Monocyte")) %>%
  dplyr::select(cell, cell_group) %>%
  mutate(across(everything(), ~replace_na(.x, "-"))) %>% 
  write_csv("~/Desktop/post-doc/Experiments/exp176-scEpi2-seq-mouse1-3-H3K27me3-prox-mid-distal/analysis/CpG/scTAPS/Bcell1_Monocyte2_groups.csv", col_names=F)


MD2 %>%
  dplyr::filter(cell %in% header_cells) %>%
  mutate(cell_group = case_when(
    subcluster == "T-cell" ~ "T-cell",
    subcluster == "Monocyte" ~ "Monocyte")) %>%
  dplyr::select(cell, cell_group) %>%
  mutate(across(everything(), ~replace_na(.x, "-"))) %>% 
  write_csv("~/Desktop/post-doc/Experiments/exp176-scEpi2-seq-mouse1-3-H3K27me3-prox-mid-distal/analysis/CpG/scTAPS/Tcell1_Monocyte2_groups.csv", col_names=F)

fread("~/Desktop/post-doc/Experiments/exp176-scEpi2-seq-mouse1-3-H3K27me3-prox-mid-distal/analysis/CpG/scTAPS/DMRs.bed") %>%
  mutate(V10 = if_else(V10 == "group1", true = "Intestinal", false = "Immune")) %>%
  write_tsv("~/Desktop/post-doc/Experiments/exp176-scEpi2-seq-mouse1-3-H3K27me3-prox-mid-distal/analysis/CpG/scTAPS/DMRs_supercluster.bed", col_names = F)


list.files("~/Desktop/post-doc/Experiments/exp176-scEpi2-seq-mouse1-3-H3K27me3-prox-mid-distal/analysis/CpG/scTAPS/", pattern = "DMRs.bed", full.names = TRUE) %>%
  map_dfr(~ fread(.x)[, filename := basename(.x)])

allDMRs <- list.files("~/Desktop/post-doc/Experiments/exp176-scEpi2-seq-mouse1-3-H3K27me3-prox-mid-distal/analysis/CpG/scTAPS/", pattern = "DMRs.bed", full.names = TRUE) %>%
  map(~ fread(.x, colClasses = list(character = 1))[, filename := basename(.x)]) %>%
  rbindlist(fill = TRUE) %>% 
  #mutate(filename = gsub(filename, "_DMRs.bed", "")) %>%
  as.data.table()

allDMRs  %>%
  dplyr::filter(V5>25 & V6>25) %>%
  #dplyr::filter()
  mutate(DMRs = ( V5 > 25 & V6 > 25) & V11 < 0.001,
         percentage = V8-V9)%>% #& (V5+V6 > 40) ) %>%
  ggplot()+
  ggrastr::rasterise( 
    geom_jitter(aes(y= -log10(V11), x = (percentage*100), col = DMRs), alpha=0.5, size=0.5), dpi = 900)+
  #facet_wrap(~(filename), ncol = 2)+
  scale_color_manual(values = c("grey60", "red"))+
  theme_bw()+xlab("Difference in Methylation within DMR (%)")+
  ylab("-log10(adj.p-value)")+
  ylim(0,15.5)

xc <- allDMRs %>% dplyr::filter(filename != "Supercluster_DMRs.bed") %>%
  dplyr::filter(V11 < 0.05 & V5 > 200 & V6 > 200) %>%
  arrange(V4) 

allDMRs  %>%
  #dplyr::filter(filename == "Supercluster_DMRs.bed") %>%
  #dplyr::filter()
  mutate(DMRs = (V4 > -5) & V12 < 0.05)  %>%
  dplyr::filter(DMRs) %>%
  summarise(n=n())

  

ggplot()+geom_point(aes(x=(V8-V9), y= -log10(V11)))+
  facet_grid(cols=vars(filename))


xb_tet <- xb %>%
  dplyr::filter(V2 == "10" & V3 > 62815000 & V3 <  62840530) %>% as.data.table() 


TETa <- merge.data.table(x = xb_tet, y = clusterTET, by.y = "cell", by.x = "cell") %>%
  #group_by(cell_group) %>%
  mutate(bin = round(V3/100)*100) %>%
  group_by(cell_group, bin) %>%
  summarize(meth_mean= mean(percentage)) %>%
  ggplot()+geom_col(aes(x= bin, y= meth_mean, fill= cell_group), width=100)+
  geom_vline(xintercept = 62817530, linetype=2, color= "blue")+ 
  geom_vline(xintercept = 62824530, linetype=2, color= "blue")+
  geom_vline(xintercept = 62831530, linetype=2, color= "blue4")+ 
  geom_vline(xintercept = 62834530, linetype=2, color= "blue4")+
  scale_fill_manual(values = c("red3", "yellow3"))+
  facet_grid(rows= vars(cell_group))+
  xlab("position (bp)")+ylab('')+
  #ylab("mean CpG methylation 100bp-1 bin")+
  theme_bw()+theme(legend.position = "none")+ggtitle("Tet1")

TETaa <- merge.data.table(x= xb_tet, y = clusterTET, by.y = "cell", by.x = "cell") %>%
  mutate(bin = round(V3/5000)*5000) %>%
  group_by(cell_group, bin) %>%
  summarize(cov= n()) %>%
  ggplot()+geom_line(aes(x= bin, y= cov, color= cell_group), linewidth = 1)+
  geom_vline(xintercept = 62817530, linetype=2, color= "blue")+ 
  geom_vline(xintercept = 62824530, linetype=2, color= "blue")+
  geom_vline(xintercept = 62831530, linetype=2, color= "blue4")+ 
  geom_vline(xintercept = 62834530, linetype=2, color= "blue4")+
  scale_color_manual(values = c("red3", "yellow3"))+
  #facet_grid(rows= vars(cell_group))+
  xlab("")+ylab('Coverage')+
  #ylab("mean CpG methylation 100bp-1 bin")+
  theme_bw()+theme(legend.position = "none")+ggtitle("Tet1")

(TETaa/TETa)+plot_layout(heights = c(1.5,6))


xb_2G <- xb %>%
  dplyr::filter(V2 == "7" & V3 > 141887227 & V3 <  141932227) %>% as.data.table() 



Gob2a <- merge.data.table(x = xb_2G, y = clusterTET, by.y = "cell", by.x = "cell") %>%
#group_by(cell_group) %>%
  mutate(bin = round(V3/100)*100) %>%
  group_by(cell_group, bin) %>%
  summarize(meth_mean= mean(percentage)) %>%
  ggplot()+geom_col(aes(x= bin, y= meth_mean, fill= cell_group), width=100)+
  geom_vline(xintercept = 141907227, linetype=2, color= "blue")+ 
  geom_vline(xintercept = 141912227, linetype=2, color= "blue")+
  scale_fill_manual(values = c("red3", "yellow3"))+
  facet_grid(rows= vars(cell_group))+
  xlab("position (bp)")+ylab('')+
  #ylab("mean CpG methylation 100bp-1 bin")+
  theme_bw()+theme(legend.position = "none")+ggtitle("Tet1")

Gob2aa <- merge.data.table(x= xb_2G, y = clusterTET, by.y = "cell", by.x = "cell") %>%
  mutate(bin = round(V3/2500)*2500) %>%
  group_by(cell_group, bin) %>%
  summarize(cov= n()) %>%
  ggplot()+geom_line(aes(x= bin, y= cov, color= cell_group), linewidth = 1)+
  geom_vline(xintercept = 141907227, linetype=2, color= "blue")+ 
  geom_vline(xintercept = 141912227, linetype=2, color= "blue")+
  scale_color_manual(values = c("red3", "yellow3"))+
  #facet_grid(rows= vars(cell_group))+
  xlab("")+ylab('Coverage')+
  #ylab("mean CpG methylation 100bp-1 bin")+
  theme_bw()+theme(legend.position = "none")+ggtitle("Tet1")

(Gob2aa/Gob2a)+plot_layout(heights = c(1.5,6))

#2  98663518  98668518

clusterHoxd3 <- MD2 %>%
  dplyr::filter(cell %in% header_cells) %>%
  mutate(cell_group = case_when(
    subcluster == "Intestinal1" ~ "Intestinal1",
    subcluster == "Intestinal2" ~ "Intestinal2")) %>%
  dplyr::select(cell, cell_group) %>% drop_na()%>% as.data.table()


xb_hoxd3 <- xb %>%
  dplyr::filter(V2 == "2" & V3 > 74732518 & V3 <  74765518) %>% as.data.table() 


#2  74742518  74755518 -7.060775   562   378   354 0.25478935 0.42227040 Intestinal1 4.103384e-12


hoxd3a <- merge.data.table(x = xb_hoxd3, y = clusterHoxd3, by.y = "cell", by.x = "cell") %>%
  #group_by(cell_group) %>%
  mutate(bin = round(V3/100)*100) %>%
  group_by(cell_group, bin) %>%
  summarize(meth_mean= mean(percentage)) %>%
  ggplot()+geom_col(aes(x= bin, y= meth_mean, fill= cell_group), width=100)+
  geom_vline(xintercept = 74742518, linetype=2, color= "blue")+ 
  geom_vline(xintercept = 74755518, linetype=2, color= "blue")+
  scale_fill_manual(values = c("yellow3", "yellow4"))+
  facet_grid(rows= vars(cell_group))+
  xlab("position (bp)")+ylab('')+
  ylab("mean CpG methylation 100bp-1 bin")+
  theme_bw()+theme(legend.position = "none")

hoxd3aa <-  merge.data.table(x = xb_hoxd3, y = clusterHoxd3, by.y = "cell", by.x = "cell") %>%
  mutate(bin = round(V3/500)*500) %>%
  group_by(cell_group, bin) %>%
  summarize(cov= n()) %>%
  ungroup()%>% 
  #group_by(cell_group) %>%
  #mutate(cov = scale(cov)) %>%
  ggplot()+geom_line(aes(x= bin, y= cov, col= cell_group))+
  geom_vline(xintercept = 74742518, linetype=2, color= "blue")+ 
  geom_vline(xintercept = 74755518, linetype=2, color= "blue")+
  scale_color_manual(values = c("yellow3", "yellow4"))+
  facet_grid(rows= vars(cell_group))+
  xlab("")+ylab('Coverage')+
  #ylab("mean CpG methylation 100bp-1 bin")+
  theme_bw()+theme(legend.position = "none")+ggtitle("Hoxd3")

(hoxd3aa/hoxd3a)+plot_layout(heights = c(2,6))



16
95451798
95455798
-
  
  
  clusterErg <- MD2 %>%
  dplyr::filter(cell %in% header_cells) %>%
  mutate(cell_group = case_when(
    subcluster == "T-cell" ~ "T-cell",
    subcluster == "B-cell" ~ "B-cell")) %>%
  dplyr::select(cell, cell_group) %>% drop_na()%>% as.data.table()


xb_Erg <- xb %>%
  dplyr::filter(V2 == "1" & V3 > 63260828-5e3 & V3 <  63268828+5e3) %>% as.data.table() 


#chr1:63,260,828-63,268,828
#CTCF 63264982-63265218
Erg3a <- merge.data.table(x = xb_Erg, y = clusterErg, by.y = "cell", by.x = "cell") %>%
  #group_by(cell_group) %>%
  mutate(bin = round(V3/100)*100) %>%
  group_by(cell_group, bin) %>%
  summarize(meth_mean= mean(percentage)) %>%
  ggplot()+geom_col(aes(x= bin, y= meth_mean, fill= cell_group), width=100)+
  geom_vline(xintercept = 63260828, linetype=2, color= "blue")+ 
  geom_vline(xintercept = 63268828, linetype=2, color= "blue")+
  geom_vline(xintercept = 63264982, linetype=2, color= "red")+ 
  geom_vline(xintercept = 63265218, linetype=2, color= "red")+
  scale_fill_manual(values = c("yellow3", "yellow4"))+
  facet_grid(rows= vars(cell_group))+
  xlab("position (bp)")+ylab('')+
  ylab("mean CpG methylation 100bp-1 bin")+
  theme_bw()+theme(legend.position = "none")

Erg3aa <-  merge.data.table(x = xb_Erg, y = clusterErg, by.y = "cell", by.x = "cell") %>%
  mutate(bin = round(V3/500)*500) %>%
  group_by(cell_group, bin) %>%
  summarize(cov= n()) %>%
  ungroup()%>% 
  group_by(cell_group) %>%
  mutate(cov = scale(cov)) %>%
  ggplot()+geom_line(aes(x= bin, y= cov, col= cell_group))+
  geom_vline(xintercept = 63260828, linetype=2, color= "blue")+ 
  geom_vline(xintercept = 63268828, linetype=2, color= "blue")+
  geom_vline(xintercept = 63264982, linetype=2, color= "red")+ 
  geom_vline(xintercept = 63265218, linetype=2, color= "red")+
  scale_color_manual(values = c("yellow3", "yellow4"))+
  facet_grid(rows= vars(cell_group))+
  xlab("")+ylab('Coverage')+
  #ylab("mean CpG methylation 100bp-1 bin")+
  theme_bw()+theme(legend.position = "none")+ggtitle("CTCF near Zdbf2")

(Erg3aa/Erg3a)+plot_layout(heights = c(2,6))

  
xe <- allDMRs %>% #dplyr::filter(filename != "Supercluster_DMRs.bed")%>%
  dplyr::filter(abs(V8-V9) > 0.25 & V11 < 0.001 & V5> 30 & V6 > 30) %>%
  arrange(-abs(V4)) %>%
  group_by(filename) %>%
  slice_head(n = 25) %>% 
  setNames(c("chr", "start",  "end", "score", "CpGcov1", "CpGcov2", "NoCpGs", "fCpGme1", "fCpGme2", "Subcluster", "pval", "qval", "original")) %>%
  as.data.table()

xb <- fread("~/Desktop/post-doc/Experiments/exp176-scEpi2-seq-mouse1-3-H3K27me3-prox-mid-distal/analysis/CpG/scTAPS/allcells_methscan_input.tsv.gz")
xb <- xb %>% rename(chr = V2, start = V3) %>% as.data.table()

setkeyv(x = xb, c('chr', 'start', 'end'))

setkeyv(x = xe,c('chr', 'start', 'end'))

xf <- foverlaps(xb, xe) %>% drop_na() 

xg <- xf %>% group_by(cell, score) %>%
  summarize(fme = mean(percentage)) %>% as.data.table()

MD3 <- chictaps.integrated@meta.data %>%
  rownames_to_column() %>%
  separate(rowname, into = c("mouse_pl", "CN"), sep = "[_]") %>%
  separate(mouse_pl, into = c("mouse", "plate"), sep = "[-]") %>%
  mutate(cell = paste0("JvB-176-SI-", mouse, "-H3K27me3-scEpi2-seq-", plate, "_", CN )) %>%
  dplyr::select(cell, supercluster, subcluster) %>% as.data.table()

merge.data.table(x=MD3, y= xg, by.y = "cell", by.x = "cell") %>%
  group_by(subcluster, score) %>%
  summarize(fme = mean(fme)) %>%
  ggplot()+
  geom_tile(aes(x=as.factor(score), y=subcluster, fill=fme))+
  scale_fill_viridis_c(option="B")+
  theme_bw()
  
xe2 <- allDMRs %>% #dplyr::filter(filename != "Supercluster_DMRs.bed")%>%
  dplyr::filter(abs(V8-V9) > 0.25 & V11 < 0.001 & V5> 30 & V6 > 30) %>%
  arrange(-abs(V4)) %>%
  group_by(filename) %>%
  slice_head(n = 25) %>% 
  setNames(c("chr", "start",  "end", "score", "CpGcov1", "CpGcov2", "NoCpGs", "fCpGme1", "fCpGme2", "Subcluster", "pval", "qval", "original")) %>%
  as.data.table()

setkeyv(x = xe2,c('chr', 'start', 'end'))

xf2 <- foverlaps(xb, xe2) %>% drop_na() 

xg2 <- xf2 %>% group_by(cell, score) %>%
  summarize(fme = mean(percentage)) %>% as.data.table()

MD3 <- chictaps.integrated@meta.data %>%
  rownames_to_column() %>%
  separate(rowname, into = c("mouse_pl", "CN"), sep = "[_]") %>%
  separate(mouse_pl, into = c("mouse", "plate"), sep = "[-]") %>%
  mutate(cell = paste0("JvB-176-SI-", mouse, "-H3K27me3-scEpi2-seq-", plate, "_", CN )) %>%
  dplyr::select(cell, supercluster, subcluster) %>% as.data.table()

merge.data.table(x=MD3, y= xg2, by.y = "cell", by.x = "cell") %>%
  group_by(subcluster, score) %>%
  summarize(fme = mean(fme)) %>%
  select(subcluster, score, fme)
  ggplot()+
  geom_tile(aes(x=as.factor(score), y=subcluster, fill=fme))+
  scale_fill_viridis_c(option="B")+
  theme_bw()


dist_matrix <- 
  merge.data.table(x=MD3, y= xg2, by.y = "cell", by.x = "cell") %>%
  group_by(subcluster, score) %>%
  summarize(fme = mean(fme)) %>%
  select(subcluster, score, fme) %>% 
  pivot_wider(names_from = subcluster, values_from = fme, values_fill = 0) %>%
  column_to_rownames("score") %>%
  as.matrix()

# Perform hierarchical clustering on rows and columns
row_hc <- hclust(dist(dist_matrix), method = "ward.D2")
col_hc <- hclust(dist(t(dist_matrix)), method = "ward.D2")

# Get hierarchical order
row_order <- row_hc$order
col_order <- col_hc$order

# Reorder row and column names based on clustering
ordered_row_ids <- rownames(dist_matrix)[row_order]
ordered_col_ids <- colnames(dist_matrix)[col_order]

heatmap_data <- as.data.frame(as.table(dist_matrix)) %>%
  rename(score = Var1, subcluster = Var2, fme = Freq) %>%
  mutate(
    score = factor(score, levels = ordered_row_ids),
    subcluster = factor(subcluster, levels = ordered_col_ids)
  )


heatmap <- ggplot(heatmap_data, aes(x = score, y = subcluster, fill = fme)) +
  geom_tile() +
  scale_fill_viridis_c(option="B", name = "Fraction \nCpG methylation") +  # Perceptual color scale
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  labs(
       x = "DMRs", y = "Subcluster")+
  theme(
    axis.text.x = element_blank(),  # Remove x-axis labels
    axis.ticks.x = element_blank()  # Remove x-axis ticks
  )

row_dendro <- ggdendro::dendro_data(row_hc)
col_dendro <- ggdendro::dendro_data(col_hc)

# Row dendrogram (left side)
row_dendrogram_plot <- ggplot(ggdendro::segment(row_dendro)) +
  geom_segment(aes(x = y, y = x, xend = yend, yend = xend)) +
  theme_void() +
  coord_flip()  # Flip to match heatmap orientation

# Column dendrogram (top)
col_dendrogram_plot <- ggplot(ggdendro::segment(col_dendro)) +
  geom_segment(aes(x = x, y = y, xend = xend, yend = yend)) +
  theme_void()
  
  
(row_dendrogram_plot / heatmap)+plot_layout(heights = c(1,4))

(col_dendrogram_plot|row_dendrogram_plot)/(col_dendrogram_plot | heatmap)+plot_layout(widths  = c(1,6), heights = c(1,4))



cd <- (col_dendrogram_plot|row_dendrogram_plot)+plot_layout(widths  = c(1,6))
rd <- (col_dendrogram_plot | heatmap)+plot_layout(widths  = c(1,6))

(cd/rd)+plot_layout( heights = c(1,4))
