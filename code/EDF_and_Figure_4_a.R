rm(list=ls())
library(Signac)
library(Seurat)
library(GenomicRanges)
library(patchwork)
library(clusterProfiler)
library(org.Mm.eg.db)
library(enrichplot)
library(tidyverse)
library(data.table)
library(AnnotationHub)
library(AnnotationHub)

#install.packages("Matrix")
#if (!require("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")

#BiocManager::install("GenomicRanges")
#BiocManager::install("BiocGenerics")
#BiocManager::install("clusterProfiler")
#BiocManager::install("fgsea")
#BiocManager::install("org.Mm.eg.db")
#BiocManager::install("enrichplot")
#BiocManager::install("Rsamtools")
#BiocManager::install("ensembldb") 
#BiocManager::install("AnnotationHub")
#BiocManager::install("biovizBase")

# load count matrix 50kb bins, cells with at least 10^2.5 CpGs
q1<-readRDS("~/Desktop/post-doc/Experiments/exp176-scEpi2-seq-mouse1-3-H3K27me3-prox-mid-distal/analysis/mm10_50kb_filtered.cov")
colnames(q1)<-gsub("sc50kb_mm10_support_JvB-176-SI-","",colnames(q1))

# prepare metadata
MD<-data.frame(matrix(0,ncol=3,nrow=ncol(q1)))
rownames(MD)<-colnames(q1)
colnames(MD)<-c("mouse","plate","AP_axis")

ss<-strsplit(rownames(MD),"-")
vec<-NULL
for (i in 1:length(ss)){
  MD[i,1]<-ss[[i]][1]
  vec[i]<-ss[[i]][2]
}
ss<-strsplit(rownames(MD),"_")
for (i in 1:length(ss)){
  MD[i,2]<-ss[[i]][1]
}

qqq<-read.csv('~/Desktop/post-doc/Experiments/exp176-scEpi2-seq-mouse1-3-H3K27me3-prox-mid-distal/analysis/chic_bc_well.tsv',sep='\t')
QQQ<-read.csv('~/Desktop/post-doc/Experiments/exp176-scEpi2-seq-mouse1-3-H3K27me3-prox-mid-distal/analysis/JvB-176-idx.tsv',sep='\t')
rownames(qqq)<-qqq[,3]

QQQ[,12]<-qqq[QQQ[,1],1]
QQQ[,13]<-paste(QQQ[,4],QQQ[,12],sep="_")
QQQ[,13]<-gsub("JvB-176-SI-","",QQQ[,13])
QQQ[,13]<-gsub("H3K27me3-scEpi2-seq-","",QQQ[,13])
rownames(QQQ)<-QQQ[,13]
QQ<-QQQ[rownames(MD),]
MD[,3]<-QQ$APaxis
MD<-MD[!is.na(MD[,3]),]
MD<-MD[MD[,1]=="mouse3",]
q1<-q1[,rownames(MD)]
rownames(q1)<-gsub("chr","",rownames(q1))


#fread("~/Desktop/post-doc/Experiments/exp176-scEpi2-seq-mouse1-3-H3K27me3-prox-mid-distal/analysis/merge/merge.bed.gz") %>%
#  distinct(V4)

chrom_assay <- Signac::CreateChromatinAssay(
  counts = q1,
  sep = c(":", "-"),
  fragments = "~/Desktop/post-doc/Experiments/exp176-scEpi2-seq-mouse1-3-H3K27me3-prox-mid-distal/analysis/mouse3.fragments.sorted.bed.gz",
  min.cells = 10,
  min.features = 200
)

chictaps <- CreateSeuratObject(
  counts = chrom_assay,
  assay = "peaks",
  meta.data = MD
)

peaks.keep <- seqnames(granges(chictaps)) %in% as.character(c(1:19))
chictaps <- chictaps[as.vector(peaks.keep), ]



ah <- AnnotationHub()
query(ah, "EnsDb.mmusculus.v99")
ensdb_v99 <- ah[["AH78811"]]

# extract gene annotations from EnsDb
annotations <- GetGRangesFromEnsDb(ensdb = ensdb_v99)
genome(annotations) <- "mm10"
Annotation(chictaps) <- annotations


# start Signac analysis
chictaps <- RunTFIDF(chictaps)
chictaps <- FindTopFeatures(chictaps, min.cutoff = 'q0')
chictaps <- RunSVD(chictaps)

# check for components that are strongly correlating with sequence depth
# for mouse3 I remove first component

DepthCor(chictaps)

chictaps <- RunUMAP(object = chictaps, reduction = 'lsi', dims = c(2:30))
chictaps <- FindNeighbors(object = chictaps, reduction = 'lsi', dims = c(2:30))
chictaps <- FindClusters(object = chictaps, verbose = FALSE, algorithm = 3, resolution = .1)
chictaps3 <- chictaps
DimPlot(chictaps, reduction = "umap",label=T)
DimPlot(chictaps, reduction = "umap",group.by = "mouse")
DimPlot(chictaps, reduction = "umap",group.by = "plate")
DimPlot(chictaps, reduction = "umap",group.by = "AP_axis")
##############mouse2##################

q1<-readRDS("~/Desktop/post-doc/Experiments/exp176-scEpi2-seq-mouse1-3-H3K27me3-prox-mid-distal/analysis/mm10_50kb_filtered.cov")
colnames(q1)<-gsub("sc50kb_mm10_support_JvB-176-SI-","",colnames(q1))

# prepare metadata
MD<-data.frame(matrix(0,ncol=3,nrow=ncol(q1)))
rownames(MD)<-colnames(q1)
colnames(MD)<-c("mouse","plate","AP_axis")

ss<-strsplit(rownames(MD),"-")
vec<-NULL
for (i in 1:length(ss)){
  MD[i,1]<-ss[[i]][1]
  vec[i]<-ss[[i]][2]
}
ss<-strsplit(rownames(MD),"_")
for (i in 1:length(ss)){
  MD[i,2]<-ss[[i]][1]
}

qqq<-read.csv('~/Desktop/post-doc/Experiments/exp176-scEpi2-seq-mouse1-3-H3K27me3-prox-mid-distal/analysis/chic_bc_well.tsv',sep='\t')
QQQ<-read.csv('~/Desktop/post-doc/Experiments/exp176-scEpi2-seq-mouse1-3-H3K27me3-prox-mid-distal/analysis/JvB-176-idx.tsv',sep='\t')
rownames(qqq)<-qqq[,3]

QQQ[,12]<-qqq[QQQ[,1],1]
QQQ[,13]<-paste(QQQ[,4],QQQ[,12],sep="_")
QQQ[,13]<-gsub("JvB-176-SI-","",QQQ[,13])
QQQ[,13]<-gsub("H3K27me3-scEpi2-seq-","",QQQ[,13])
rownames(QQQ)<-QQQ[,13]
QQ<-QQQ[rownames(MD),]
MD[,3]<-QQ$APaxis
MD<-MD[!is.na(MD[,3]),]
MD<-MD[MD[,1]=="mouse2",]
q1<-q1[,rownames(MD)]
rownames(q1)<-gsub("chr","",rownames(q1))


#fread("~/Desktop/post-doc/Experiments/exp176-scEpi2-seq-mouse1-3-H3K27me3-prox-mid-distal/analysis/merge/merge.bed.gz") %>%
#  distinct(V4)

chrom_assay <- Signac::CreateChromatinAssay(
  counts = q1,
  sep = c(":", "-"),
  fragments = "~/Desktop/post-doc/Experiments/exp176-scEpi2-seq-mouse1-3-H3K27me3-prox-mid-distal/analysis/mouse2.fragments.sorted.bed.gz",
  min.cells = 10,
  min.features = 200
)

chictaps <- CreateSeuratObject(
  counts = chrom_assay,
  assay = "peaks",
  meta.data = MD
)

peaks.keep <- seqnames(granges(chictaps)) %in% as.character(c(1:19))
chictaps <- chictaps[as.vector(peaks.keep), ]




ah <- AnnotationHub()
query(ah, "EnsDb.mmusculus.v99")
ensdb_v99 <- ah[["AH78811"]]

# extract gene annotations from EnsDb
annotations <- GetGRangesFromEnsDb(ensdb = ensdb_v99)
genome(annotations) <- "mm10"
Annotation(chictaps) <- annotations


# start Signac analysis
chictaps <- RunTFIDF(chictaps)
chictaps <- FindTopFeatures(chictaps, min.cutoff = 'q0')
chictaps <- RunSVD(chictaps)

# check for components that are strongly correlating with sequence depth
# for mouse3 I remove first component

DepthCor(chictaps)

chictaps <- RunUMAP(object = chictaps, reduction = 'lsi', dims = c(2:20))
chictaps <- FindNeighbors(object = chictaps, reduction = 'lsi', dims = c(2:20))
chictaps <- FindClusters(object = chictaps, verbose = FALSE, algorithm = 3, resolution = .1)
chictaps2 <- chictaps

DimPlot(chictaps2, reduction = "umap",label=T)
DimPlot(chictaps2, reduction = "umap",group.by = "mouse")
DimPlot(chictaps2, reduction = "umap",group.by = "plate")
DimPlot(chictaps2, reduction = "umap",group.by = "AP_axis")

##############mouse1##################

q1<-readRDS("~/Desktop/post-doc/Experiments/exp176-scEpi2-seq-mouse1-3-H3K27me3-prox-mid-distal/analysis/mm10_50kb_filtered.cov")
colnames(q1)<-gsub("sc50kb_mm10_support_JvB-176-SI-","",colnames(q1))

# prepare metadata
MD<-data.frame(matrix(0,ncol=3,nrow=ncol(q1)))
rownames(MD)<-colnames(q1)
colnames(MD)<-c("mouse","plate","AP_axis")

ss<-strsplit(rownames(MD),"-")
vec<-NULL
for (i in 1:length(ss)){
  MD[i,1]<-ss[[i]][1]
  vec[i]<-ss[[i]][2]
}
ss<-strsplit(rownames(MD),"_")
for (i in 1:length(ss)){
  MD[i,2]<-ss[[i]][1]
}

qqq<-read.csv('~/Desktop/post-doc/Experiments/exp176-scEpi2-seq-mouse1-3-H3K27me3-prox-mid-distal/analysis/chic_bc_well.tsv',sep='\t')
QQQ<-read.csv('~/Desktop/post-doc/Experiments/exp176-scEpi2-seq-mouse1-3-H3K27me3-prox-mid-distal/analysis/JvB-176-idx.tsv',sep='\t')
rownames(qqq)<-qqq[,3]

QQQ[,12]<-qqq[QQQ[,1],1]
QQQ[,13]<-paste(QQQ[,4],QQQ[,12],sep="_")
QQQ[,13]<-gsub("JvB-176-SI-","",QQQ[,13])
QQQ[,13]<-gsub("H3K27me3-scEpi2-seq-","",QQQ[,13])
rownames(QQQ)<-QQQ[,13]
QQ<-QQQ[rownames(MD),]
MD[,3]<-QQ$APaxis
MD<-MD[!is.na(MD[,3]),]
MD<-MD[MD[,1]=="mouse1",]
q1<-q1[,rownames(MD)]
rownames(q1)<-gsub("chr","",rownames(q1))


#fread("~/Desktop/post-doc/Experiments/exp176-scEpi2-seq-mouse1-3-H3K27me3-prox-mid-distal/analysis/merge/merge.bed.gz") %>%
#  distinct(V4)

chrom_assay <- Signac::CreateChromatinAssay(
  counts = q1,
  sep = c(":", "-"),
  fragments = "~/Desktop/post-doc/Experiments/exp176-scEpi2-seq-mouse1-3-H3K27me3-prox-mid-distal/analysis/mouse1.fragments.sorted.bed.gz",
  min.cells = 10,
  min.features = 200
)

chictaps <- CreateSeuratObject(
  counts = chrom_assay,
  assay = "peaks",
  meta.data = MD
)

peaks.keep <- seqnames(granges(chictaps)) %in% as.character(c(1:19))
chictaps <- chictaps[as.vector(peaks.keep), ]


ah <- AnnotationHub()
query(ah, "EnsDb.mmusculus.v99")
ensdb_v99 <- ah[["AH78811"]]

# extract gene annotations from EnsDb
annotations <- GetGRangesFromEnsDb(ensdb = ensdb_v99)
genome(annotations) <- "mm10"
Annotation(chictaps) <- annotations


# start Signac analysis
chictaps <- RunTFIDF(chictaps)
chictaps <- FindTopFeatures(chictaps, min.cutoff = 'q0')
chictaps <- RunSVD(chictaps)

# check for components that are strongly correlating with sequence depth
# for mouse3 I remove first component

DepthCor(chictaps)

chictaps <- RunUMAP(object = chictaps, reduction = 'lsi', dims = c(2:20))
chictaps <- FindNeighbors(object = chictaps, reduction = 'lsi', dims = c(2:20))
chictaps <- FindClusters(object = chictaps, verbose = FALSE, algorithm = 3, resolution = .1)
chictaps1 <- chictaps

DimPlot(chictaps1, reduction = "umap",label=T)
DimPlot(chictaps1, reduction = "umap",group.by = "mouse")
DimPlot(chictaps1, reduction = "umap",group.by = "plate")
DimPlot(chictaps1, reduction = "umap",group.by = "AP_axis")

#### merge mouse datasets ######

chictaps1$dataset <- 'mouse1'
chictaps2$dataset <- 'mouse2'
chictaps3$dataset <- 'mouse3'

combined <- merge(
  x = chictaps1,
  y = list(chictaps2, chictaps3),
  add.cell.ids = c("mouse1", "mouse2", "mouse3")
)

combined <- RunTFIDF(combined)
combined <- FindTopFeatures(combined, min.cutoff = )
combined <- RunSVD(combined)
combined <- RunUMAP(combined, dims = 2:20, reduction = 'lsi', min.dist = 0.1)
combined <- FindNeighbors(object = combined, reduction = 'lsi', dims = c(2:20))
combined <- FindClusters(object = combined, verbose = FALSE, algorithm = 3, resolution = .2)


DimPlot(combined, reduction = "umap",label=T)
DimPlot(combined, reduction = "umap",group.by = "mouse")
DimPlot(combined, reduction = "umap",group.by = "plate")
DimPlot(combined, reduction = "umap",group.by = "AP_axis")


cf <- intersect(rownames(chictaps1), rownames(chictaps2))

cf2 <- intersect(cf, rownames(chictaps3))

chictaps1 <-subset(chictaps1, features = cf2)
chictaps2 <-subset(chictaps2, features = cf2)
chictaps3 <-subset(chictaps3, features = cf2)

chictaps1 <- RunTFIDF(chictaps1)
chictaps1 <- FindTopFeatures(chictaps1, min.cutoff = 10)
chictaps1 <- RunSVD(chictaps1)

chictaps2 <- RunTFIDF(chictaps2)
chictaps2 <- FindTopFeatures(chictaps2, min.cutoff = 10)
chictaps2 <- RunSVD(chictaps2)

chictaps3 <- RunTFIDF(chictaps3)
chictaps3 <- FindTopFeatures(chictaps3, min.cutoff = 10)
chictaps3 <- RunSVD(chictaps3)

anchors <- FindIntegrationAnchors(
  object.list = list(chictaps1, chictaps2, chictaps3),
  anchor.features = 2000,
  reduction = "rlsi" 
)


chictaps.integrated <- IntegrateData(anchorset = anchors, dims = 1:20)


saveRDS(chictaps.integrated, file = "~/Desktop/post-doc/Experiments/exp176-scEpi2-seq-mouse1-3-H3K27me3-prox-mid-distal/analysis/merge/cchictaps.chr1.19.integrated.rds")



#chictaps.integrated <- readRDS("~/Desktop/post-doc/Experiments/exp176-scEpi2-seq-mouse1-3-H3K27me3-prox-mid-distal/analysis/merge/chictaps.integrated.rds")
DefaultAssay(chictaps.integrated) <- "peaks"

chictaps.integrated <- RunTFIDF(chictaps.integrated)
chictaps.integrated <- FindTopFeatures(chictaps.integrated, min.cutoff = 10)
chictaps.integrated <- RunSVD(chictaps.integrated) 
chictaps.integrated <- RunUMAP(chictaps.integrated, dims = 2:20, reduction = 'lsi', min.dist = 0.2, k.param=10)
chictaps.integrated <- FindNeighbors(object = chictaps.integrated, reduction = 'lsi', dims = c(2:20), k.param=9)
chictaps.integrated <- FindClusters(object = chictaps.integrated, verbose = FALSE, algorithm = 3, resolution = .1)

chictaps.integrated@meta.data <- chictaps.integrated@meta.data %>%
  mutate(sex = if_else(mouse == "mouse3", true = "female", false = "male"))

DimPlot(chictaps.integrated, reduction = "umap",label=T, shuffle=T)+scale_color_manual(values= c("yellow2", "orange2", "red2", "black"))

aa <- DimPlot(chictaps.integrated, reduction = "umap",label=T, shuffle=T)+scale_color_manual(values= c("yellow2","red2", "orange2"))
#DimPlot(chictaps.integrated, reduction = "umap",group.by = "mouse")
#DimPlot(chictaps.integrated, reduction = "umap",group.by = "plate")
bb <-DimPlot(chictaps.integrated, reduction = "umap",group.by = "AP_axis")+scale_color_manual(values= c("#5664BE", "#58BEBB", "#5AB83D" ))
cc <- DimPlot(chictaps.integrated, reduction = "umap",group.by = "sex", shuffle = T)+scale_color_manual(values= c("brown4", "black"))
#FeaturePlot(chictaps.integrated, features = "nCount_peaks", reduction = "umap")+scale_color_viridis_c( trans="log10", option="B")
#5AB83D, #58BEBB, #5664BE
aa|bb|cc

m3 <- fread("~/Desktop/post-doc/Experiments/exp176-scEpi2-seq-mouse1-3-H3K27me3-prox-mid-distal/analysis/mouse3.fragments.sorted.bed.gz" ) %>%
  mutate(total = sum(V5)) %>%
  group_by(V1) %>%
  summarise(n=(sum(V5)/total)*100) %>%
  ungroup() %>%
  distinct(V1, n) %>%
  arrange(-n)


m2 <- fread("~/Desktop/post-doc/Experiments/exp176-scEpi2-seq-mouse1-3-H3K27me3-prox-mid-distal/analysis/mouse2.fragments.sorted.bed.gz" ) %>%
  mutate(total = sum(V5)) %>%
  group_by(V1) %>%
  summarise(n=(sum(V5)/total)*100) %>%
  ungroup() %>%
  distinct(V1, n) %>%
  arrange(-n)

m1 <- fread("~/Desktop/post-doc/Experiments/exp176-scEpi2-seq-mouse1-3-H3K27me3-prox-mid-distal/analysis/mouse1.fragments.sorted.bed.gz" ) %>%
  mutate(total = sum(V5)) %>%
  group_by(V1) %>%
  summarise(n=(sum(V5)/total)*100) %>%
  ungroup() %>%
  distinct(V1, n) %>%
  arrange(-n)

m1 %>%
  dplyr::filter(V1 %in% c(1:19, "X", "Y")) %>%
  ggplot()+geom_col(aes(x=V1, y=n))+theme_bw()+
  ylab("reads (%)")+ggtitle("mouse1")+xlab("chromosomes")+
m2 %>%
  dplyr::filter(V1 %in% c(1:19, "X", "Y")) %>%
  ggplot()+geom_col(aes(x=V1, y=n))+
  ylab("reads (%)")+ggtitle("mouse2")+xlab("chromosomes")+theme_bw()+
m3 %>%
  dplyr::filter(V1 %in% c(1:19, "X", "Y")) %>%
  ggplot()+geom_col(aes(x=V1, y=n))+
  ylab("reads (%)")+ggtitle("mouse3")+xlab("chromosomes")+theme_bw()
  
m1 %>%
  dplyr::filter(V1 %in% c("X", "Y")) %>%
  ggplot()+geom_col(aes(x=V1, y=n))+theme_bw()+ylim(0,8)+
  ylab("reads (%)")+ggtitle("mouse1")+xlab("Sex chromosomes")+
  m2 %>%
  dplyr::filter(V1 %in% c("X", "Y")) %>%
  ggplot()+geom_col(aes(x=V1, y=n))+xlab("Sex chromosomes")+
  ylab("reads (%)")+ggtitle("mouse2")+theme_bw()+ylim(0,8)+
  m3 %>%
  dplyr::filter(V1 %in% c("X", "Y")) %>%
  ggplot()+geom_col(aes(x=V1, y=n))+xlab("Sex chromosomes")+ylim(0,8)+
  ylab("reads (%)")+ggtitle("mouse3")+theme_bw()




#find differentially digested genome locations
da_peaks <- FindAllMarkers(
  object = chictaps.integrated,
  test.use = 'wilcox',
  min.pct = 0.1
)

ok_down<-da_peaks$avg_log2FC < -1 & da_peaks$p_val_adj<0.05
ok_up<-da_peaks$avg_log2FC>1 & da_peaks$p_val_adj<0.05
da_peaks_down<-da_peaks[ok_down,]
da_peaks_up<-da_peaks[ok_up,]




#plot top 20 up and down digested locations
pdf("~/Desktop/post-doc/Experiments/exp176-scEpi2-seq-mouse1-3-H3K27me3-prox-mid-distal/analysis/merge/top20_up_in_proximal_cluster0vs1.pdf")
for (k in 1:20){
  q<-CoveragePlot(
  object = chictaps.integrated,
  region = rownames(da_peaks_up)[k],
  extend.upstream = 100000,
  extend.downstream = 100000,
  idents = c(0,1),
  fontsize = 6,
)
plot(q)
}
dev.off()

pdf("~/Desktop/post-doc/Experiments/exp176-scEpi2-seq-mouse1-3-H3K27me3-prox-mid-distal/analysis/merge/top20_down_in_proximal_cluster0vs1.pdf")
for (k in 1:20){
  q<-CoveragePlot(
    object = chictaps.integrated,
    region = rownames(da_peaks_down)[k],
    extend.upstream = 100000,
    extend.downstream = 100000,
    idents = c(0,1),
    fontsize = 6,
  )
  plot(q)
}
dev.off()


#####cluster0###

da_peaks1 <- FindMarkers(
  object = chictaps.integrated,
  ident.1 = 0,
  ident.2 = 2,
  test.use = 'wilcox',
  min.pct = 0.1
)
ok_down1<-da_peaks1$avg_log2FC < -1 & da_peaks1$p_val_adj<0.05
ok_up1<-da_peaks1$avg_log2FC>1 & da_peaks1$p_val_adj<0.05
da_peaks_down1<-da_peaks1[ok_down1,]
da_peaks_up1<-da_peaks1[ok_up1,]


#plot top 20 up and down digested locations
pdf("~/Desktop/post-doc/Experiments/exp176-scEpi2-seq-mouse1-3-H3K27me3-prox-mid-distal/analysis/merge/top20_up_in_proximal_cluster0V2.pdf")
for (k in 1:30){
  q<-CoveragePlot(window = 1500,
    object = chictaps.integrated,
    region = rownames(da_peaks_up1)[k],
    extend.upstream = 100000,
    extend.downstream = 100000,
    idents = c(0,1,2),
    fontsize = 6,
  )
  plot(q)+  scale_fill_manual(values= c("yellow2","red2", "orange2"))+
    scale_color_manual(values= c("yellow2","red2", "orange2"))+theme_bw()
}
dev.off()

pdf("~/Desktop/post-doc/Experiments/exp176-scEpi2-seq-mouse1-3-H3K27me3-prox-mid-distal/analysis/merge/top20_down_in_proximal_cluster0V2.pdf")
for (k in 1:30){
  q<-CoveragePlot(window = 1500,
    object = chictaps.integrated,
    region = rownames(da_peaks_down1)[k],
    extend.upstream = 100000,
    extend.downstream = 100000,
    idents = c(0,1,2),
    fontsize = 6,
  )
  plot(q)+  scale_fill_manual(values= c("yellow2","red2", "orange2"))+
    scale_color_manual(values= c("yellow2","red2", "orange2"))+theme_bw()
}
dev.off()


da_peaks1 <- FindMarkers(
  object = chictaps.integrated,
  ident.1 = 0,
  test.use = 'wilcox',
  min.pct = 0.1
)
ok_down1<-da_peaks1$avg_log2FC < -1 & da_peaks1$p_val_adj<0.05
ok_up1<-da_peaks1$avg_log2FC>1 & da_peaks1$p_val_adj<0.05
da_peaks_down1<-da_peaks1[ok_down1,]
da_peaks_up1<-da_peaks1[ok_up1,]


#plot top 20 up and down digested locations
pdf("~/Desktop/post-doc/Experiments/exp176-scEpi2-seq-mouse1-3-H3K27me3-prox-mid-distal/analysis/merge/top20_up_in_proximal_cluster0.pdf")
for (k in 1:20){
  q<-CoveragePlot(
    object = chictaps.integrated,
    region = rownames(da_peaks_up1)[k],
    extend.upstream = 100000,
    extend.downstream = 100000,
    idents = c(0,1,2,3,4,5),
    fontsize = 6,
  )
  plot(q)
}
dev.off()

pdf("~/Desktop/post-doc/Experiments/exp176-scEpi2-seq-mouse1-3-H3K27me3-prox-mid-distal/analysis/merge/top20_down_in_proximal_cluster0.pdf")
for (k in 1:20){
  q<-CoveragePlot(
    object = chictaps.integrated,
    region = rownames(da_peaks_down1)[k],
    extend.upstream = 100000,
    extend.downstream = 100000,
    idents = c(0,1,2,3,4,5),
    fontsize = 6,
  )
  plot(q)
}
dev.off()


#####cluster1###

da_peaks1 <- FindMarkers(
  object = chictaps.integrated,
  ident.1 = 2,
  test.use = 'wilcox',
  min.pct = 0.1
)
ok_down1<-da_peaks1$avg_log2FC < -1 & da_peaks1$p_val_adj<0.05
ok_up1<-da_peaks1$avg_log2FC>1 & da_peaks1$p_val_adj<0.05
da_peaks_down1<-da_peaks1[ok_down1,]
da_peaks_up1<-da_peaks1[ok_up1,]


#plot top 20 up and down digested locations
pdf("~/Desktop/post-doc/Experiments/exp176-scEpi2-seq-mouse1-3-H3K27me3-prox-mid-distal/analysis/merge/top20_up_in_proximal_cluster2.pdf")
for (k in 1:10){
  q<-CoveragePlot(
    object = chictaps.integrated,
    region = rownames(da_peaks_up1)[k],
    extend.upstream = 100000,
    extend.downstream = 100000,
    idents = c(0,1,2,3,4,5),
    fontsize = 6,
  )
  plot(q)
}
dev.off()

pdf("~/Desktop/post-doc/Experiments/exp176-scEpi2-seq-mouse1-3-H3K27me3-prox-mid-distal/analysis/merge/top20_down_in_proximal_cluster2.pdf")
for (k in 1:10){
  q<-CoveragePlot(
    object = chictaps.integrated,
    region = rownames(da_peaks_down1)[k],
    extend.upstream = 100000,
    extend.downstream = 100000,
    idents = c(0,1,2,3,4,5),
    fontsize = 6,
  )
  plot(q)
}
dev.off()




da_peaks_all <- FindMarkers(
  object = chictaps.integrated,
  test.use = 'wilcox',
  min.pct = 0.1
)

ok_down1<-da_peaks1$avg_log2FC < -1 & da_peaks1$p_val_adj<0.05
ok_up1<-da_peaks1$avg_log2FC>1 & da_peaks1$p_val_adj<0.05
da_peaks_down1<-da_peaks1[ok_down1,]
da_peaks_up1<-da_peaks1[ok_up1,]


chictaps.integrated <- readRDS("~/Desktop/post-doc/Experiments/exp176-scEpi2-seq-mouse1-3-H3K27me3-prox-mid-distal/analysis/merge/cchictaps.chr1.19.integrated.rds")
all_markers <- FindAllMarkers(
  object = chictaps.integrated,
  test.use = 'wilcox',
  min.pct = 0.1,      
  logfc.threshold = 0.2  
) 

all_markers1 <- all_markers%>% separate(gene, into = c("chr", "start", "end"), sep = "[-]") %>% 
  mutate(start = as.numeric(start),
         end = as.numeric(end)) %>% as.data.table()



# Read the GTF file
gtf <- fread("~/Desktop/post-doc/Experiments/exp176-scEpi2-seq-mouse1-3-H3K27me3-prox-mid-distal/analysis/gencode.vM10.basic.annotation.gtf.gz", sep = "\t", header = FALSE, quote = "")

# Assign standard GTF column names
colnames(gtf) <- c("seqname", "source", "feature", "start", "end", 
                   "score", "strand", "frame", "attributes")

# View the first few rows
head(gtf)
gtf %>%
  group_by(feature) %>%
  summarize(n=n())

fread("~/Desktop/post-doc/Experiments/exp176-scEpi2-seq-mouse1-3-H3K27me3-prox-mid-distal/analysis/gencode.vM10.basic.annotation.gtf.gz") 

gtf_parsed <- gtf %>%
  mutate(
    gene_name = str_extract(attributes, 'gene_name "[^"]+"') %>% str_remove_all('gene_name |"'),
    gene_id = str_extract(attributes, 'gene_id "[^"]+"') %>% str_remove_all('transcript_id |"'),
    chr = gsub("chr", "", seqname),
    start = as.numeric(start),
    end = as.numeric(end)
  ) %>% dplyr::filter(feature == "gene") %>%
  dplyr::select(-attributes, -seqname) %>% as.data.table() # Remove original attributes column if needed

head(gtf_parsed)


setkeyv(x = all_markers1, c('chr', 'start', 'end'))

setkeyv(x = gtf_parsed,c('chr', 'start', 'end'))

foverlaps(all_markers1, gtf_parsed) %>% 
  dplyr::filter(gene_name == "Lgr5" #&
                #p_val_adj <0.05 & 
                #(avg_log2FC < -0.75 | avg_log2FC > 0.75)
                ) %>%
  write_tsv("~/Desktop/post-doc/Experiments/exp176-scEpi2-seq-mouse1-3-H3K27me3-prox-mid-distal/analysis/merge/lgr5.txt")

library(readxl)
library(tidyverse)

# Define the Excel file path
file_path <- "~/Desktop/post-doc/Experiments/exp176-scEpi2-seq-mouse1-3-H3K27me3-prox-mid-distal/analysis/merge/Grun.xlsx"

# Get all sheet names
sheets <- excel_sheets(file_path)

# Read all sheets and combine them with bind_rows()
combined_df <- bind_rows(
  lapply(sheets, function(sheet) {
    read_excel(file_path, sheet = sheet) %>%
      mutate(SheetName = sheet)  # Add sheet name as a new column
  })
)

# View combined data
head(combined_df)
genes <- combined_df %>%
  dplyr::filter(pv <0.05 & fc > 6) %>%
  separate(GENEID, into = c("gene_name", "chr"), sep = "[_]") %>% 
  distinct(gene_name) %>%
  select(gene_name) %>%
  pull()


 Intestine <- foverlaps(all_markers1, gtf_parsed) %>% 
   dplyr::filter(  p_val_adj <0.05)  %>% 
   dplyr::filter(gene_name %in% c("Tead1", "Tcf4", "Lgr5", "Kank1", "Ehf", "Rbfox2", "Ephb3", "Hopx", "Lrig", "Prom1",
                                                   "Pou2f3", "Sox9", "Dclk1","Xbp1", "Nfib",#,
                                  "Foxa1", "Foxa2" , "Zfhx3",  "Hunk", "Bmpr1a", "Notch2", "Bmp7"
   )
   )%>%
   ggplot()+
   geom_boxplot(aes(x=cluster, y=avg_log2FC), outlier.shape=NA)+
   geom_jitter(aes(x=cluster, y=avg_log2FC, col=-log10(p_val_adj)), size=2)+
   geom_hline(yintercept = 0, linetype=2, color="blue")+
   scale_color_viridis_c(option="H",limits = c(0,15), oob = scales::oob_squish)+ggtitle("H3K27me3 on \n Intestinal Genes")+
 theme_bw()+theme(legend.position = "none")+
   xlab("Cluster")+
   ylab("average log2FC")
   
 
Enterocrine <- foverlaps(all_markers1, gtf_parsed) %>%  
   
    dplyr::filter(gene_name %in% c(
      "Nkx2-2", 
      "Lhx5",
      "Lhx1", 
      "Nkx2-2os",
      "Nptx1", 
      "Gata6",
      "Shc2",
      "Lmx1b",
      "Fev", 
      "Cdk5r2",
      "Pax6os1"
 )) %>%dplyr::filter(p_val_adj <0.05) %>%
  ggplot()+geom_boxplot(aes(x=cluster, y=avg_log2FC), outlier.shape=NA)+
  geom_jitter(aes(x=cluster, y=avg_log2FC, col=-log10(p_val_adj)), size=2)+
  geom_hline(yintercept = 0, linetype=2, color="blue")+
  scale_color_viridis_c(option="H",limits = c(0,15), oob = scales::oob_squish)+ggtitle("H3K27me3 on \n Secretory Genes")+
  theme_bw()+
  ylab("average log2FC")+
  xlab("Cluster")

Immune <- foverlaps(all_markers1, gtf_parsed) %>%  
  
  dplyr::filter(gene_name %in% c(
"Enthd1",
"Gli2",
"Gm44536",
"Hoxb5os",
"Hoxb6",
"Ikzf3",
"Kcnb1",
"Ntn1",
"Pax5",
"Plxnc1",
"Pou2f2",
"Rcsd1", "Rasal1",
"Pax5", "Skap1", "Ikzf3", "Hoxb6", "Pde2a")
)%>%
  ggplot()+
  geom_boxplot(aes(x=cluster, y=avg_log2FC), outlier.shape=NA)+
  geom_jitter(aes(x=cluster, y=avg_log2FC, col=-log10(p_val_adj)), size=2)+
  #geom_boxplot(aes(x=cluster, y=avg_log2FC), outlier.shape=NA)+
  geom_hline(yintercept = 0, linetype=2, color="blue")+
  scale_color_viridis_c(option="H",limits = c(0,15), oob = scales::oob_squish)+ggtitle("H3K27me3 on \n Immune Genes")+
  theme_bw()+theme(legend.position = "none")+
  xlab("Cluster")+
  ylab("average log2FC")


Intestine|Immune|Enterocrine
aa
all_markers1 %>%
  ggplot()+
  geom_point(aes(x= avg_log2FC, y=-log10(p_val_adj), col = cluster), alpha=0.5, size=0.5)+
  theme_bw()+
  facet_grid(cols=vars(cluster))+
  scale_color_manual(values = c("yellow3", "orange2", "red2"))
  
 
xa <- gtf_parsed %>%
  dplyr::filter((gene_name %in% c(  "Nkx2-2","Fev", "Shc2", "Gata6", "Gata4", "Cdk5r2")))  %>% 
  mutate(pos = paste0(chr, "-", start, "-", end))
da_peaks_up1 <- xa$pos

pdf("~/Desktop/post-doc/Experiments/exp176-scEpi2-seq-mouse1-3-H3K27me3-prox-mid-distal/analysis/merge/Enteroendocrine_genes.pdf")
for (k in 1:6){
  q<-CoveragePlot(
    object = chictaps.integrated,
    region = (da_peaks_up1)[k],
    extend.upstream = 300000,
    extend.downstream = 300000,
    idents = c(0,1,2),
    fontsize = 6,
  )
  plot(q)
}
dev.off()


chictaps_split <- SplitObject(chictaps.integrated, split.by = "seurat_clusters")

intestinalcells <- chictaps_split$`0`
immunecells <- chictaps_split$`1`
EEcells <- chictaps_split$`2`

intestinalcells <- intestinalcells %>%
                      RunTFIDF() %>%
                      FindTopFeatures(min.cutoff = 20) %>%
                      RunSVD() %>%
                      RunUMAP(dims= 2:20, reduction = "lsi", min.dist = 0.2, k.param=5) %>%
                      FindNeighbors(reduction = 'lsi', dims = c(2:20), k.param=5) %>%
                      FindClusters(verbose = FALSE, algorithm = 3, resolution = .2)
  


aaa <- DimPlot(intestinalcells, reduction = "umap",label=T, shuffle=T)+scale_color_manual(values= c("yellow2","yellow4"))+ggtitle("Intestinal Cells")
#DimPlot(chictaps.integrated, reduction = "umap",group.by = "mouse")
#DimPlot(chictaps.integrated, reduction = "umap",group.by = "plate")
bbb <-DimPlot(intestinalcells, reduction = "umap",group.by = "AP_axis",shuffle = T)+scale_color_manual(values= c("magenta", "blue", "cyan"))
ccc <- DimPlot(intestinalcells, reduction = "umap",group.by = "sex", shuffle = T)+scale_color_manual(values= c("brown4", "black"))
#FeaturePlot(chictaps.integrated, features = "nCount_peaks", reduction = "umap")+scale_color_viridis_c( trans="log10", option="B")

aaa|bbb|ccc


da_peaks1 <- FindAllMarkers(intestinalcells,logfc.threshold = 3, test.use = "wilcox")

ok_down1<-da_peaks1$avg_log2FC < -1 & da_peaks1$p_val_adj<0.05
ok_up1<-da_peaks1$avg_log2FC>1 & da_peaks1$p_val_adj<0.05
da_peaks_down1<-da_peaks1[ok_down1,]
da_peaks_up1<-da_peaks1[ok_up1,]


#plot top 20 up and down digested locations
pdf("~/Desktop/post-doc/Experiments/exp176-scEpi2-seq-mouse1-3-H3K27me3-prox-mid-distal/analysis/merge/top20_up_in_Int.pdf")
for (k in 1:20){
  q<-CoveragePlot(window = 5000,
                  object = intestinalcells,
                  region = rownames(da_peaks_up1)[k],
                  extend.upstream = 100000,
                  extend.downstream = 100000,
                  idents = c(0,1),
                  fontsize = 6,
  )
  plot(q)+  scale_fill_manual(values= c("yellow2","red2", "orange2"))+
    scale_color_manual(values= c("yellow2","red2", "orange2"))+theme_bw()
}
dev.off()

pdf("~/Desktop/post-doc/Experiments/exp176-scEpi2-seq-mouse1-3-H3K27me3-prox-mid-distal/analysis/merge/top20_down_in_Int.pdf")
for (k in 1:15){
  q<-CoveragePlot(window = 5000,
                  object = intestinalcells,
                  region = rownames(da_peaks_down1)[k],
                  extend.upstream = 100000,
                  extend.downstream = 100000,
                  idents = c(0,1),
                  fontsize = 6,
  )
  plot(q)+  scale_fill_manual(values= c("yellow2","red2", "orange2"))+
    scale_color_manual(values= c("yellow2","red2", "orange2"))+theme_bw()
}
dev.off()









immunecells <- immunecells %>%
  RunTFIDF() %>%
  FindTopFeatures(min.cutoff = 20) %>%
  RunSVD() %>%
  RunUMAP(dims= 3:20, reduction = "lsi", min.dist = 0.2, k.param=7) %>%
  FindNeighbors(reduction = 'lsi', dims = c(3:20), k.param=10) %>%
  FindClusters(verbose = FALSE, algorithm = 3, resolution = .2)




aaaa <- DimPlot(immunecells, reduction = "umap",label=T, shuffle=T)+scale_color_manual(values= c("red1","red4", "maroon4", "black"))+ggtitle("Immune Cells")
#DimPlot(chictaps.integrated, reduction = "umap",group.by = "mouse")
#DimPlot(chictaps.integrated, reduction = "umap",group.by = "plate")
bbbb <-DimPlot(immunecells, reduction = "umap",group.by = "AP_axis",shuffle = T)+scale_color_manual(values= c("magenta", "blue", "cyan"))
cccc <- DimPlot(immunecells, reduction = "umap",group.by = "sex", shuffle = T)+scale_color_manual(values= c("brown4", "black"))
#FeaturePlot(chictaps.integrated, features = "nCount_peaks", reduction = "umap")+scale_color_viridis_c( trans="log10", option="B")

aaaa|bbbb|cccc

da_peaks1 <- FindMarkers(immunecells,logfc.threshold = 3, test.use = "wilcox")

ok_down1<-da_peaks1$avg_log2FC < -1 & da_peaks1$p_val_adj<0.05
ok_up1<-da_peaks1$avg_log2FC>1 & da_peaks1$p_val_adj<0.05
da_peaks_down1<-da_peaks1[ok_down1,]
da_peaks_up1<-da_peaks1[ok_up1,]

da_peaks1 <- FindMarkers(immunecells,logfc.threshold = 2, ident.1 = 0, ident.2 = 2, test.use = "wilcox")
ok_down1<-da_peaks1$avg_log2FC < -1 & da_peaks1$p_val_adj<0.05
ok_up1<-da_peaks1$avg_log2FC>1 & da_peaks1$p_val_adj<0.05
da_peaks_down1<-da_peaks1[ok_down1,]
da_peaks_up1<-da_peaks1[ok_up1,]

#plot top 20 up and down digested locations
pdf("~/Desktop/post-doc/Experiments/exp176-scEpi2-seq-mouse1-3-H3K27me3-prox-mid-distal/analysis/merge/up_in_Immune_cl0_cl2.pdf")
for (k in 1:10){
  q<-CoveragePlot(window = 5000,
                  object = immunecells,
                  region = rownames(da_peaks_up1)[k],
                  extend.upstream = 100000,
                  extend.downstream = 100000,
                  idents = c(0,2),
                  fontsize = 6,
  )
  plot(q)+  scale_fill_manual(values= c("yellow2","red2", "orange2"))+
    scale_color_manual(values= c("yellow2","red2", "orange2"))+theme_bw()
}
dev.off()

pdf("~/Desktop/post-doc/Experiments/exp176-scEpi2-seq-mouse1-3-H3K27me3-prox-mid-distal/analysis/merge/down_in_Immune_cl0_cl2.pdf")
for (k in 1:10){
  q<-CoveragePlot(window = 5000,
                  object = immunecells,
                  region = rownames(da_peaks_down1)[k],
                  extend.upstream = 100000,
                  extend.downstream = 100000,
                  idents = c(0,1,2),
                  fontsize = 6,
  )
  plot(q)+  scale_fill_manual(values= c("yellow2","red2", "orange2"))+
    scale_color_manual(values= c("yellow2","red2", "orange2"))+theme_bw()
}
dev.off()









EEcells <- EEcells %>%
  RunTFIDF() %>%
  FindTopFeatures(min.cutoff = 20) %>%
  RunSVD() %>%
  RunUMAP(dims= 2:20, reduction = "lsi", min.dist = 0.2, k.param=10) %>%
  FindNeighbors(reduction = 'lsi', dims = c(2:20), k.param=10) %>%
  FindClusters(verbose = FALSE, algorithm = 3, resolution = .1)



aaaaa <- DimPlot(EEcells, reduction = "umap",label=T, shuffle=T, pt.size = 2)+ggtitle("Secretory cells")+scale_color_manual(values= c("orange1", "orange4"))
#DimPlot(chictaps.integrated, reduction = "umap",group.by = "mouse")
#DimPlot(chictaps.integrated, reduction = "umap",group.by = "plate")
bbbbb <-DimPlot(EEcells, reduction = "umap",group.by = "AP_axis",shuffle = T, pt.size = 2)+scale_color_manual(values= c("magenta", "blue", "cyan"))
ccccc <- DimPlot(EEcells, reduction = "umap",group.by = "sex", shuffle = T, pt.size = 2)+scale_color_manual(values= c("brown4", "black"))
#FeaturePlot(chictaps.integrated, features = "nCount_peaks", reduction = "umap")+scale_color_viridis_c( trans="log10", option="B")

aaaaa|bbbbb|ccccc

da_peaks1 <- FindAllMarkers(EEcells,logfc.threshold = 1, test.use = "wilcox")

ok_down1<-da_peaks1$avg_log2FC < -1 & da_peaks1$p_val_adj<0.05
ok_up1<-da_peaks1$avg_log2FC>1 & da_peaks1$p_val_adj<0.05
da_peaks_down1<-da_peaks1[ok_down1,]
da_peaks_up1<-da_peaks1[ok_up1,]


#plot top 20 up and down digested locations
pdf("~/Desktop/post-doc/Experiments/exp176-scEpi2-seq-mouse1-3-H3K27me3-prox-mid-distal/analysis/merge/top20_up_in_EEC.pdf")
for (k in 1:15){
  q<-CoveragePlot(window = 5000,
                  object = EEcells,
                  region = rownames(da_peaks_up1)[k],
                  extend.upstream = 100000,
                  extend.downstream = 100000,
                  idents = c(0,1),
                  fontsize = 6,
  )
  plot(q)+  scale_fill_manual(values= c("yellow2","red2", "orange2"))+
    scale_color_manual(values= c("yellow2","red2", "orange2"))+theme_bw()
}
dev.off()

pdf("~/Desktop/post-doc/Experiments/exp176-scEpi2-seq-mouse1-3-H3K27me3-prox-mid-distal/analysis/merge/top20_down_in_EEC.pdf")
for (k in 1:15){
  q<-CoveragePlot(window = 5000,
                  object = EEcells,
                  region = rownames(da_peaks_down1)[k],
                  extend.upstream = 100000,
                  extend.downstream = 100000,
                  idents = c(0,1),
                  fontsize = 6,
  )
  plot(q)+  scale_fill_manual(values= c("yellow2","red2", "orange2"))+
    scale_color_manual(values= c("yellow2","red2", "orange2"))+theme_bw()
}
dev.off()


(aa|bb|cc)/(aaa|bbb|ccc)/(aaaa|bbbb|cccc)/(aaaaa|bbbbb|ccccc)+plot_layout(heights = c(3,2,2,1))

bind_rows(
  
  (immunecells@meta.data %>% mutate(supercluster = "Immune",
                                 subcluster  = case_when(
                                   seurat_clusters == 0 ~ "B-cell",
                                   seurat_clusters == 1 ~ "T-cell",
                                   seurat_clusters == 2 ~ "Monocyte",
                                   TRUE ~ as.character(seurat_clusters)  # Default case (if other values exist)
                                 ) )),

(intestinalcells@meta.data %>% mutate(supercluster = "Intestinal",
                                 subcluster  = case_when(
                                   seurat_clusters == 0 ~ "Intestinal1",
                                   seurat_clusters == 1 ~ "Intestinal2",
                                   TRUE ~ as.character(seurat_clusters)  # Default case (if other values exist)
                                 ) )),

(EEcells@meta.data %>% mutate(supercluster = "Secretory",
                                     subcluster  = case_when(
                                       seurat_clusters == 0 ~ "Goblet",
                                       seurat_clusters == 1 ~ "EECs",
                                       TRUE ~ as.character(seurat_clusters)  # Default case (if other values exist)
                                     ) ))

) %>% group_by(subcluster) %>%
  summarise(n=n()) %>%
  ggplot()+geom_col(aes(x=subcluster, y=n))



bind_rows(
  
  (immunecells@meta.data %>% mutate(supercluster = "Immune",
                                    subcluster  = case_when(
                                      seurat_clusters == 0 ~ "B-cell",
                                      seurat_clusters == 1 ~ "T-cell",
                                      seurat_clusters == 2 ~ "Monocyte",
                                      TRUE ~ as.character(seurat_clusters)  # Default case (if other values exist)
                                    ) )),
  
  (intestinalcells@meta.data %>% mutate(supercluster = "Intestinal",
                                        subcluster  = case_when(
                                          seurat_clusters == 0 ~ "Intestinal1",
                                          seurat_clusters == 1 ~ "Intestinal2",
                                          TRUE ~ as.character(seurat_clusters)  # Default case (if other values exist)
                                        ) )),
  
  (EEcells@meta.data %>% mutate(supercluster = "Secretory",
                                subcluster  = case_when(
                                  seurat_clusters == 0 ~ "Goblet",
                                  seurat_clusters == 1 ~ "EECs",
                                  TRUE ~ as.character(seurat_clusters)  # Default case (if other values exist)
                                ) ))
  
) %>% group_by(supercluster) %>%
  summarise(n=n()) %>%
  ggplot()+geom_col(aes(x=supercluster, y=n))

MD2 <- bind_rows(
  
  (immunecells@meta.data %>% mutate(supercluster = "Immune",
                                    subcluster  = case_when(
                                      seurat_clusters == 0 ~ "B-cell",
                                      seurat_clusters == 1 ~ "T-cell",
                                      seurat_clusters == 2 ~ "Monocyte",
                                      TRUE ~ as.character(seurat_clusters)  # Default case (if other values exist)
                                    ) )),
  
  (intestinalcells@meta.data %>% mutate(supercluster = "Intestinal",
                                        subcluster  = case_when(
                                          seurat_clusters == 0 ~ "Intestinal1",
                                          seurat_clusters == 1 ~ "Intestinal2",
                                          TRUE ~ as.character(seurat_clusters)  # Default case (if other values exist)
                                        ) )),
  
  (EEcells@meta.data %>% mutate(supercluster = "Secretory",
                                subcluster  = case_when(
                                  seurat_clusters == 0 ~ "Goblet",
                                  seurat_clusters == 1 ~ "EECs",
                                  TRUE ~ as.character(seurat_clusters)  # Default case (if other values exist)
                                ) ))
  
) 

chictaps.integrated <- AddMetaData(chictaps.integrated,metadata = MD2)


a <- DimPlot(chictaps.integrated, reduction = "umap",group.by = "supercluster", repel = T, label=T)+scale_color_manual(values= c("red2", "yellow2", "orange2"))
b <- DimPlot(chictaps.integrated, reduction = "umap",group.by = "subcluster", repel = T, label=T)+scale_color_manual(values= c("red", "brown","orange1", "yellow2","yellow4",  "red4", "maroon" ))

a|b

#DimPlot(chictaps.integrated, reduction = "umap",group.by = "mouse")
#DimPlot(chictaps.integrated, reduction = "umap",group.by = "plate")
bb <-DimPlot(chictaps.integrated, reduction = "umap",group.by = "AP_axis",shuffle = T)+scale_color_manual(values= c("magenta", "blue", "cyan"))
cc <- DimPlot(chictaps.integrated, reduction = "umap",group.by = "sex", shuffle = T)+scale_color_manual(values= c("brown4", "black"))
#FeaturePlot(chictaps.integrated, features = "nCount_peaks", reduction = "umap")+scale_color_viridis_c( trans="log10", option="B")

bb|cc

###countspercelll###



chictaps.integrated@meta.data %>%
  ggplot()+
  geom_density(aes( x= nCount_peaks, fill=seurat_clusters),alpha=0.3, adjust = 1/4, position = position_identity())+
  scale_color_manual(values= c("red2", "yellow2", "orange2"))+
  scale_fill_manual(values= c("red2", "yellow2", "orange2"))+
  scale_x_log10()+theme_bw()+
  facet_grid(rows=vars(seurat_clusters))

chictaps.integrated@meta.data %>%
  ggplot()+
  geom_density(aes( x= nCount_peaks, fill=seurat_clusters),alpha=0.3, adjust = 1/4, position = position_identity())+
  scale_color_manual(values= c("red2", "yellow2", "orange2"))+
  scale_fill_manual(values= c("red2", "yellow2", "orange2"))+
  scale_x_log10()+theme_bw()+
  facet_grid(rows=vars(seurat_clusters))


chictaps.integrated@meta.data %>%
  ggplot()+
  geom_density(aes( x= nCount_peaks, fill=mouse),alpha=0.3, adjust = 1/4, position = position_identity())+
  scale_fill_manual(values= c("blue","cyan", "green2"))+
  scale_x_log10()+theme_bw()+
  facet_grid(rows=vars(mouse))

chictaps.integrated@meta.data %>%
  group_by(mouse) %>%
  #mutate(nCount_peaks = scale(nCount_peaks)) %>%
  ggplot()+
  geom_violin(aes( y= nCount_peaks, x=paste0(supercluster), fill=supercluster),outlier.shape = NA,  alpha=0.3)+
  geom_boxplot(aes( y= nCount_peaks, x=paste0(supercluster), fill=supercluster),outlier.shape = NA, width =0.25,  alpha=0.3)+
  geom_jitter(aes( y= nCount_peaks, x=paste0(supercluster), col=supercluster), width =0.2, alpha=0.3)+
  scale_color_manual(values= c("red2", "yellow2", "orange2"))+
  scale_fill_manual(values= c("red2", "yellow2", "orange2"))+theme_bw()+
  scale_y_log10()+#facet_grid(cols=vars(mouse))+
  ylab("reads per cell (log10")+
  xlab("")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) 


chictaps.integrated@meta.data %>%
  group_by(mouse) %>%
  #mutate(nCount_peaks = scale(nCount_peaks)) %>%
  ggplot()+
  geom_violin(aes( y= nCount_peaks, x=paste0(subcluster), fill=subcluster),outlier.shape = NA,  alpha=0.3, trim=F)+
  geom_boxplot(aes( y= nCount_peaks, x=paste0(subcluster), fill=subcluster),outlier.shape = NA, width =0.25,  alpha=0.3)+
  geom_jitter(aes( y= nCount_peaks, x=paste0(subcluster), col=subcluster), width =0.2, alpha=0.2, size=0.75)+
  scale_color_manual(values= c("red", "brown","orange1", "yellow2","yellow4",  "red4", "maroon" ))+
  scale_fill_manual(values= c("red", "brown","orange1", "yellow2","yellow4",  "red4", "maroon" ))+
  scale_y_log10()+#facet_grid(rows=vars(supercluster))+
  ylab("reads per cell (log10")+
  xlab("")+theme_bw()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+theme(legend.position = "none")


chictaps.integrated <- chictaps.combined
chictaps.integrated@meta.data %>%
  group_by(mouse) %>%
  #mutate(nCount_peaks = scale(nCount_peaks)) %>%
  ggplot()+
  geom_violin(aes( y= nCount_peaks, x=factor(subcluster, order), fill=factor(subcluster, order)),outlier.shape = NA,  alpha=0.3, trim=F)+
  geom_boxplot(aes( y= nCount_peaks, x=factor(subcluster, order), fill=factor(subcluster, order)),outlier.shape = NA, width =0.25,  alpha=0.3)+
  geom_jitter(aes( y= nCount_peaks, x=factor(subcluster, order), col=factor(subcluster, order)), width =0.2, alpha=0.2, size=0.75)+
  scale_color_manual(values= c( "yellow2","yellow4", "orange1","brown", "red", "red4", "maroon" ))+
  scale_fill_manual(values= c( "yellow2","yellow4", "orange1","brown", "red", "red4", "maroon" ))+
  scale_y_log10()+#facet_grid(rows=vars(supercluster))+
  ylab("reads per cell (log10)")+theme_bw()+theme(legend.position = "bottom")+
  xlab("")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

order = c("Intestinal1",  "Intestinal2", "Goblet", "EECs", "B-cell", "T-cell", "Monocyte")
chictaps.integrated@meta.data %>%
  ggplot()+
  geom_boxplot(aes( y= nCount_peaks, x=seurat_clusters, fill=seurat_clusters),width=0.5, alpha=0.3)+
  geom_jitter(aes( y= nCount_peaks, x=seurat_clusters, col=seurat_clusters), width =0.2, alpha=0.3)+
  scale_fill_manual(values= c("yellow2","red2", "orange2"))+
  scale_color_manual(values= c("yellow2","red2", "orange2"))+
  scale_y_log10()+theme_bw()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) 

chictaps.integrated@meta.data %>%
  count(seurat_clusters, AP_axis, mouse) %>%
  group_by(seurat_clusters) %>%
  mutate(proportion = n / sum(n)) %>%
  ggplot()+
  geom_col(aes(x=seurat_clusters, y=proportion, fill=AP_axis))+
  scale_fill_manual(values= c("magenta", "blue", "cyan"))+
  scale_color_manual(values= c("magenta", "blue", "cyan"))+theme_bw()

chictaps.integrated@meta.data %>%
  ggplot()+
  geom_boxplot(aes( y= nCount_peaks, x=seurat_clusters, fill=seurat_clusters),width=0.5, alpha=0.3)+
  geom_jitter(aes( y= nCount_peaks, x=seurat_clusters, col=seurat_clusters), width =0.2, alpha=0.3)+
  scale_fill_manual(values= c("yellow2","red2", "orange2"))+
  scale_color_manual(values= c("yellow2","red2", "orange2"))+
  facet_grid(cols=vars(sex))+
  scale_y_log10()+theme_bw()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) 
