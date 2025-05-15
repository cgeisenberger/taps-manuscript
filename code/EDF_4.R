
library(tidyverse)
library(data.table)
library(Rcpp)
library(patchwork)


fread('~/Desktop/post-doc/ChiC-TAPS/Fig2/all_cells_sum.tsv')%>%
  mutate(mod = str_split(V1, "-", simplify = T)[, 6],
         pl = str_split(V1, "-", simplify = T)[, 7],
         pl = str_split(pl, "_", simplify = T)[, 1])%>% dplyr::filter(mod =='k27me3')%>%
  ggplot()+geom_point(aes(x=log10(N), y=fme, col=mod),size =.5, alpha=0.8)+
  facet_grid(cols= vars(mod), rows = vars(pl), space = 'free')+
  scale_color_manual(values=c("#D674BA","#86BF88","#7A416A"))+theme_bw()+theme(legend.position = 'none')+
  ylab('CpG methylation per cell')+xlab('unique ChIC reads per cell')+ggtitle('RPE-1 hTERT')+
  geom_vline(xintercept = 4, linetype=2)+
 # geom_hline(yintercept = 0.65, linetype=2)+
  geom_hline(yintercept = 0.45, linetype=2)
  

filter<-fread('~/Desktop/post-doc/ChiC-TAPS/Fig2/all_cells_sum.tsv')%>%
  mutate(mod = str_split(V1, "-", simplify = T)[, 6])%>%
  group_split(mod)
  
filt_cells<-list(
     filter[[1]]%>% dplyr::filter(N>1e4, fme<0.45),
     filter[[2]]%>% dplyr::filter(N>1e4, fme>0.65),
     filter[[3]]%>% dplyr::filter(N>1e4, fme<0.45))%>%
  rbindlist()
     

encode <- readRDS('~/Desktop/post-doc/ChiC-TAPS/ExtendedDataFigures/Supp2/rpe-jvl.rds')%>%
  #dplyr::filter(cov>0)%>% 
  select(chr, start, cov, beta)%>% dplyr::filter(cov>1)%>%
  group_by(chr)%>%
  mutate(roundbp = round(start/1e4)*1e4)%>% group_by(chr, roundbp) %>% summarize(me = mean(beta))
encode <- encode%>% ungroup()%>% mutate(chr = str_replace(chr, "chr", ""), location = paste0(chr, ":", roundbp))

meth_cuts_reads <- readRDS("~/Desktop/post-doc/ChiC-TAPS/Fig2/meth_cuts_reads.RDS")#%>% dplyr::filter(filt_cells$V1 %in% V1)#[[1]][, max(bp) - min(bp), .(V1, chr, bin)]

rbindlist(meth_cuts_reads)%>% summarize(n= length(unique(V1)))


# which bin in which mod
mod_bin <- rbindlist(map(meth_cuts_reads, function(x) {
  
  x[V1 %in% filt_cells$V1
  ][, .(sumtot = sum(tot), summe=sum(me)), .(chr, bin = (round(bp / 1e5)*1e5),  mod = str_extract(V1, "k[2-9]{1,2}me3"))]})
)[, .(sumtot = sum(sumtot), summe = sum(summe)), .(chr, bin)]


CT <- mod_bin %>% mutate(beta_ChICTAPS = summe/sumtot) %>% mutate(chr = str_replace(chr, "chr", ""), location = paste0(chr, ":", bin))


x = merge.data.frame(CT, encode, by='location')%>% mutate(beta_encode = me/100)%>%
  ggplot()+geom_point(aes(x=beta_encode, y=beta_ChICTAPS), alpha=0.2, size=0.5)+theme_bw()



y= encode%>%  mutate(beta_encode = me/100)%>%  ggplot()+geom_histogram(aes(x= beta_encode), bins=200)+theme_bw()
                                                              
z = CT %>% ggplot()+geom_histogram(aes(y= beta_ChICTAPS), bins=200)+theme_bw()
                                                                    
(y/x)|(y/z)




library(tidyverse)
library(data.table)
library(Rcpp)
library(patchwork)

cppFunction({'List dist_1d(IntegerVector x, int MAX) {

  int n = x.size();
  IntegerVector from ((n * n - n) / 2);
  IntegerVector to   ((n * n - n) / 2);
  int k = 0;

  for(int i = 0; i < n; ++i) {
    for(int j = i + 1; j < n; ++j){
    
      from[k] = x[i];
      to[k]   = x[j];
      
      ++k;
    }
  }

  IntegerVector dist = to - from; 
  LogicalVector filter = dist < MAX;
  IntegerVector from_f = from[filter];
  IntegerVector dist_f = dist[filter];
  
  List final =  List::create(Named("from") = from_f, Named("dist") = dist_f);

 return final;
}'})

# per cell summary
meth_cuts_sum <- list.files("~/Desktop/post-doc/ChiC-TAPS/RPE-1_rawdata/rawdata/",  full.names = T) %>% #"pl1.tsv|pl9.tsv|pl16.tsv",
  map(function(x){ 
    fread(x)[V2 %in% c(1:22, "X") #& V5 %in% c("TA", "TT")
    ][, data.table(str_split_fixed(V1, ":", 6), chr = V2, bp = as.integer(V3), me = (V4 == "Z"))
    ][, .(N = uniqueN(paste(V3, V4, chr)), fme = sum(me) / .N ), .(V1)] 
  }) %>% rbindlist()

#meth_cuts_sum %>% write_tsv("Documents/Projects/chic_taps/data/RPE/all_cells_sum.tsv")
meth_cuts_sum <- fread("~/Desktop/post-doc/ChiC-TAPS/Fig3/all_cells_sum.tsv")

# cutoffs
filterlim <- data.table(modd = c("k27me3","k36me3", "k9me3"),
                        Nmin = c(4.2, 4, 4),
                        Nmax = c(5.5, 5, 5),
                        fmin = c(0.32, 0.7, 0.25),
                        fmax = c(0.41, 0.8, 0.35))


meth_cuts_sum[, mod := str_extract(V1, "k[2-9]{1,2}me[1-3]{1}")]
meth_cuts_sum[, pass := N > 10^filterlim[modd == mod]$Nmin & 
                N < 10^filterlim[modd == mod]$Nmax & 
                fme > filterlim[modd == mod]$fmin & 
                fme < filterlim[modd == mod]$fmax, .(mod)]
meth_cuts_sum[(pass), rank := frank(-N, ties.method = "first"), .(mod)]

# plotting cutoffs
meth_cuts_sum %>% ggplot() + 
  geom_point(aes(x = log10(N), y = fme, col = pass), size = 0.1) + facet_wrap(vars(mod)) 

# data
meth_cuts_raw <- list.files("~/Desktop/post-doc/ChiC-TAPS/RPE-1_rawdata/rawdata/",  "tsv.gz", full.names = T) %>% #"pl1.tsv|pl9.tsv|pl16.tsv",
  map(function(x){ 
    gc()
    fread(x)[V2 %in% c(1:22, "X") #& V5 %in% c("TA", "TT")
    ][, data.table(str_split_fixed(V1, ":", 6), chr = V2, bp = as.integer(V3), me = (V4 == "Z"))
    ][V1 %in% meth_cuts_sum[(pass) & rank < 251]$V1
    ][,.(V1, V3 = as.integer(V3), chr, bp, me)]
  })

# which bin in which mod
mod_bin <- rbindlist(map(meth_cuts_reads, function(x) {
  x[, .(N = as.double(.N)), .(V1, chr, bin = round(bp / 1e5))][, N := N / sum(N)]})
)[, mod := str_extract(V1, "k[2-9]{1,2}me3")
][, .(N = sum(N)), .(mod, chr, bin)]

mod_bin[chr == "1"][, N := N / quantile(N, 0.999), .(mod)][bin > 5 & bin < 200] %>%
  ggplot() + geom_col(aes(x = bin, y = N, fill = mod), position = position_identity(), alpha = 0.33333)


# pre processing for chic part
cuts_filt <- map(meth_cuts_raw, function(x){ 
  print(x$V1[1])
  x[,.(me = sum(me), tot = .N), .(V1, V3, chr)]
})

sourceCpp("~/Desktop/post-doc/ForkDOT/dist_1d.cpp")
cust_dist <- map(cuts_filt, function(x){
  print(x$V1[1])
  x[order(V3)][, dist_1d(V3, 2000), .(V1, chr)] })

# pre processing for methylation in relation to chic cuts
me_filt <- 
  map(meth_cuts_raw, function(x){ 
    print(x$V1[1])
    x[,.(N = as.double(.N)), .(V1, me, dist = abs(bp - as.numeric(V3)))]
  })



# figure
A <- rbindlist(map(1:12, function(x) cuts_filt[[x]][chr == "2" ])#& as.numeric(V3) > 17e7 & as.numeric(V3) < 26e7])
)[,.(N = as.double(.N)),.(V1, bp_round = round(as.numeric(V3) / 2e5) * 2e5)
][, N := N / sum(N), .(V1)
][, mod := str_extract(V1, "k[2-9]{1,2}me3")
][bp_round > 1e7 & bp_round < 5e7
][, N := N / quantile(N, 0.99), .(mod)
][N > 1, N := 1
][, cell := as.numeric(as.factor(V1)), .(mod)
] %>%
  ggplot() + geom_raster(aes(x = bp_round, y = cell, fill = N)) +
  scale_fill_viridis_c(option ="B", limits = c(0.2, 0.8), oob = scales::squish) +
  facet_grid(rows = vars(mod) ) + 
  coord_cartesian(expand = F)+
  theme_bw()+
  theme(panel.background = element_rect(fill = "black"), panel.grid = element_blank())

B1 <- rbindlist(cust_dist)[, mod := str_extract(V1, "k[2-9]{1,2}me3")
][, .(N = as.double(.N)), .(mod, V1, dist = round((dist + 0.0001) / 20) * 20)
  #][, tot := sum(N), .(V1)
][, cell := as.numeric(as.factor(V1)), .(mod)
][, N := N / sum(N), .(V1)
][dist > 0 & dist < 2000
][, N := N / max(N), .(mod)  
  #][N > 0.05, N := 0.05
] %>% dplyr::filter(mod == 'k27me3')%>%
  ggplot() +
  geom_raster(aes(x = dist, y = cell, fill = N)) +
  scale_fill_viridis_c(option = 'B', 
                       name= "k27me3", 
                       limits = c(0, 0.9), oob = scales::squish, direction=1
  )+
  facet_grid(rows = vars(mod), scales = "free") +
  coord_cartesian(expand = F)+
  theme_bw()

B2 <- rbindlist(cust_dist)[, mod := str_extract(V1, "k[2-9]{1,2}me3")
][, .(N = as.double(.N)), .(mod, V1, dist = round((dist + 0.0001) / 20) * 20)
  #][, tot := sum(N), .(V1)
][, cell := as.numeric(as.factor(V1)), .(mod)
][, N := N / sum(N), .(V1)
][dist > 0 & dist < 2000
][, N := N / max(N), .(mod)  
  #][N > 0.05, N := 0.05
] %>% dplyr::filter(mod == 'k36me3')%>%
  ggplot() +
  geom_raster(aes(x = dist, y = cell, fill = N)) +
  scale_fill_viridis_c(option = 'B', 
                       name= "k36me3", 
                       limits = c(0, 0.9), oob = scales::squish, direction=1
  )+
  facet_grid(rows = vars(mod), scales = "free") +
  coord_cartesian(expand = F)+
  theme_bw()

B3 <- rbindlist(cust_dist)[, mod := str_extract(V1, "k[2-9]{1,2}me3")
][, .(N = as.double(.N)), .(mod, V1, dist = round((dist + 0.0001) / 25) * 25)
  #][, tot := sum(N), .(V1)
][, cell := as.numeric(as.factor(V1)), .(mod)
][, N := N / sum(N), .(V1)
][dist > 0 & dist < 2000
][, N := N / max(N), .(mod)  
  #][N > 0.05, N := 0.05
] %>% dplyr::filter(mod == 'k9me3')%>%
  ggplot() +
  geom_raster(aes(x = dist, y = cell, fill = N)) +
  scale_fill_viridis_c(option = 'B', 
                       name= "k9me3", 
                       limits = c(0, 0.9), oob = scales::squish, direction=1
  )+
  facet_grid(rows = vars(mod), scales = "free") +
  coord_cartesian(expand = F)+
  theme_bw()

B1/B2/B3


D <- dcast(rbindlist(me_filt), V1 + dist ~ me, value.var = "N"
)[dist < 500
][is.na(`FALSE`), `FALSE` := 0
][is.na(`TRUE`), `TRUE` := 0
][, mod := str_extract(V1, "k[2-9]{1,2}me3")
][, tot := sum(`TRUE` + `FALSE`), .(mod, V1)
][, f := `TRUE` / (`TRUE` + `FALSE`)
][dist > 3 
][,.(mean = median(f), upper = quantile(f, 0.975), lower = quantile(f, 0.025)), .(mod, dist)
] %>%
  ggplot() +
  #geom_ribbon(aes(x = dist, ymin = lower, ymax = upper, fill = mod), alpha = 0.5) +
  geom_line(aes(x = dist, y = mean, col = mod)) +
  #geom_smooth(aes(x = dist, y = mean, col = mod), span =0.1, alpha=0.5, se=F, linetype =2)+
  #facet_grid(rows = vars(mod), scales = "free")+
  theme_bw()+  scale_color_manual(values=c("#D674BA","#86BF88","#7A416A"))+theme_bw()+
  scale_fill_manual(values=c("#D674BA","#86BF88","#7A416A"))+theme_bw()+ylab('Fraction CpG Methylation')+xlab('Distance from ChIC cut (bp)')

C <- D + geom_ribbon(aes(x = dist, ymin = lower, ymax = upper, fill = mod), alpha = 0.5) 




(A | B) / (C | D)



#RPE_WGBS <- readRDS('~/Desktop/post-doc/ChiC-TAPS/ExtendedDataFigures/Supp2/rpe-jvl.rds')%>%
#  dplyr::filter(chr %in% c((paste0('chr', 1:22)), "chrX"), cov>0)%>% 
 # select(chr, start, cov, beta)%>%
#  group_by(chr)%>%
#  mutate(roundbp = round(start/1e5)*1e5)%>% group_by(chr, roundbp) %>% summarize(me = mean(beta))%>% mutate(source = "WGBS")

#RPE_CT <- rbindlist(map(1:3, function(x) cuts_filt[[x]][])#& as.numeric(V3) > 17e7 & as.numeric(V3) < 26e7])
#)[, mod := str_extract(V1, "K[2-9]{1,2}m3")
#][, cell := as.numeric(as.factor(V1)), .(mod)
#]%>% group_by(chr)%>%mutate(roundbp = round(V3/1e5)*1e5)%>% group_by(chr, roundbp)%>%
#  summarize(me = sum(me)/sum(tot)*100)%>% mutate(source = "TAPS", chr = paste0('chr', chr))


encode <- readRDS('~/Desktop/post-doc/ChiC-TAPS/ExtendedDataFigures/Supp2/rpe-jvl.rds')%>%
  #dplyr::filter(cov>0)%>% 
  select(chr, start, cov, beta)%>% dplyr::filter(cov>1)%>%
  group_by(chr)%>%
  mutate(roundbp = round(start/2.5e4)*2.5e4)%>% group_by(chr, roundbp) %>% summarize(me = mean(beta))
encode <- encode%>% ungroup()%>% mutate(chr = str_replace(chr, "chr", ""), location = paste0(chr, ":", roundbp))


mod_bin <- rbindlist(map(meth_cuts_reads, function(x) {
  
  x[V1 %in% filt_cells$V1
  ][, .(sumtot = sum(tot), summe=sum(me)), .(chr, bin = (round(bp / 2.5e4)*2.5e4),  mod = str_extract(V1, "k[2-9]{1,2}me3"))]})
)[, .(sumtot = sum(sumtot), summe = sum(summe)), .(chr, bin)]


CT <- mod_bin %>% mutate(beta_ChICTAPS = summe/sumtot) %>% mutate(chr = str_replace(chr, "chr", ""), location = paste0(chr, ":", bin))
merge.data.frame(CT, encode, by='location')%>% mutate(beta_encode = me/100)%>%
  ggplot()+geom_bin2d(aes(x=beta_encode, y=beta_ChICTAPS), bins=200)+
  geom_smooth(aes(x=beta_encode, y=beta_ChICTAPS), method = 'lm', linetype =2, col ="green1")+coord_fixed()+
  theme_bw()+ylab('TET-assisted Pyrdine Borane sequencing [TAPS] ')+
  scale_fill_viridis_c(option = "C", trans='log10')+
  xlab('Whole Genome Bisulfite Sequencing [WGBS]')+ggtitle('average CpG methylation [25kb-1]')

tidy 
corrme <-merge.data.frame(CT, encode, by='location')%>% mutate(beta_encode = me/100)

cor.test(corrme$beta_encode, corrme$beta_ChICTAPS, method = 'pearson')

gc()

