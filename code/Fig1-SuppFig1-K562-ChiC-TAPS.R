

library(tidyverse)
library(data.table)
library(Rcpp)
library(patchwork)


gc()


# per cell summary
meth_cuts_sum <- list.files("~/Desktop/post-doc/ChiC-TAPS/K562-rawdata/",  "tsv.gz", full.names = T) %>% #"pl1.tsv|pl9.tsv|pl16.tsv",
  map(function(x){ 
    fread(x)[V2 %in% c(1:22, "X") #& V5 %in% c("TA", "TT")
    ][, data.table(str_split_fixed(V1, ":", 6), chr = V2, bp = as.integer(V3), me = (V4 == "Z"))
    ][, .(N = uniqueN(paste(V3, V4, chr)), fme = sum(me) / .N ), .(V1)] 
  }) %>% rbindlist()


meth_cuts_sum %>% write_tsv("~/Desktop/post-doc/ChiC-TAPS/K562-rawdata/all_cells_sum.tsv")
meth_cuts_sum <- fread("~/Desktop/post-doc/ChiC-TAPS/K562-rawdata/all_cells_sum.tsv")
#meth_cuts_sum %>% group_by(mod)%>% summarize(N= n())
# cutoffs
filterlim <- data.table(modd = c("K27m3","K36m3", "K9m3"),
                        Nmin = c(4.2, 3.8, 4),
                        Nmax = c(5.5, 5, 5),
                        fmin = c(0.05, 0.45, 0.05),
                        fmax = c(0.4, 0.8, 0.4))


meth_cuts_sum[, mod := str_extract(V1, "K[2-9]{1,2}m[1-3]{1}")]
meth_cuts_sum[, pass := N > 10^filterlim[modd == mod]$Nmin & 
                N < 10^filterlim[modd == mod]$Nmax & 
                fme > filterlim[modd == mod]$fmin & 
                fme < filterlim[modd == mod]$fmax, .(mod)]
meth_cuts_sum[(pass), rank := frank(-N, ties.method = "first"), .(mod)]

# plotting cutoffs
meth_cuts_sum %>% ggplot() + 
  geom_point(aes(x = log10(N), y = fme, col = mod
                 ), size = 0.5) + facet_wrap(vars(mod))+theme_bw()+
  geom_vline(xintercept = c(3.8,4,4.2), linetype=2)+
  geom_hline(yintercept = c(0.05, 0.45, 0.05), linetype=2)+
  scale_color_manual(values=c("#D674BA","#86BF88","#7A416A"))+theme_bw()+theme(legend.position = 'none')+
  ylab('CpG methylation per cell')+xlab('unique ChIC reads per cell')+ggtitle('K562')
  
gc()
# data
meth_cuts_raw <- list.files("~/Desktop/post-doc/ChiC-TAPS/K562-rawdata/",  "tsv.gz", full.names = T) %>% #"pl1.tsv|pl9.tsv|pl16.tsv",
  map(function(x){ 
    gc()
    fread(x)[V2 %in% c(1:22, "X") #& V5 %in% c("TA", "TT")
    ][, data.table(str_split_fixed(V1, ":", 6), chr = V2, bp = as.integer(V3), me = (V4 == "Z"))
    ][V1 %in% meth_cuts_sum[(pass) & rank < 251]$V1
    ][,.(V1, V3 = as.integer(V3), chr, bp, me)]
  })

# which bin in which mod
#mod_bin <- rbindlist(map(meth_cuts_reads, function(x) {
 # x[, .(N = as.double(.N)), .(V1, chr, bin = round(bp / 1e5))][, N := N / sum(N)]})
##)[, mod := str_extract(V1, "k[2-9]{1,2}me3")
#][, .(N = sum(N)), .(mod, chr, bin)]

#mod_bin[chr == "1"][, N := N / quantile(N, 0.999), .(mod)][bin > 5 & bin < 200] %>%
#  ggplot() + geom_col(aes(x = bin, y = N, fill = mod), position = position_identity(), alpha = 0.33333)
#readRDS('~/Desktop/post-doc/ChiC-TAPS/Fig2/meth_cuts_reads.RDS')
sourceCpp("~/Desktop/post-doc/ForkDOT/dist_1d.cpp")
# pre processing for chic part
cuts_filt <- map(meth_cuts_raw, function(x){ 
  print(x$V1[1])
  x[,.(me = sum(me), tot = .N), .(V1, V3, chr)]
})

cust_dist <- map(cuts_filt, function(x){
  print(x$V1[1])
  x[order(V3)][, dist_1d(V3, 2000), .(V1, chr)] })

# pre processing for methylation in relation to chic cuts
me_filt <- 
  map(meth_cuts_raw, function(x){ 
    print(x$V1[1])
    x[,.(N = as.double(.N)), .(V1, me, dist = abs(bp - as.numeric(V3)))]
  })


fme1 <- rbindlist(map(1:3, function(x) cuts_filt[[x]][chr == "2" ])#& as.numeric(V3) > 17e7 & as.numeric(V3) < 26e7])
)[, mod := str_extract(V1, "K[2-9]{1,2}m3")
][V3 > 1e7 & V3 < 5e7
][, cell := as.numeric(as.factor(V1)), .(mod)
]%>% mutate(roundbp = round(V3/1e5)*1e5)%>% group_by(mod, roundbp, cell)%>%
  summarize(fme = sum(me)/sum(tot)) %>%dplyr::filter(cell<201)%>%
  ggplot() + geom_raster(aes(x = roundbp, y = cell, fill = fme)) +
  scale_fill_distiller(palette = 'RdYlBu', limits = c(0.1, 0.7), oob = scales::squish) +
  facet_grid(rows = vars(mod), space= 'free') + 
  coord_cartesian(expand = F)+
  theme_bw()+
  theme(panel.background = element_rect(fill = "lightgrey"), panel.grid = element_blank())+theme(legend.position = "bottom")


rbindlist(map(1:3, function(x) cuts_filt[[x]][chr == "2" ])#& as.numeric(V3) > 17e7 & as.numeric(V3) < 26e7])
)[, mod := str_extract(V1, "K[2-9]{1,2}m3")
][V3 > 1e7 & V3 < 5e7
][, cell := as.numeric(as.factor(V1)), .(mod)
]%>% mutate(roundbp = round(V3/1e5)*1e5)%>% group_by(mod, roundbp)%>%
  summarize(fme = sum(me)/sum(tot), sumtot = n())%>%
  ggplot+geom_histogram(aes(x=sumtot))+facet_grid(rows=vars(mod))+scale_y_log10()


fme0 <- rbindlist(map(1:3, function(x) cuts_filt[[x]][chr == "2" ])#& as.numeric(V3) > 17e7 & as.numeric(V3) < 26e7])
)[, mod := str_extract(V1, "K[2-9]{1,2}m3")
][V3 > 1e7 & V3 < 5e7
][, cell := as.numeric(as.factor(V1)), .(mod)
]%>% mutate(roundbp = round(V3/1e5)*1e5)%>% group_by(mod, roundbp)%>%
  summarize(fme = sum(me)/sum(tot), sumtot = sum(tot)) %>%# dplyr::filter(sumtot>1000)%>%
  ggplot() + geom_line(aes(x = roundbp, y = fme) )+
  scale_color_manual(values=c("black"))+theme_bw()+theme(legend.position = 'none')+ ylim(0,1)+
  #scale_fill_distiller(palette = 'RdYlBu', limits = c(0.0, 0.3), oob = scales::squish) +
 # facet_grid(rows = vars(mod), scales = 'free') +# scale_y_log10()+# ylim(0,100)+
  coord_cartesian(expand = F)#+
  #theme_bw()+
  #theme(panel.background = element_rect(fill = "#4575b4"), panel.grid = element_blank())

(fme0/fme1)+plot_layout(heights = c(1,5))


# figure
rbindlist(map(1:3, function(x) cuts_filt[[x]][chr == "2" ])#& as.numeric(V3) > 17e7 & as.numeric(V3) < 26e7])
)[,.(N = as.double(.N)),.(V1, bp_round = round(as.numeric(V3) / 1e5) * 1e5)
][, N := N / sum(N), .(V1)
][, mod := str_extract(V1, "K[2-9]{1,2}m3")
][bp_round > 1e7 & bp_round < 5e7
][, N := N / quantile(N, 0.99), .(mod)
][N > 1, N := 1
][, cell := as.numeric(as.factor(V1)), .(mod)
][] %>% dplyr::filter(cell<201)%>%
  ggplot() + geom_raster(aes(x = bp_round, y = cell, fill = N)) +
  scale_fill_viridis_c(option ="B", limits = c(0.2, 0.8), oob = scales::squish) +
  facet_grid(rows = vars(mod), space= 'free') + 
  coord_cartesian(expand = F)+
  theme_bw()+
  theme(panel.background = element_rect(fill = "black"), panel.grid = element_blank())


rbindlist(map(1:3, function(x) cuts_filt[[x]][chr == "2" ])#& as.numeric(V3) > 17e7 & as.numeric(V3) < 26e7])
)[,.(N = as.double(.N)),.(V1, bp_round = round(as.numeric(V3) / 1e5) * 1e5)
][, N := N / sum(N), .(V1)
][, mod := str_extract(V1, "K[2-9]{1,2}m3")
][bp_round > 1e7 & bp_round < 5e7
][, N := N / quantile(N, 0.99), .(mod)
][N > 1, N := 1
][, cell := as.numeric(as.factor(V1)), .(mod)
][] %>% dplyr::filter(cell<201)%>% group_by(mod, bp_round)%>% summarize(N = mean(N))%>%
  ggplot() + geom_line(aes(x = bp_round, y = N, color = mod), size=1) +
  #scale_fill_viridis_c(option ="B", limits = c(0.2, 0.8), oob = scales::squish) +
  facet_grid(rows = vars(mod), space= 'free') + 
  coord_cartesian(expand = F)+
  theme_bw()+scale_color_manual(values=c("#D674BA","#86BF88","#7A416A"))

library(colorspace)

B1 <- rbindlist(cust_dist)[, mod := str_extract(V1, "K[2-9]{1,2}m3")
][, .(N = as.double(.N)), .(mod, V1, dist = round((dist + 0.0001) / 20) * 20)
  #][, tot := sum(N), .(V1)
][, cell := as.numeric(as.factor(V1)), .(mod)
][, N := N / sum(N), .(V1)
][dist > 50 & dist < 1500
][, N := N / max(N), .(mod)  
  #][N > 0.05, N := 0.05
] %>% dplyr::filter(mod == 'K27m3' & dist>39)%>%group_by(cell)%>% mutate(N=scale(N))%>%
  ggplot() +
  geom_raster(aes(x = dist, y = cell, fill = N)) +
  scale_fill_viridis_c(option = 'B', 
                       name= "ChIC reads Pair Correlation density", 
                       limits = c(-1.5, 3), oob = scales::squish, direction=1
  )+
  facet_grid(rows = vars(mod), scales = "free") +
  coord_cartesian(expand = F)+
  theme_bw()+theme(legend.position = 'none')


B2 <- rbindlist(cust_dist)[, mod := str_extract(V1, "K[2-9]{1,2}m3")
][, .(N = as.double(.N)), .(mod, V1, dist = round((dist + 0.0001) / 20) * 20)
  #][, tot := sum(N), .(V1)
][, cell := as.numeric(as.factor(V1)), .(mod)
][, N := N / sum(N), .(V1)
][dist > 50 & dist < 1500
][, N := N / max(N), .(mod)  
  #][N > 0.05, N := 0.05
] %>% dplyr::filter(mod == 'K36m3'  & dist>39)%>% group_by(cell)%>% mutate(N=scale(N))%>%
  ggplot() +
  geom_raster(aes(x = dist, y = cell, fill = N)) +
  scale_fill_viridis_c(option = 'B', 
                       name= "ChIC reads Pair Correlation density", 
                      limits = c(-1.5, 3), oob = scales::squish, direction =1
                       )+
  facet_grid(rows = vars(mod), scales = "free") +
  coord_cartesian(expand = F)+
  theme_bw()+theme(legend.position = 'none')

B3 <- rbindlist(cust_dist)[, mod := str_extract(V1, "K[2-9]{1,2}m3")
][, .(N = as.double(.N)), .(mod, V1, dist = round((dist + 0.0001) / 20) * 20)
  #][, tot := sum(N), .(V1)
][, cell := as.numeric(as.factor(V1)), .(mod)
][, N := N / sum(N), .(V1)
][dist > 50 & dist < 1500
][, N := N / max(N), .(mod)  
  #][N > 0.05, N := 0.05
] %>% dplyr::filter(mod == 'K9m3')%>%group_by(cell)%>% mutate(N=scale(N))%>%
  ggplot() +
  geom_raster(aes(x = dist, y = cell, fill = N)) +
  scale_fill_viridis_c(option = 'B', 
                       name= "ChIC reads Pair Correlation density", 
                       limits = c(-1.5, 3), oob = scales::squish,
                       direction = 1
  )+
  facet_grid(rows = vars(mod), scales = "free") +
  coord_cartesian(expand = F)+
  theme_bw()+theme(legend.position = 'bottom')


B1/B2/B3


D <- dcast(rbindlist(me_filt), V1 + dist ~ me, value.var = "N"
)[dist < 500
][is.na(`FALSE`), `FALSE` := 0
][is.na(`TRUE`), `TRUE` := 0
][, mod := str_extract(V1, "K[2-9]{1,2}m3")
][, tot := sum(`TRUE` + `FALSE`), .(mod, V1)
][, f := `TRUE` / (`TRUE` + `FALSE`)
][dist > 3 
][,.(mean = median(f), upper = quantile(f, 0.975), lower = quantile(f, 0.025)), .(mod, dist)
][] %>%
  ggplot() +
  #geom_ribbon(aes(x = dist, ymin = lower, ymax = upper, fill = mod), alpha = 0.5) +
  geom_line(aes(x = dist, y = mean, col = mod)) +
  facet_grid(rows = vars(mod), scales = "free")

C <- D + geom_ribbon(aes(x = dist, ymin = lower, ymax = upper, fill = mod), alpha = 0.5) 

meth_cuts_raw


(A | B) / (C | D)

###extended data figure 1####


K562WGBS <- readRDS('~/Desktop/post-doc/ChiC-TAPS/ExtendedDataFigures/k562-wgbs/k562-encode-wgbs-1-ENCFF721JMB.rds')%>%
  dplyr::filter(chr %in% c((paste0('chr', 1:22)), "chrX"), cov>0)%>% 
  select(chr, start, cov, beta)%>%
  group_by(chr)%>%
  mutate(roundbp = round(start/1e4)*1e4)%>% group_by(chr, roundbp) %>% summarize(me = mean(beta))%>% mutate(source = "WGBS")

K562_CT <- rbindlist(map(1:3, function(x) cuts_filt[[x]][])#& as.numeric(V3) > 17e7 & as.numeric(V3) < 26e7])
)[, mod := str_extract(V1, "K[2-9]{1,2}m3")
][, cell := as.numeric(as.factor(V1)), .(mod)
]%>% group_by(chr)%>%mutate(roundbp = round(V3/1e4)*1e4)%>% group_by(chr, roundbp)%>%
  summarize(me = sum(me)/sum(tot)*100)%>% mutate(source = "TAPS", chr = paste0('chr', chr))


bind_rows(K562_CT, K562WGBS)%>% mutate(loc = paste0(chr, "-", roundbp))%>% ungroup()%>% select(loc, source, me)%>% 
  pivot_wider(names_from = source, values_from = me )%>%
  ggplot()+geom_bin2d(aes(y=TAPS, x=WGBS), bins=100)+
  geom_smooth(aes(y=TAPS, x=WGBS), method = 'lm', linetype =2, col ="green1")+coord_fixed()+
  theme_bw()+ylab('TET-assisted Pyrdine Borane sequencing [TAPS] ')+
  scale_fill_viridis_c(option = "C", trans='log10')+
  xlab('Whole Genome Bisulfite Sequencing [WGBS]')+ggtitle('average CpG methylation [10kb-1]')

corrme <- bind_rows(K562_CT, K562WGBS)%>% mutate(loc = paste0(chr, "-", roundbp))%>% ungroup()%>% select(loc, source, me)%>% 
  pivot_wider(names_from = source, values_from = me )

cor.test(corrme$TAPS, corrme$WGBS, method = 'pearson')

gc()

#encode <- NULL

K562_ChIP<-list((fread('~/Desktop/post-doc/ChiC-TAPS/ExtendedDataFigures/k562-chip/h3k27me3/H3K27me3-ENCFF190OWE.tsv.gz')%>%
                mutate(mod = 'K27m3', method = 'chip')),
             (fread('~/Desktop/post-doc/ChiC-TAPS/ExtendedDataFigures/k562-chip/h3k27me3/H3K27me3-ENCFF692KQZ.tsv.gz')%>%
                mutate(mod = 'K27m3', method = 'chip')),
             
             (fread('~/Desktop/post-doc/ChiC-TAPS/ExtendedDataFigures/k562-chip/h3k36me3/H3K36me3-ENCFF639PLN.tsv.gz')%>%
                mutate(mod = 'K36m3', method = 'chip')),
             (fread('~/Desktop/post-doc/ChiC-TAPS/ExtendedDataFigures/k562-chip/h3k36me3/H3K36me3-ENCFF673KBG.tsv.gz')%>%
                mutate(mod = 'K36m3', method = 'chip')),
             
             (fread('~/Desktop/post-doc/ChiC-TAPS/ExtendedDataFigures/k562-chip/h3k9me3/H3K9me3-ENCFF146NLP.tsv.gz')%>%
                mutate(mod = 'K9m3', method = 'chip')),
             (fread('~/Desktop/post-doc/ChiC-TAPS/ExtendedDataFigures/k562-chip/h3k9me3/H3K9me3-ENCFF559DHZ.tsv.gz')%>%
                mutate(mod = 'K9m3', method = 'chip'))
)%>% rbindlist()%>% setNames(c('chr', 'bp', 'mod', 'method'))%>% group_by(chr, mod)%>%
  mutate(roundbp = round(bp/1e5)*1e5)%>% group_by(chr, mod, roundbp) %>% summarize(count = n())%>% mutate(source = 'CHIP')
  

K562_ChIC <- rbindlist(map(1:3, function(x) cuts_filt[[x]][])#& as.numeric(V3) > 17e7 & as.numeric(V3) < 26e7])
)[, mod := str_extract(V1, "K[2-9]{1,2}m3")
]%>% select(mod, chr, V3) %>% group_by(chr, mod)%>%
  mutate(roundbp = round(V3/1e5)*1e5,  chr = paste0('chr', chr))%>% group_by(chr, mod, roundbp) %>% summarize(count = n())%>% mutate(source = 'CHIC')

merged <- bind_rows(K562_ChIC, K562_ChIP)%>% mutate(loc = paste0(chr, "-", roundbp))%>% ungroup()%>% group_by(mod)%>% group_split()


k27 <- merged[[1]]
k36 <- merged[[2]]
#k9 <- merged[[3]]

comp1 <- k27 %>%  ungroup()%>% select(loc, source, count)%>% 
  pivot_wider(names_from = source, values_from = count )%>%
  ggplot()+geom_bin2d(aes(y=CHIC, x=CHIP), bins=100)+
  geom_smooth(aes(y=CHIC, x=CHIP), method = 'lm', linetype =2, col ="green")+coord_equal()+
  theme_bw()+ylab('ChIC')+ ylim(0,5e3)+ xlim(0,5e3)+
  scale_fill_viridis_c(option = "C", trans='log10')+
  xlab('CHIP')+ggtitle('H3K27me3 signal \n per 100kb')
k27.test<-k27 %>%  ungroup()%>% select(loc, source, count)%>% 
  pivot_wider(names_from = source, values_from = count )
 
cor.test(k27.test$CHIC, k27.test$CHIP, method = 'pearson')


comp2 <-  k36 %>%  ungroup()%>% select(loc, source, count)%>% 
    pivot_wider(names_from = source, values_from = count )%>%
    ggplot()+geom_bin2d(aes(y=CHIC, x=CHIP), bins=100)+
    geom_smooth(aes(y=CHIC, x=CHIP), method = 'lm', linetype =2, col ="green")+#coord_equal()+
    theme_bw()+ylab('ChIC')+ ylim(0,5e3)+ xlim(0,5e3)+
  scale_fill_viridis_c(option = "C", trans='log10')+
  xlab('CHIP')+ggtitle('H3K36me3 signal \n per 100kb')

k36.test<-k36 %>%  ungroup()%>% select(loc, source, count)%>% 
  pivot_wider(names_from = source, values_from = count)

cor.test(k36.test$CHIC, k36.test$CHIP, method = 'pearson')

comp3 <-  k9 %>%  ungroup()%>% select(loc, source, count)%>% 
    pivot_wider(names_from = source, values_from = count )%>%
    ggplot()+geom_point(aes(y=CHIC, x=CHIP), alpha=0.15, size=0.2)+
    geom_smooth(aes(y=CHIC, x=CHIP), method = 'lm', linetype =2, col ="red")+coord_fixed()+
    theme_bw()+ylab('ChIC')+ ylim(0,5e3)+ xlim(0,5e3)+
  xlab('CHIP')+ggtitle('H3K9me3 signal per 100kb')
  
comp1|comp2

cor.test(k27$count, corrme$WGBS, method = 'pearson')


k9me3_geisenberger <- rbindlist(map(1:3, function(x) cuts_filt[[x]][])#& as.numeric(V3) > 17e7 & as.numeric(V3) < 26e7])
)[, mod := str_extract(V1, "K[2-9]{1,2}m3")
]%>% select(mod, chr, V3) %>% group_by(chr, mod)%>% dplyr::filter(mod == 'K9m3')%>%
  mutate(roundbp = round(V3/5e4)*5e4)%>% group_by(chr, mod, roundbp) %>% summarize(count = n())%>% mutate(source = 'CHIC_Geisenberger',  loc = paste0(chr, "-", roundbp))%>%
 ungroup()%>%  select(loc, count , source)

k9me3_Zeller <- fread('~/Desktop/post-doc/ChiC-TAPS/ExtendedDataFigures/k562-chip/GSE164779_RAW/GSM5018608_K562-EtOH-H3K9me3.G1sorted.merged.sorted.tagged.countTable.csv.gz')%>%
  mutate(loc = paste0(sampleName, "-", V2))%>% select(-c(sampleName,V2, V3))%>% pivot_longer(!loc, names_to = 'cell' ) %>% dplyr::filter(loc != 'reference_name-start')%>%
  group_by(loc)%>% drop_na()%>%summarize(count= sum(value))%>% mutate(source = 'CHIC_Zeller')


bind_rows(k9me3_geisenberger, k9me3_Zeller)%>%
  ungroup()%>% select(loc, source, count)%>% 
  pivot_wider(names_from = source, values_from = count )%>%
  ggplot()+geom_bin2d(aes(y=CHIC_Geisenberger, x=CHIC_Zeller), bins =100)+
  geom_smooth(aes(y=CHIC_Geisenberger, x=CHIC_Zeller), method = 'lm', linetype =2, col ="green")+coord_fixed()+
  scale_fill_viridis_c(option = "C", trans='log10')+
  theme_bw()+ylab('ChIC (This study)')+ ylim(0,3e3)+ xlim(0,3e3)+
  xlab('ChIC (Zeller et al.,)')+ggtitle('H3K9me3 signal per 50kb')


k9.test<-bind_rows(k9me3_geisenberger, k9me3_Zeller)%>%
  ungroup()%>% select(loc, source, count)%>% 
  pivot_wider(names_from = source, values_from = count )


cor.test(k9.test$CHIC_Geisenberger, k9.test$CHIC_Zeller, method = 'pearson')


