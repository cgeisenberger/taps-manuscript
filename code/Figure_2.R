library(tidyverse)
library(data.table)
library(patchwork)


setwd('~/Desktop/post-doc/ChiC-TAPS/')

loess_wrap <- function(x, y, sspan){
  
  X <- data.table(x2 = c(x - max(x), x, x + max(x)), 
                  y2 = rep(y, 3))
  
  predict(loess(y2 ~ x2, data = X, span = sspan), x)
}


roll_median <- function(x, windowsize){
  
  
  frollapply(c(rep(NA, windowsize), x, rep(NA, windowsize)), 
             n = windowsize * 2 + 1, 
             align = "center",
             FUN = function(y){
               xx <- y[!is.na(y)]
               
               approx((1:length(xx)) / length(xx), sort(xx), 0.5)$y}
  )[(1:length(x)) + windowsize]
}


rollfun            <- function(x, w, FUN = "se", window = 0.25){
  
  was <- ceiling(length(x) * window / 2) 
  
  rowe  <- frollsum(c(rep(NA, was), w, rep(NA, was)),  
                    n = was*2+1, 
                    align = "center",
                    na.rm = T, hasNA = T)[1:length(x) + was]
  
  roweme  <- frollsum(c(rep(NA, was), x*w, rep(NA, was)),  
                      n = was*2+1, 
                      align = "center",
                      na.rm = T, hasNA = T)[1:length(x) + was] / rowe
  
  if(FUN == "mean"){return(roweme)}
  
  if(FUN == "se"){
    return(
      sqrt(frollsum(c(rep(NA, was), (x - roweme)^2 * w, rep(NA, was)),  
                    n = was*2+1, 
                    align = "center",
                    na.rm = T)[1:length(x) + was] / rowe))}
}






# per cell summary

meth_cuts_sum <- fread("~/Desktop/post-doc/ChiC-TAPS/Fig2/all_cells_sum.tsv")

facs <- fread("~/Desktop/post-doc/ChiC-TAPS/Fig2/facs.csv"
)[, mod := str_extract(V1, "k[2-9]{1,2}me3")
][, rank := frank(progression + runif(.N, max = 0.00001), ties.method = "dense"), .(mod)]

sphase <- readRDS("~/Desktop/post-doc/ChiC-TAPS/Fig2/JvB-exp049-15m-EdU-beads_sphase.RDS")$final_order %>% distinct(cell, rank)

RT <- fread("~/Desktop/post-doc/ChiC-TAPS/Fig2/JvB-exp049-15m-EdU-beads_single_fork.tsv.gz"
)[, rank := setNames(sphase$rank, sphase$cell)[as.character(cell)]
][state == 2
][, .(tmed =  median(rank, na.rm = T)), .(bin = round(center / 5e5), chr)
][, rank_bin := frank(round(tmed, 2), ties.method = "dense")]

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
{meth_cuts_sum %>% ggplot() + 
    geom_point(aes(x = log10(N), y = fme, col = pass), size = 0.1) + facet_wrap(vars(mod)) 
}

# per read collapsing

# only once 

meth_cuts_reads <- readRDS("~/Desktop/post-doc/ChiC-TAPS/RPE-1_rawdata/meth_cuts_reads.RDS")#[[1]][, max(bp) - min(bp), .(V1, chr, bin)]

# which bin in which mod
mod_bin <- rbindlist(map(meth_cuts_reads, function(x) {
  x[, .(N = as.double(.N)), .(V1, chr, bin = round(bp / 2e4))][, N := N / sum(N)]})
)[, mod := str_extract(V1, "k[2-9]{1,2}me3")
][, .(N = sum(N)), .(mod, chr, bin)]


me_avg <- merge.data.table(
  rbindlist(map(meth_cuts_reads, function(x){
    
    x[, .(me = sum(me), tot = sum(tot)), .(V1, chr, bin = round(bp / 2e4))
    ][mod_bin[,.SD[N == max(N)], .(chr, bin)], on = c("chr" , "bin")
    ][str_detect(V1, mod)
    ][, .(me = as.double(sum(me)), tot = as.double(sum(tot))), .(mod, V1, chr, bin = round(bin / 25))]
    
  })),
  RT, by = c("chr", "bin")
)[, rank := setNames(facs$rank, facs$V1)[V1]
][][, c("me", "tot") :=  .(me / sum(tot), tot / sum(tot)), .(V1)
][,.(me = mean(me / tot), tot = sum(tot)), .(mod, rank_bin, rank)
][, rank2 := frank(rank, ties.method = "dense"), .(mod)
][, mez := scale(me)[,1],.(mod, rank_bin)
][order(rank2)]

me_avg[abs(mez) > 2, mez := 2 * sign(mez)]

#me_avg[, me_lo := predict(loess(mez ~ rank_bin, data = .SD, span = 0.1), .SD$rank_bin), .(mod, rank2)]
#me_avg[, me_lo := loess_wrap(rank_bin, mez),  .(mod, rank2)]
me_avg[rank_bin < 86, me_lo2 := loess_wrap(rank2, mez, sspan = 0.2),  .(mod, rank_bin)]

#me_avg[, me_lo2 := NULL]

genome_example <- merge.data.table(
  merge.data.table(
    rbindlist(map(meth_cuts_reads, function(x){
      
      x[chr == "1", .(me = sum(me), tot = sum(tot)), .(V1, chr, bin = round(bp / 2e4))
      ][mod_bin[,.SD[N == max(N)], .(chr, bin)], on = c("chr" , "bin")
      ][str_detect(V1, mod)
      ][, merge_bin := round(bin / 25)]
      
    })),
    RT, by.x = c("chr", "merge_bin"), by.y = c("chr", "bin")
  )[chr == "1"
  ][, rank := setNames(facs$rank, facs$V1)[V1]],
  me_avg, by = c("mod", "rank", "rank_bin"), all.x = T
)

genome_example[, ccrank := rank + runif(1, 0, 0.001), .(V1)]
genome_example[, ccrank2 := frank(ccrank, ties.method = "dense"), .(mod)]

#genome_example[, mez_bin := scale(me.x / tot.x)[,1], .(mod, bin, chr)]
#genome_example[abs(mez_bin) > 2, mez_bin := 2 * sign(mez_bin) ]
#genome_example[, N := .N, .(mod, bin, chr)]
#genome_example[N > 10, me_lo_bin := predict(loess(mez_bin ~ ccrank, data = .SD), .SD$ccrank), .(mod, bin, chr)]

########
#######loading_BG_FACSdata####
rep1 <- read.FCS(filename = "~/Desktop/post-doc/ForkDOT/Fig1-cowplot/RPE_FUCCI-6x.fcs", transformation="linearize")# %>% as_data_frame(.@exprs)
rep2 <- read.FCS(filename = "~/Desktop/post-doc/ForkDOT/Fig1-cowplot/RPE_FUCCI-8x.fcs", transformation="linearize")# %>% as_data_frame(.@exprs)
rep3 <- read.FCS(filename = "~/Desktop/post-doc/ForkDOT/Fig1-cowplot/RPE_FUCCI-10x.fcs", transformation="linearize")# %>% as_data_frame(.@exprs)



rep1 <- flowCore::read.FCS(filename = "~/Desktop/post-doc/ForkDOT/Fig1-cowplot/RPE_FUCCI-6x.fcs", transformation="linearize")# %>% as_data_frame(.@exprs)
rep2 <- flowCore::read.FCS(filename = "~/Desktop/post-doc/ForkDOT/Fig1-cowplot/RPE_FUCCI-8x.fcs", transformation="linearize")# %>% as_data_frame(.@exprs)
rep3 <- flowCore::read.FCS(filename = "~/Desktop/post-doc/ForkDOT/Fig1-cowplot/RPE_FUCCI-10x.fcs", transformation="linearize")# %>% as_data_frame(.@exprs)
rep1 <- as.data.frame(rep1@exprs)%>%  mutate(log_green = log10(`FL1-H`),
                                             log_red =  log10(`FL11-H`),
                                             blue10e4 =  (`FL8-H`/10e4),
                                             Z_green = scale(log_green),
                                             Z_red= scale(log_red),
                                             Z_blue = scale(blue10e4)) %>%dplyr::filter(
                                               abs(Z_green)<2.2,
                                               abs(Z_red)<2.2,
                                               abs(Z_blue)<2.2)%>%  drop_na()

rep2 <- as.data.frame(rep2@exprs) %>% mutate(log_green = log10(`FL1-H`),
                                             log_red =  log10(`FL11-H`),
                                             blue10e4 =  (`FL8-H`/10e4),
                                             Z_green = scale(log_green),
                                             Z_red= scale(log_red),
                                             Z_blue = scale(blue10e4)) %>%dplyr::filter(
                                               abs(Z_green)<2.2,
                                               abs(Z_red)<2.2,
                                               abs(Z_blue)<2.2)%>%  drop_na()

rep3 <- as.data.frame(rep3@exprs) %>% mutate(log_green = log10(`FL1-H`),
                                             log_red =  log10(`FL11-H`),
                                             blue10e4 =  (`FL8-H`/10e4),
                                             Z_green = scale(log_green),
                                             Z_red= scale(log_red),
                                             Z_blue = scale(blue10e4)) %>%dplyr::filter(
                                               abs(Z_green)<2.2,
                                               abs(Z_red)<2.2,
                                               abs(Z_blue)<2.2)%>%  drop_na()


dta <-bind_rows(rep1, rep2, rep3) %>% select(Z_green, Z_red, Z_blue) %>% mutate(Z_green = Z_green-0.3, 
                                                                                Z_red=Z_red+0.2)
#dta%>%
#  ggplot()+
#  geom_density2d( aes(y=Z_red, x= Z_green),alpha =0.5, size=0.6, color="black", binwidth = 0.03 )
library(infer)
dta_small=dta%>% slice_sample(n=1e4) %>% as.data.frame()%>% select(starts_with("Z"))

x= dta_small%>%
  ggplot()+
  geom_density2d( aes(y=Z_red, x= Z_green),alpha =0.5, size=0.6, color="black", binwidth = 0.03 )


facs_me <- rbindlist(map(meth_cuts_reads, function(x){
  
  x[, .(me = sum(me), tot = sum(tot)), .(V1, chr, bin = round(bp / 2e4))
  ][mod_bin[,.SD[N == max(N)], .(chr, bin)], on = c("chr" , "bin")
  ][str_detect(V1, mod)]
  
})
)[,.(me = mean(me / tot)), .(V1)]

facs_me[, green := setNames(log(facs$`*[488] 530/40`), facs$V1)[V1]
][, red := setNames(log(facs$`*[561] 585/29`), facs$V1)[V1]
][, rank := setNames((facs$rank), facs$V1)[V1]
][, mod := str_extract(V1, "k[2-9]{1,2}me3")
]

facs_me[me > quantile(me, 0.95)]


facs_me1 <- facs_me %>% mutate(Z_green = scale(green),
                   Z_red = scale(red))

facs_me2 <- facs_me1 %>% 
  arrange(rank)%>%
  group_by(mod)%>%
  mutate(smooth_me = roll_median(me, 25))%>% ungroup()%>% as.data.table()
A1 <- facs_me2[mod == "k27me3"][order(abs(scale(smooth_me)[,1]))
] %>%
  ggplot() + 
  geom_point(data=dta_small, aes(x=Z_green, y= Z_red),alpha =0.33, size=0.1, color="grey80")+
  geom_point(aes(x=Z_green, y= Z_red, fill=smooth_me),shape =21, stroke=0.25, size =1.5)+
  scale_fill_distiller(palette = "RdYlBu", name= "CpG Methylation"#, limits = c(0.6, 0.8), oob = scales::squish
                       )+ 
  coord_cartesian()+ggtitle('Fig2a,\n10 NN, Winsorized, H3K27me3')+theme_bw()#+facet_grid(cols = vars(mod))+theme_bw()
  

A2 <- facs_me2[mod == "k9me3"][order(abs(scale(smooth_me)[,1]))
] %>%
  ggplot() + 
  geom_point(data=dta_small, aes(x=Z_green, y= Z_red),alpha =0.33, size=0.1, color="grey80")+
  geom_point(aes(x=Z_green, y= Z_red, fill=smooth_me),shape =21, stroke=0.25, size =1.5)+
  scale_fill_distiller(palette = "RdYlBu", name= "Methylation")+ coord_cartesian()+ggtitle('Fig2a,\n10NN, Winsorized, H3K9me3')+theme_bw()#+facet_grid(cols = vars(mod))+theme_bw()

A3 <- facs_me2[mod == "k36me3"][order(abs(scale(smooth_me)[,1]))
] %>%
  ggplot() + 
  geom_point(data=dta_small, aes(x=Z_green, y= Z_red),alpha =0.33, size=0.1, color="grey80")+
  geom_point(aes(x=Z_green, y= Z_red, fill=smooth_me),shape =21, stroke=0.25, size =1.5)+
  scale_fill_distiller(palette = "RdYlBu", name= "Methylation", limits = c(0.78, 0.82), oob = scales::squish
                       )+ coord_cartesian()+ggtitle('Fig2a,\n10NN, Winsorized, H3K36me3')+theme_bw()#+facet_grid(cols = vars(mod))+theme_bw()

A <- A1 / A3 / A2

A
library(ggh4x)

facs_me2 %>%
  ggplot()+
  geom_point(aes(y=me, x=rank/3, color=mod), alpha=0.5, size=1.2)+
  geom_line(aes(y=smooth_me, rank/3, color=mod), size=1.5)+
  theme_bw()+
  scale_color_manual(values=c("#D674BA","#86BF88","#7A416A"))+
  ggh4x::facet_grid2(. ~ mod, scales = "free_y", independent = "y")+
  xlab("Integrated Cell Cycle Progression")+
  ylab("CpG methylation per cell")+
  theme(legend.position = "top")

A1nn <- A1nn+xlab("hGem-mAG (log10 A.U.)")+ylab("hCdt1-mKO (log10 A.U.)")
A2NN <- A2NN+xlab("hGem-mAG (log10 A.U.)")+ylab("hCdt1-mKO (log10 A.U.)")
A5NN <- A5NN+xlab("hGem-mAG (log10 A.U.)")+ylab("hCdt1-mKO (log10 A.U.)")
A10NN <-  A10NN+xlab("hGem-mAG (log10 A.U.)")+ylab("hCdt1-mKO (log10 A.U.)")
A10_winsor <- A10_winsor+xlab("hGem-mAG (log10 A.U.)")+ylab("hCdt1-mKO (log10 A.U.)")
A1nn|A2NN|A5NN|A10NN|A10_winsor

facs_me

B1 <- me_avg[rank_bin < 86] %>%
  ggplot() +
  geom_raster(aes(x = rank2, y = rank_bin, fill = me_lo2)) +
  scale_fill_distiller(palette = "RdYlBu", name= "Methylation")+
  facet_grid(rows = vars(mod)) +
  coord_cartesian(expand = F)+
  theme(legend.position = "bottom")+
  theme(panel.background = element_rect(fill = "white"), panel.grid = element_blank()) 

B2 <- me_avg %>% group_by(mod)%>%
  mutate(tot1 = tot/max(tot))%>%
  group_by(rank_bin, mod) %>%
  summarize(tot1 = sum(tot1))%>%
  group_by(mod) %>%
  mutate(tot1 = tot1/max(tot1)) %>%
  ggplot() +
  geom_col(aes(x = tot1, y = rank_bin, fill = mod), WIDTH=1, orientation = "y")  +
  facet_grid(rows = vars(mod)) +
  scale_fill_manual(values=c("#D674BA","#86BF88","#7A416A"))+
  coord_cartesian(expand = F)+theme_bw()

B <- (( B1 | B2) + plot_layout(widths = c(4,1)))

B

library(scales)

momics        <- fread("~/Desktop/post-doc/ForkDOT/Fig1-cowplot/scEU_rpf_chic_Meth_psb_gw_smooth_2kb.tsv.gz")

C1 <- momics[chr == "1" & round < 4e7 #& cluster2 ==3
             ] %>%
  dplyr::filter(variable %in% c("k27me3", "k9me3", "k36me3") & value > 0 & smooth == 2)%>%
  mutate(round = round(round / 4e4) * 4e4)%>%
  group_by(round)%>% 
  dplyr::filter(value == max(value))%>% ungroup()%>% 
  #pivot_longer(cols = c(k27me3, k9me3, k36me3), names_to = "sample") %>% 
  ggplot()+
  # geom_segment(data = gene_location[chr == "2" & to < 5e7], 
  #               aes(x = from, xend = to, y = 0, yend = 0), size = 1.5) +
  geom_raster(aes(x = round, y = 1, fill = variable), alpha = 1, position= position_identity())+
  theme_bw()+
  ylab("Z-score")+
  scale_fill_manual(values=c("#D674BA","#86BF88","#7A416A"))+
  xlab("Chromosome 1 (bp)")+
  theme(legend.position="none")+coord_cartesian(expand=F)




C2 <- genome_example[bin > 5 & bin < 2000#][, cc_rank_cor := cc_rank2 - round((tmed * 0.4 + 0.5) * 500)][
] %>%
  {
    ggplot(data = .) +
      geom_raster(aes(x = bin, y = -ccrank2, fill = me_lo2)) + 
      scale_fill_distiller(palette = "RdYlBu", 
                           name= "Methylation", 
                           limits = c(-0.55, 0.35), oob = scales::squish)+
      #geom_point(aes(x = bin, y= 300, col = mod), size = 1)+
      geom_line(aes(x = merge_bin * 25, y = (tmed * 0.4 + 0.5) * -500), col = "black", alpha = 0.6) + 
      theme(panel.background = element_rect(fill = "white"), panel.grid = element_blank()) +
      coord_cartesian(expand = F) 
    
  }

(C1/C2)+plot_layout(heights = c(1,10))


A / B / C




RepTiming = list(
  fread("~/Desktop/post-doc/ForkDOT/Fig1-cowplot/4D_early_8R_UC.tsv.gz")%>% setNames(c("chr", "bp"))%>% mutate(round = round(bp / 10e4) * 10e4,
                                                                                                               RT= "early") %>% 
    group_by(chr, round, RT) %>%
    summarize(count = dplyr::n()),
  
  fread("~/Desktop/post-doc/ForkDOT/Fig1-cowplot/4D_late_3O_JH.tsv.gz")%>% setNames(c("chr", "bp"))%>% mutate(round = round(bp / 10e4) * 10e4,
                                                                                                              RT= "late") %>% 
    group_by(chr, round, RT) %>%
    summarize(count = dplyr::n())
)%>% rbindlist(idcol = F)

RT1=RepTiming %>% dplyr::filter(chr == "chr1", round>1.4e7 & round<2.5e7)%>%
  pivot_wider(names_from = RT, values_from = count, values_fill = 0)%>% mutate(total = early+late)%>%
  dplyr::filter(total>10)%>%
  mutate(reptiming= log2(early/late))%>% ggplot()+geom_raster(aes(x=round, y=0, fill=reptiming))+theme_bw()+
  scale_fill_distiller(palette = "RdYlBu", name= "Replication Timing")+coord_cartesian(expand = F)

RT2=fread("~/Desktop/post-doc/Experiments/datasets_RPE1/Halazonetis_Nature_RPE1_origins/RPE1.tsv.gz") %>% 
  rename(chr= V1, bp = V2) %>% mutate(RT = "origins") %>% group_by(chr) %>%
  mutate(round = round(bp / 5e4) * 5e4) %>% 
  group_by(chr, round) %>%
  summarize(count = dplyr::n())%>% dplyr::filter(chr==1, round>1.4e7, round<2.5e7)%>%
  mutate(count = (count - quantile(count,0.02)) / (quantile(count, 0.98) - quantile(count,0.02)))%>%
  ggplot()+geom_col(aes(x=round, y=count), fill="red")+theme_bw()+
  coord_cartesian(expand = F)+ylim(0,1.5)

library(round)
cuttracks = rbindlist(map(list.files("~/Desktop/post-doc/ChiC-TAPS/Fig2/sphase_order/"), 
                          function(x){
                            fread(paste0("~/Desktop/post-doc/ChiC-TAPS/Fig2/filtered_cells/", sub("_.*", ".tsv.gz", x))
                            )[chr == 1 &  posterior > 975
                            ][, rank := as.data.table(readRDS(paste0("~/Desktop/post-doc/ChiC-TAPS/Fig2/sphase_order/", x)
                            )$final_order)[, .(rank = unique(rank)), .(cell)][,setNames(rank, cell)][as.character(cell)]
                            ][, .(N = as.double(.N), exp = sub("_.*", "", x)), .( rank, bp = round_r3(bp, d = -5))
                            ][, N := (N) / max((N)), .(rank)][bp>1.4e7 & bp < 2.5e7]
                          })
)[str_detect(exp, "075", T), rank := abs(rank - 1)][
  str_detect(exp, "30WO", T), rank := abs(rank - 1)]%>% 
  ggplot() + 
  geom_raster(aes(y = rank, x = bp, fill = N) ) +
  #geom_point(aes(y = rank, x = bp), size = 0.5) +
  facet_grid(rows=vars(exp), scales = "free") +
  #scale_fill_viridis_c(option="rocket")+
  #scale_fill_distiller(palette = "Blues", name= "Tour Order") +
  scale_fill_gradientn( colors = c("white", "blue3", "blue3"), 
                        trans = "log", 
                        # limits = c(0.05, 0.15), oob = scales::squish,
                        name = "Cuts density") +
  theme_bw()+
  theme(   panel.background = element_rect(fill = "white"), 
           legend.position="none",
           panel.grid = element_blank()) +
  coord_cartesian(expand = F)+scale_y_reverse()
(RT1/RT2/cuttracks)+plot_layout(heights = unit(c(1,1,10), 'cm'))

(RT1/RT2/cuttracks)+plot_layout(heights = c(1,1,10))



facs <- fread("~/Desktop/post-doc/ChiC-TAPS/Fig2/facs.csv")

facs %>%
  ggplot() + geom_point(aes(y = log(`*[561] 585/29`), x = log(`*[488] 530/40`), 
                            fill=progression), size=2, alpha=0.66, shape=21)+theme_bw()+
  scale_fill_viridis_c(option="H", name = "Integrated \nCell Cycle\nProgression")+
  xlab("mAG-Geminin Fluorescent Intensity (log2))")+
  ylab("mKO-Cdt1 Fluorescent Intensity (log2)")
  



facs %>% 
  mutate(cellcycle = case_when(
    log(`*[488] 530/40`) < 3.4 ~ "G1",
    log(`*[561] 585/29`) > 2.66 &   log(`*[488] 530/40`) > 3.4 ~ "S-phase",
    log(`*[561] 585/29`) < 2.66 &   log(`*[488] 530/40`) > 3.4 ~ "G2")) %>%
  ggplot() + geom_point(aes(y = log(`*[561] 585/29`), x = log(`*[488] 530/40`),fill= cellcycle), size=2,alpha=0.75, shape=21)+
  theme_bw()+
  scale_fill_manual(values=c("red2","green3", "yellow2"))+
  theme(legend.position = "none")+
  xlab("mAG-Geminin Fluorescent Intensity (log2))")+
  ylab("mKO-Cdt1 Fluorescent Intensity (log2)")

facs %>% 
  mutate(cellcycle = case_when(
    log(`*[488] 530/40`) < 3.4 ~ "G1",
    log(`*[561] 585/29`) > 2.7 &   log(`*[488] 530/40`) > 3.4 ~ "S-phase",
    log(`*[561] 585/29`) < 2.7 &   log(`*[488] 530/40`) > 3.4 ~ "G2")) %>%
  ggplot() + geom_histogram(aes(x = `*[405] 460/50`, fill = cellcycle),
                            position = position_identity(), bins = 75, alpha = 0.8)+theme_bw()+
  scale_fill_manual(values=c("red2","green3", "yellow2"))+
  xlab("DAPI Fluorescent Intensity (A.U.)")+
  theme(legend.position = "none")




facs %>% 
  mutate(cellcycle = case_when(
    log(`*[488] 530/40`) < 3.4 ~ "G1",
    log(`*[561] 585/29`) > 2.7 &   log(`*[488] 530/40`) > 3.4 ~ "S-phase",
    log(`*[561] 585/29`) < 2.7 &   log(`*[488] 530/40`) > 3.4 ~ "G2")) %>%
  drop_na() %>%
  mutate(progression = frank(progression, ties.method = "dense")) %>%
  ggplot() + geom_tile(aes(x = progression/3, y=1, fill = cellcycle), position= position_dodge())+
  theme_bw()+
  scale_fill_manual(values=c("red2","green3", "yellow2"))



fread("~/Desktop/post-doc/ChiC-TAPS/Fig2/all_cells_sum.tsv")[, mod := str_extract(V1, "k[2-9]{1,2}me[1-3]{1}")][] %>%
  mutate(cell = str_split(V1, "_", simplify = TRUE)[, 2]) %>%
  mutate(emptywell = if_else((cell %in% c("1", "25", "49", "73")), "empty",  "cell")) %>%
  ggplot()+
  geom_boxplot(aes(x=emptywell, y=N, fill=mod), alpha=0.8, outlier.shape=NA)+
  scale_fill_manual(values=c("#D674BA","#86BF88","#7A416A"))+
  theme_bw()+
  facet_grid(cols=vars(mod))+
  ylim(0,1e5)+scale_y_log10()+
  ylab("Unique ChIC cuts \nper cell (log10)")
  

