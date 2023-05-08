library(tidyverse)
library(data.table)
library(patchwork)


setwd('~/Desktop/post-doc/ChiC-TAPS/')

loess_wrap <- function(x, y, sspan){
  
  X <- data.table(x2 = c(x - max(x), x, x + max(x)), 
                  y2 = rep(y, 3))
  
  predict(loess(y2 ~ x2, data = X, span = sspan), x)
}

# per cell summary

meth_cuts_sum <- fread("Documents/Projects/chic_taps/data/RPE/all_cells_sum.tsv")

facs <- fread("Downloads/facs.csv"
)[, mod := str_extract(V1, "k[2-9]{1,2}me3")
][, rank := frank(progression + runif(.N, max = 0.00001), ties.method = "dense"), .(mod)]

sphase <- readRDS("Documents/Projects/Jerry/JvB-exp049-15m-EdU-beads_sphase.RDS")$final_order %>% distinct(cell, rank)

RT <- fread("Documents/Projects/Jerry/JvB-exp049-15m-EdU-beads_single_fork.tsv.gz"
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

meth_cuts_reads <- readRDS("Documents/Projects/chic_taps/data/RPE/meth_cuts_reads.RDS")#[[1]][, max(bp) - min(bp), .(V1, chr, bin)]

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

A1 <- facs_me[mod == "k27me3"][order(abs(scale(me)[,1]))
] %>%
  ggplot() + geom_point(aes(x = green, y = red, fill = me), size = 2, col = "black", shape = 21) +
  scale_fill_gradientn(colours = c("blue", "white", "red")) +
  facet_grid(cols = vars(mod)) + theme_bw()

A2 <- facs_me[mod == "k9me3"][order(abs(scale(me)[,1]))
] %>%
  ggplot() + geom_point(aes(x = green, y = red, fill = me), size = 2, col = "black", shape = 21) +
  scale_fill_gradientn(colours = c("blue", "white", "red")) +
  facet_grid(cols = vars(mod)) + theme_bw()

A3 <- facs_me[mod == "k36me3"][order(abs(scale(me)[,1]))
] %>%
  ggplot() + geom_point(aes(x = green, y = red, fill = me), size = 2, col = "black", shape = 21) +
  scale_fill_gradientn(colours = c("blue", "white", "red")) +
  facet_grid(cols = vars(mod)) + theme_bw()

A <- A1 | A2 | A3





B1 <- me_avg[rank_bin < 86] %>%
  ggplot() +
  geom_raster(aes(x = rank2, y = rank_bin, fill = me_lo2)) +
  scale_fill_viridis_c() +
  facet_grid(rows = vars(mod)) +
  coord_cartesian(expand = F) 

B2 <- me_avg %>%
  ggplot() +
  geom_col(aes(x = tot, y = rank_bin), orientation = "y")  +
  facet_grid(rows = vars(mod)) +
  coord_cartesian(expand = F)

B <- B1 | B2

C <- genome_example[bin > 5 & bin < 2000#][, cc_rank_cor := cc_rank2 - round((tmed * 0.4 + 0.5) * 500)][
] %>%
  {
    ggplot(data = .) +
      geom_raster(aes(x = bin, y = -ccrank2, fill = me_lo2)) + 
      scale_fill_viridis_c() +
      #geom_point(aes(x = bin, y= 300, col = mod), size = 1)+
      geom_line(aes(x = merge_bin * 25, y = (tmed * 0.4 + 0.5) * -500), col = "red", alpha = 1) + 
      theme(panel.background = element_rect(fill = "black"), panel.grid = element_blank()) +
      coord_cartesian(expand = F) 
    
  }
A / B / C