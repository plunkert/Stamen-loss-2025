# This script produces a figure of the candidate genes in and around the QTL region on 
# chromosome 5

library(readxl)
library(gggenes)
library(tidyverse)

setwd("~/Documents/GitHub/Stamen-loss-2025/Data/")

dat <- read_excel("cand_genes_gggenes_input.xlsx", sheet="Genes")
#vars <- read_excel("cand_genes_gggenes_input.xlsx", sheet="Variants")

theme_set(theme_bw(base_size=10))
theme_update(panel.grid.major = element_blank(),
             panel.grid.minor = element_blank(), 
             strip.background = element_blank(),
             legend.position = "none")

ggplot(dat, aes(xmin = start, xmax = end, y = molecule, fill = type, color = type))+
 geom_gene_arrow(arrowhead_width=grid::unit(0, "mm"), arrowhead_height=grid::unit(0, "mm")) + 
#  geom_feature(
#    data = vars,
#    aes(x = position, y = molecule))+
  geom_blank(data=dat)+
  scale_x_continuous(limits=c(800000, 3.5e6))+
  ggrepel::geom_text_repel(aes(x = end - ((end-start)/2), y = molecule, label = feature, angle = c(rep(90, 13),0,0,0)), size = 5, color = "black", 
                           nudge_y = c(rep(c(0.3,-0.3), 6),0.3,0.15,0.15,0.15))+
  guides(fill = F) +
  scale_color_manual(aesthetics = c("fill", "color"), values = c("gene" = "red", "qtl" = "lightseagreen")) +
  theme_genes() + theme(legend.position = "none",
                        axis.text = element_text(color="black", 
                                                   size=12))+
  labs(x="Genomic coordinate on TAIR10 chromosome 5 (bp)", y="", colour="black")

