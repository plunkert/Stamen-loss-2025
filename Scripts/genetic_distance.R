# Process pixy dxy output to get single values for the whole genome

# Load packages
library(ggplot2)
library(dplyr)
library(ggrepel)
setwd("~/Documents/GitHub/Stamen-loss-2025/")

dat <- read.delim("./Data/all_filtered.50k.out_dxy.txt", header=TRUE)
dat <- read.delim("./Data/all_filtered.50k.out.nobqsr_dxy.txt", header=TRUE)


dat$pop_pair <- paste(dat$pop1, dat$pop2, sep="_") %>% as.factor()

sums <- dat %>% group_by(pop_pair) %>% summarize(Sum.count.diff = sum(count_diffs, na.rm = T), Sum.count.comp = sum(count_comparisons, na.rm = T))

sums <- dat %>% filter(!is.na(count_comparisons) & !is.na(count_diffs)) %>% group_by(pop_pair) %>% 
  summarize(Sum.count.diff = sum(count_diffs, na.rm = T), Sum.count.comp = sum(count_comparisons, na.rm = T))

sums$Mean.dxy <- sums$Sum.count.diff/sums$Sum.count.comp
sums

means <- dat %>% group_by(pop_pair) %>% summarize(genome_wide_dxy = mean(avg_dxy, na.rm=T))
means


# Try PCA instead

eigenValues_test <- read.delim("./Data/plinkPCA_Var_invar.eigenval", sep = " ", header = F)
eigenVectors_test <- read.delim("./Data/plinkPCA_Var_invar.eigenvec", sep = " ", header = F)
# each row is a genotype, each column is a PC.
colnames(eigenVectors_test) <- c("Line", "SeqID", paste0("PC", 1:(ncol(eigenVectors_test)-2)))

# calculate the PVE
pve_test <- round((eigenValues_test / (sum(eigenValues_test))*100), 2)
pve_test$PC <- c(1:4)

# add pop column
eigenVectors_test$Pop <- toupper(substr(eigenVectors_test$Line, 1, 3))
ggplot(pve_test)+ geom_bar(aes(y = V1,x = PC), stat = "identity")

pca_plot <- ggplot(data = eigenVectors_test) +
  geom_jitter(mapping = aes(x = -(PC1), y = PC2, shape = Pop), col = "black", size = 1.3, stroke = 0.5, width = 0.025, height=0.025, show.legend = TRUE ) +
  geom_hline(yintercept = 0, linetype="dotted") +
  geom_vline(xintercept = 0, linetype="dotted") +
  labs(x = paste0("PC 1 (",sprintf(pve_test[1,1], fmt='%#.1f'),"%)"),
       y = paste0("PC 2 (",sprintf(pve_test[2,1], fmt='%#.1f'),"%)")) +
  scale_shape_manual(name = "Population",
                     labels = eigenVectors_test$Pop,
                     values = rep(c(22, 21, 24, 23, 25), times = 4)) +
  theme_classic()+
  theme(
    legend.title = element_text(color = "black", size = 12),
    legend.text = element_text(color = "black", size = 12),
    axis.title.y = element_text(color = "black", size = 10),
    axis.title.x = element_text(color = "black", size = 10),
    axis.text = element_text(color = "black", size = 10),
    legend.spacing.y = unit(0.03, "cm"))+
  geom_text_repel(data = eigenVectors_test,aes(x = -(PC1), y = PC2, label=Pop), nudge_x=0.15, fontface = "bold", segment.color = 'transparent', size=4.5, inherit.aes = FALSE)
  


