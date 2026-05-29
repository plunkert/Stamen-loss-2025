# Investigate putative linked QTL on chromosomes 1 and 5 in Tsu Kas

setwd("~/Documents/Github/Stamen-loss-2025/") # Change this

require(dplyr)
require(ggplot2)
require(ggpubr)

# plotPXG function is being wacky so I'm just plotting phenotype ~ genotype
# for these recombinants using tidyverse

# read in Tsu Kas genotype and phenotype data
tk.dat <- read.csv("./Data/tsu_kas_pheno_geno/Tsu_Kas_Rqtl_format.csv")[-c(1,2),]

# relevant markers from stepwise analysis:
# 1@80.0 : "C1_24316033"
# 1@108.0 : "X1_27650244"
# 5@2.8 : "X5.1213621"
# 5@9.8 : "X5_2799359"

tk.dat$mean_ss_num <- as.numeric(tk.dat$mean_ss_num)

# create column indicating genotype at both chr1 QTLs
tk.dat$chr1_geno <- paste(tk.dat$C1_24316033, tk.dat$X1_27650244, sep="_") %>% as.factor()

# create column indicating genotype at both chr5 QTLs
tk.dat$chr5_geno <- paste(tk.dat$X5.1213621, tk.dat$X5_2799359, sep = "_") %>% as.factor()

# check out factor levels and remove those that include unknown genotypes
table(tk.dat$chr1_geno)
table(tk.dat$chr5_geno)

# how many have an unknown genotype?
sum(!grepl(pattern="-", x=tk.dat$chr1_geno)) # both genotypes known for 129 RILs
sum(grepl(pattern="-", x=tk.dat$chr1_geno)) # at least one genotype unknown for 212 RILs

sum(!grepl(pattern="-", x=tk.dat$chr5_geno)) # both genotypes known for 306 RILs
sum(grepl(pattern="-", x=tk.dat$chr5_geno)) # at least one genotype unknown for 35 RILs

# make columns that have no unknown genotypes

tk.dat$chr1_geno_known <- case_when(!grepl(pattern="-", x=tk.dat$chr1_geno) ~ tk.dat$chr1_geno,
                                    .default = NA)

tk.dat$chr5_geno_known <- case_when(!grepl(pattern="-", x=tk.dat$chr5_geno) ~ tk.dat$chr5_geno,
                                    .default = NA)

# change AA BB to KK TT for ease of reading

tk.dat$chr1_geno_known <- tk.dat$chr1_geno_known %>% gsub(pattern="A", replacement="K")
tk.dat$chr1_geno_known <- tk.dat$chr1_geno_known %>% gsub(pattern="B", replacement="T")

tk.dat$chr5_geno_known <- tk.dat$chr5_geno_known %>% gsub(pattern="A", replacement="K")
tk.dat$chr5_geno_known <- tk.dat$chr5_geno_known %>% gsub(pattern="B", replacement="T")

# plot short stamen number as a function of chr1 genotype
chr1_boxplot = tk.dat %>%
  filter(!is.na(chr1_geno_known)) %>%
  ggplot() +
  aes(x = chr1_geno_known, fill = chr1_geno_known, y = mean_ss_num) +
  geom_boxplot(outliers = F) +
  geom_jitter(position=position_jitter(0.1)) +
  # Labels
  scale_x_discrete(name = 'Genotype at 1@80 and 1@108 QTLs') +
  scale_y_continuous(
    name = 'RIL Short Stamen Number', limits=c(0.3, 2.2)
  ) +
  # Style
  scale_fill_manual(values = c("#88CCEE", "#44AA99", "#DDCC77", "#882255")) +
  theme_minimal() +
  theme(
    legend.position = "none",
    axis.text = element_text(size=16),
    axis.title = element_text(size=16),
  )

chr5_boxplot = tk.dat %>%
  filter(!is.na(chr5_geno_known)) %>%
  ggplot() +
  aes(x = chr5_geno_known, fill = chr5_geno_known, y = mean_ss_num) +
  geom_boxplot(outliers = F) +
  geom_jitter(position=position_jitter(0.1)) +
  # Labels
  scale_x_discrete(name = 'Genotype at 5@2.8 and 5@9.8 QTLs') +
  scale_y_continuous(
    name = 'RIL Short Stamen Number', limits=c(0.3, 2.2)
  ) +
  # Style
  scale_fill_manual(values = c("#88CCEE", "#44AA99", "#DDCC77", "#882255")) +
  theme_minimal() +
  theme(
    legend.position = "none",
    axis.text = element_text(size=16),
    axis.title = element_text(size=16),
  )

# let's investigate the effects of these loci using ANOVA and Tukey's test
TukeyHSD(aov(data=tk.dat, mean_ss_num ~ chr1_geno_known))

# letters on chr1 boxplot are a, b, b, b
TukeyHSD(aov(data=tk.dat, mean_ss_num ~ chr5_geno_known))
# letters on chr5 boxplot are a, b, b, b

svg("./Results/linked_qtl_boxplots.svg", width=12, height=7)
ggarrange(chr1_boxplot, chr5_boxplot)
dev.off()
