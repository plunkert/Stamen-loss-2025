# This script checks for differences in short stamen number between mutant lines
# and the Col-0 found in the same block. Third round of phenotyping done in 
# Spring 2025 with Esther.

rm(list=ls())
library(GGally)
library(dplyr)
library(ggplot2)

# read raw phenotype data into R
setwd("~/Documents/GitHub/Stamen-loss-2025/")
dat <- read_excel("./Data/stamen_cand_knockouts_pheno_2025.xlsx")
head(dat)

# check that each plant ID has 9 flowers
table(dat$plant_id)

# ones to check in paper datasheet bc far off from expected 9: 13 T2Col4, 6 T1A5, 5 T1Col4, 1 T2D6 (shouldn't exist)
# 7 T2D3, 10 T3E4

# make column indicating whether flower was from main stalk or side stalk
dat$stalk_type <- case_when(grepl("S", dat$flower_pos) ~ "S",
                            .default = "C")
# then remove S from flower_pos and make into an integer
dat$flower_pos <- gsub("S", "", dat$flower_pos) %>% as.numeric()

# get a mean short stamen number for each plant
dat_means <- aggregate(data=dat, short_stamens ~ plant_id, FUN=mean)

# Make new column for tray (i.e. block)
dat_means$tray <- substr(dat_means$plant_id, 1, 2)

# Make genotype column
dat_means$genotype <- case_when(grepl("Col", dat_means$plant_id) ~ "WT",
                          grepl("A", dat_means$plant_id) ~ "SALK_128849", # replace lowercase with actual line
                          grepl("B", dat_means$plant_id) ~ "SALK_103498",
                          grepl("C", dat_means$plant_id) ~ "SALK_138334",
                          grepl("D", dat_means$plant_id) ~ "CS926856",
                          grepl("E", dat_means$plant_id) ~ "CS923713",
                          grepl("F", dat_means$plant_id) ~ "SALK_098887C",
                          grepl("G", dat_means$plant_id) ~ "CS857772",
                          grepl("H", dat_means$plant_id) ~ "CS856723")

#dat$genotype <- as.factor(dat$genotype) %>% relevel(ref="WT")

filter(dat_means, dat_means$tray=="T2") %>%
  ggplot(aes(x=genotype, y=short_stamens, group=genotype)) + geom_boxplot() + 
  geom_jitter(width=0.1, height=0) +
  geom_hline(yintercept=0)+
  ggtitle("Short Stamen Number")

dat_means$genotype <- as.factor(dat_means$genotype) %>% relevel(ref="WT")

# fit linear model to test effect of genotype in each tray
m1 <- lm(data=filter(dat_means, dat_means$tray=="T1"), short_stamens ~ genotype)
summary(m1)

m2 <- lm(data=filter(dat_means, dat_means$tray=="T2"), short_stamens ~ genotype)
summary(m2)
# genotype SALK_138334 has marginally significant effect but it's increasing short stamen number
# other line for same gene (CS926856) didn't show same result

m3 <- lm(data=filter(dat_means, dat_means$tray=="T3"), short_stamens ~ genotype)
summary(m3)

m4 <- lm(data=filter(dat_means, dat_means$tray=="T4"), short_stamens ~ genotype)
summary(m4)

# What was the mean short stamen number for Col-0?
mean(dat_means[which(dat_means$genotype == "WT"), "short_stamens"])


# Do this again with LSM short stamen number so we consider flower position and such!

ss.out <- lm(data = dat, short_stamens ~ stalk_type + flower_pos + date + tray + genotype)
summary(ss.out)

# Residual diagnostics to see if data meet assumptions of lm
res.vals <- resid(ss.out)
pred.vals <- fitted(ss.out) # Fitted values
# Is there a pattern in the residuals?
plot(pred.vals, res.vals, main = "Residuals vs. pred.vals values", 
     las = 1, xlab = "Predicted values", ylab = "Residuals", pch = 19)
abline(h = 0)
# yes, a very wacky one, probably stamen num = 2 follows one pattern and stamen num=1 follows another

# Make qqnorm plot and see if points are along diagonal
qqnorm(res.vals, pch = 19)
qqline(res.vals, col = 'red') # Off diagonal at lowest quantiles, makes sense bc
# we suddenly jump from 2 to 1 short stamens

# calculate reference grid and estimated marginal means for stamen number model
ss.rg <- ref_grid(ss.out)
ss.e <- emmeans(ss.rg, "genotype")
summary(ss.e)

# okay this gives us LSMs for each genotype, but residual diagnostics are wacky
# enough that I don't trust it. Let's try fitting binomial dist. model of 
# whether stamen loss occurred instead.

dat$loss <- dat$short_stamens < 2

out.ss.bin <- glm(family="binomial", data=dat, dat$loss ~ genotype + tray)
summary(out.ss.bin)

# model containing only tray 2
out.ss.bin.t2 <- glm(family="binomial", data=dat[which(dat$tray=="T2"),], 
                     loss ~ genotype + flower_pos + stalk_type)
summary(out.ss.bin.t2)

# Conclusion that line C has more short stamens (by general linear model) is not robust to
# changes in analysis (e.g. modeling using binomial distribution)
