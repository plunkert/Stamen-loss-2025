# 9/23/25. This script checks for differences in short stamen number between mutant lines
# and the Col-0 found in the same block. Fourth round of phenotyping done in 
# Summer/Fall 2025 with Esther.

# This experiment included estradiol-inducible FLC lines, but RT-qPCR showed that we failed
# to induce FLC expression so these lines are not discussed in the manuscript since not informative
# about whether FLC induction would have led to stamen loss.

rm(list=ls())
library(GGally)
library(dplyr)
library(ggplot2)
library(emmeans)
library(readxl)

# read raw phenotype data into R
setwd("~/Documents/GitHub/Stamen-loss-2025/")
dat <- read_excel("./Data/knockout_phenos/FLC_stamen_count_data_August_2025.xlsx", sheet="Stamen Counts")
trt <- read_excel("./Data/FLC_stamen_count_data_August_2025.xlsx", sheet="Treatments")
head(dat)


# check that each plant ID has 9 flowers
table(dat$plant_id)
# G12 and G3 had estradiol applied only after bolting and G14 is missing one tube

# make column indicating whether flower was from main stalk or side stalk
dat$stalk_type <- case_when(grepl("S", dat$flower_pos, ignore.case=TRUE) ~ "S",
                            .default = "C")
# then remove S from flower_pos and make into an integer
dat$flower_pos <- gsub("S", "", dat$flower_pos, ignore.case=TRUE) %>% as.numeric()

# indicate treatment by merging with treatment sheet
dat <- merge(dat, trt, by.x="plant_id", by.y="Plant", all.x=TRUE)


# Make genotype column
dat$genotype <- case_when(grepl("A", dat$plant_id) ~ "Col-0",
                                grepl("C", dat$plant_id) ~ "SALK_0725906C",
                                grepl("D", dat$plant_id) ~ "SALK_041126C",
                                grepl("E", dat$plant_id) ~ "CS2102347",
                                grepl("F", dat$plant_id) ~ "CS2102348",
                                grepl("G", dat$plant_id) ~ "CS2102349")


# What was the mean short stamen number for Col-0 in this run?
mean(dat[which(dat$genotype == "Col-0"), "short_stamens"])

# Indicate block (tray) for KO lines not included in treatment spreadsheet
dat$block[is.na(dat$block)] <- "pink"

# Indicate genotype/treatment combos
dat$geno_trt <- paste(dat$genotype, dat$Treatment, sep="_") %>% as.factor()

# remove bolt estradiol treatment, didn't intend to phenotype
dat <- filter(dat, !grepl("bolt_estradiol", geno_trt))

# calculate mean short stamen number for each plant
dat_means <- dat %>% group_by(plant_id) %>% 
  summarize(mean_ss_num = mean(short_stamens), genotype = genotype[1], geno_trt=geno_trt[1], block=block[1])

dat_means$genotype <- relevel(as.factor(dat_means$genotype), ref="Col-0")
dat_means$geno_trt <- relevel(as.factor(dat_means$geno_trt), ref="Col-0_mock")

# subset dataframe of means into each block to be analyzed separately. Each line gets compared to its
# respective Col-0 in the same block
dat_means_blue <- dat_means[which(dat_means$block=="blue"),]
dat_means_red <- dat_means[which(dat_means$block=="red"),]
dat_means_yellow <- dat_means[which(dat_means$block=="yellow"),]
dat_means_pink <- dat_means[which(dat_means$block=="pink"),]

# boxplot of short stamen number in each block

pink_plot <- dat_means_pink %>% ggplot(aes(x=genotype, y=mean_ss_num, group=genotype)) + geom_boxplot(outliers=F) + 
  geom_jitter(width=0.1, height=0) +
  ylab("short stamen number")+ ylim(c(0,2.5))+ xlab("")+
  theme(text=element_text(size=14), axis.text.x = element_text(angle = 30, hjust = 1))

blue_plot <- dat_means_blue %>% ggplot(aes(x=geno_trt, y=mean_ss_num, group=geno_trt)) + geom_boxplot(outliers=F) + 
  geom_jitter(width=0.1, height=0) +
  ylab("short stamen number")+ ylim(c(0,2.5))+ xlab("")+
  theme(text=element_text(size=14), axis.text.x = element_text(angle = 30, hjust = 1))


red_plot <- dat_means_red %>% ggplot(aes(x=geno_trt, y=mean_ss_num, group=geno_trt)) + geom_boxplot(outliers=F) + 
  geom_jitter(width=0.1, height=0) +
  ylab("short stamen number")+ ylim(c(0,2.5))+ xlab("")+
  theme(text=element_text(size=14), axis.text.x = element_text(angle = 30, hjust = 1))


yellow_plot <- dat_means_yellow %>% ggplot(aes(x=geno_trt, y=mean_ss_num, group=geno_trt)) + geom_boxplot(outliers=F) + 
  geom_jitter(width=0.1, height=0) +
  ylab("short stamen number")+ ylim(c(0,2.5))+ xlab("")+
  theme(text=element_text(size=14), axis.text.x = element_text(angle = 30, hjust = 1))

ggarrange(pink_plot, blue_plot, red_plot, yellow_plot, ncol=2, nrow=2)



m_pink <- lm(data=dat_means_pink, mean_ss_num ~ genotype)
summary(m_pink)

m_blue <- lm(data=dat_means_blue, mean_ss_num ~ geno_trt)
summary(m_blue)

m_red <- lm(data=dat_means_red, mean_ss_num ~ geno_trt)
summary(m_red)

m_yellow <- lm(data=dat_means_yellow, mean_ss_num ~ geno_trt)
summary(m_yellow)
