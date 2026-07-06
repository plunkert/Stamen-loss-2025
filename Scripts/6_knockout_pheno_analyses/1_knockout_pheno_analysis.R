# Arabidopsis round 1 stamen phenotype analysis
# 05/17/2022 - This analysis is a simple ANOVA, like stamen-gen2-analysis.R
# but with slight modification because of different experimental design

#load required packages
rm(list=ls())
library(dplyr)
library(lme4)
library(lubridate)
library(readxl)


setwd("~/Documents/GitHub/Stamen-loss-2025/")

# read in stamen phenotyping data
dat <- read_excel("./Data/knockout_phenos/1_knockout_stamen_pheno.xlsx",sheet="Sheet1")

# since SALK_089563 wasn't genotyped properly and was done in next generation, 
# remove it from dataset
dat <- dat[dat$Line!="SALK_089563",]

# create plant_id (genotype-replicate number) for each plant
dat$plant_id <- paste(as.character(dat$Line), as.character(dat$Plant), sep="-")

# SALK_012835-9 was included in phenotyping by accident so remove it
dat <- dat[dat$plant_id != "SALK_012835-9",]
table(dat$Line)

# should be 9 flowers for each plant_id, true for all except Col-14 which has only 6
table(dat$plant_id)

# number the flower collections since I didn't indicate this during data collection
# add collection numbers (1, 2, 3 indicating first, second, third collection)

dat_coll <- dat %>% mutate(CollectionDate = as.Date(CollectionDate)) %>%
  group_by(plant_id) %>%
  mutate(CollectionNum = case_when(CollectionDate == min(CollectionDate) ~ 1, 
                                CollectionDate == max(CollectionDate) ~ 3,
                                .default = 2)) %>% ungroup()

dat_coll$plant_id <- as.factor(dat_coll$plant_id)
# function to pull out only the plants corresponding to one block (4 WT + 4 mutant)
blockSubset <- function(mutantLine, x, data){
  wt <- filter(data, Line=="Col" & (Plant == x[1] | Plant == x[2] | Plant == x[3] | Plant == x[4]))
  mut <- filter(data, Line==mutantLine)
  
  df <- rbind(wt, mut)
}

# subset dat2 into each block of mutants and the corresponding 4 wild-types
block1 <- blockSubset(mutantLine="SALK_012835", x = c(1,9,12,16), data=dat_coll)
block2 <- blockSubset(mutantLine="SALK_028815", x=c(2,10,13,17), data=dat_coll)
block3 <- blockSubset(mutantLine = "SALK_055242", x=c(3,6,14,18), data=dat_coll)
block4 <- blockSubset(mutantLine = "SALK_065777", x=c(4,7,15,19), data=dat_coll)



# Get mean short stamen number for each plant ID in each block and put in dataframe
# blockx_means. Add a genotype column to dataframes for mean values of each plant,
# which is just S (mutant) or C (WT)
block1_means <- aggregate(data=block1, ShortStamenNum ~ plant_id, FUN=mean)
block1_means$genotype <- as.factor(substr(as.character(block1_means$plant_id), 1, 1))
block2_means <- aggregate(data=block2, ShortStamenNum ~ plant_id, FUN=mean)
block2_means$genotype <- as.factor(substr(as.character(block2_means$plant_id), 1, 1))
block3_means <- aggregate(data=block3, ShortStamenNum ~ plant_id, FUN=mean)
block3_means$genotype <- as.factor(substr(as.character(block3_means$plant_id), 1, 1))
block4_means <- aggregate(data=block4, ShortStamenNum ~ plant_id, FUN=mean)
block4_means$genotype <- as.factor(substr(as.character(block4_means$plant_id), 1, 1))

# Linear models for each block (m1 for block 1, etc.) where each genotype is a level
# of the predictor and each plant's mean short stamen # is the response.
m1 <- lm(data=block1_means, ShortStamenNum ~ genotype)
summary(m1)

m2 <- lm(data=block2_means, ShortStamenNum ~ genotype)
summary(m2)

m3 <- lm(data=block3_means, ShortStamenNum ~ genotype)
summary(m3)

m4 <- lm(data=block4_means, ShortStamenNum ~ genotype)
summary(m4)

# mean Col-0 short stamen number for this experiment
wt <- dat_coll %>% filter(Line == "Col")
mean(wt$ShortStamenNum) # 1.95




