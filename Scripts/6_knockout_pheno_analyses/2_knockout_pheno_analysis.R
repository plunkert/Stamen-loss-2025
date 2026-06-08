# Arabidopsis generation 2 stamen phenotype analysis
rm(list=ls())
library(dplyr)
library(lme4)
library(lubridate)
library(readxl)

setwd("~/Documents/GitHub/Stamen-loss-2025/")

dat <- read_excel("./Data/knockout_phenos/2_knockout_stamen_pheno.xlsx", sheet="Sheet1")
dat$genotype <- as.factor(substr(dat$plant_id, 1, 1))
dat$plant_id <- as.factor(dat$plant_id)
dat$collection <- as.factor(dat$collection)
# Check that each plant id has the expected number of rows (9 rows, 3 for each collection)
table(dat$plant_id) # looks good

# Check that each genotype has the expected number of rows (36, 3 for each collection,
# 3 collections per plant, 4 plants per genotype but more for Col-0 which is A)
table(dat$genotype) # looks good

# function blockSubset subsets the dat dataframe to get rows corresponding to 
# each experimental block (4 mutant plants and their corresponding 4 WT plants 
# they were randomized with). Takes a mutant line genotype (SALK_XXXXXX) and a vector of 
# 4 numbers indicating which Col plants were used

blockSubset <- function(mutantLine1, mutantLine2, mutantLine3, x, data){
  WT_stamen_counts <- filter(data, genotype=="A")
  wt <- filter(data, genotype=="A" & plant_id == x[1] | plant_id == x[2] | plant_id == x[3] | plant_id == x[4])
  mut1 <- filter(data, genotype==mutantLine1)
  mut2 <- filter(data, genotype==mutantLine2)
  mut3 <- filter(data, genotype==mutantLine3)

  df <- rbind(wt, mut1, mut2, mut3)
}

# Check whether flower number or collection affect short stamen number across the dataset
out.flower <- lm(data=dat, short ~ flower_pos)
summary(out.flower)
out.collection <- lm(data=dat, short~collection)
summary(out.collection)

# Make separate dataframes for each block (i.e. each tray). For blocks that don't
# have 3 mutant lines, set mutantLine3 = "" so that function still works
block1 <- blockSubset(mutantLine1="B", mutantLine2="C", mutantLine3="D", x = c("A01", "A02", "A03", "A04"), data=dat)
block2 <- blockSubset(mutantLine1="E", mutantLine2="F", mutantLine3="", x = c("A05", "A06", "A07", "A08"), data=dat)
block3 <- blockSubset(mutantLine1="H", mutantLine2="I", mutantLine3="", x = c("A09", "A10", "A11", "A12"), data=dat)
block4 <- blockSubset(mutantLine1="K", mutantLine2="L", mutantLine3="", x = c("A13", "A14", "A15", "A16"), data=dat)

# Get mean short stamen number for each plant ID in each block and put in dataframe
# blockx_means. Add a genotype column to dataframes for mean values of each plant.
block1_means <- aggregate(data=block1, short ~ plant_id, FUN=mean)
block1_means$genotype <- as.factor(substr(as.character(block1_means$plant_id), 1, 1))
block2_means <- aggregate(data=block2, short ~ plant_id, FUN=mean)
block2_means$genotype <- as.factor(substr(as.character(block2_means$plant_id), 1, 1))
block3_means <- aggregate(data=block3, short ~ plant_id, FUN=mean)
block3_means$genotype <- as.factor(substr(as.character(block3_means$plant_id), 1, 1))
block4_means <- aggregate(data=block4, short ~ plant_id, FUN=mean)
block4_means$genotype <- as.factor(substr(as.character(block4_means$plant_id), 1, 1))

# Linear models for each block (m1 for block 1, etc.) where each genotype is a level
# of the predictor and each plant's mean short stamen # is the response.
m1 <- lm(data=block1_means, short ~ genotype)
summary(m1)

m2 <- lm(data=block2_means, short ~ genotype)
summary(m2)

m3 <- lm(data=block3_means, short ~ genotype)
summary(m3)

m4 <- lm(data=block4_means, short ~ genotype)
summary(m4)

# mean Col-0 short stamen number for this experiment
wt <- dat %>% filter(genotype == "A")
mean(wt$short) # 1.94
