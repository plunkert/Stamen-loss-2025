# Wrangle phenotype data, normalize, and produce histograms
# Make correlation tables for phenotypes in Tsu Kas, Belm Roda RILs

# load required packages
require(ggplot2)
require(corrtable)
require(readxl)
require(qtl)
require(dplyr)
require(bestNormalize)

setwd("~/Documents/Github/Stamen-loss-2025/") # Change this

# Read in phenotype data

# read in ovule data (LSMs and raw means) from Jeff's JMP analysis
ov.dat <- read.csv("./Data/tsu_kas_pheno_geno/ovule_lsm.csv")

# read in plant id table
ids <- read.csv("./Data/tsu_kas_pheno_geno/TK_RILs_IDs.csv")

# add ID to ids w/o CS so they can be compared to ov.dat
ids$ril <- substr(ids$stock_num, 3, 7)
ov.dat$ril <- as.character(ov.dat$ril)

ov.dat.ids <- left_join(x = ov.dat, y = ids, by = "ril")

# Read in cross data and add raw and LSM ovule data
dat.tk <-read.cross("csv", dir="./Data/tsu_kas_pheno_geno/", file = "TKrilsGenos_map55_pheno.csv",
                 genotypes=c("a","b"), na.strings=c("NA","-"))


# read in flowering time data and merge by Line with ovule data
ft <- read_excel("./Data/tsu_kas_pheno_geno/FT_pheno_TK_Lovell.xlsx")
colnames(ft) <- c("Line", "wue", "ft")
ft.ov.dat <- merge(x=ft, y = ov.dat.ids, by="Line")
ft.ov.dat$ft <- as.numeric(ft.ov.dat$ft)
ft.ov.dat$wue <- as.numeric(ft.ov.dat$wue)

# merge flowering time and ovule data with dat$pheno
dat.tk$pheno <- merge(x = dat.tk$pheno, y = ft.ov.dat, by = "id", all.x = TRUE, all.y = FALSE)
dat.tk$pheno <- dat.tk$pheno[, c("id", "lsm_ss_num", "mean_ss_num", "lsm_ov_num", "mean_ov_num", "ft", "wue")]

# quantile normalization of mean ss num
(normObj.ss <- orderNorm(dat.tk$pheno$mean_ss_num))
dat.tk$pheno$norm_mean_ss <- predict(normObj.ss, newdata = dat.tk$pheno$mean_ss_num)

# get relevant TK phenotypes for correlation table
dat.tk.tab <- cbind(dat.tk$pheno[c(3,6,5,7)])

# read in Royer cross data for Belm Roda

dat.si <- read.csv("./Data/belm_roda_pheno_geno/Royer final QTL data for analysis.csv")

# remove irrelevant columns
dat.si <- dat.si[,4:5]

(normObj.ss <- orderNorm(dat.si$mean_stamen_number))
dat.si$norm_mean_ss <- predict(normObj.ss, newdata = dat.si$mean_ss_num)

# read in flowering time data for Swe It
si.ft <- read_excel("./Data/belm_roda_pheno_geno/Chamber_FloweringTime.xlsx", sheet="Flowering Time_RILs")

# merge short stamen number and flowering time into one dataframe
dat.si <- merge(dat.si, si.ft, by.x="id", by.y = "RIL")

# read in SW IT ovule number
si.ov <- read_excel("./Data/belm_roda_pheno_geno/SW IT OvuleNo.xlsx")
colnames(si.ov) <- c("RIL", "ril_ovule_num", "spacer", "parent", "line", "parent_ovule_num")

dat.si <- merge(dat.si, si.ov[,1:2], by.x="id", by.y="RIL")

# remove irrelevant columns
dat.si <- dat.si[,c(1:3, 5, 7:8)]

# make days to flower a number
dat.si$IT_meandaystoflr <- as.numeric(dat.si$IT_meandaystoflr)
dat.si$SW_meandaystoflr <- as.numeric(dat.si$SW_meandaystoflr)


correlation_matrix(
  dat.si[,c(2, 4:6)],
  type = "pearson",
  digits = 3,
  decimal.mark = ".",
  use = "lower",
  show_significance = TRUE,
  replace_diagonal = TRUE,
  replacement = ""
) %>% write.csv("./Results/Belm_Roda_corr_matrix.csv")


correlation_matrix(
  dat.tk.tab,
  type = "pearson",
  digits = 3,
  decimal.mark = ".",
  use = "lower",
  show_significance = TRUE,
  replace_diagonal = TRUE,
  replacement = ""
) %>% write.csv("./Results/Tsu_Kas_corr_matrix.csv")


# Make histograms 

