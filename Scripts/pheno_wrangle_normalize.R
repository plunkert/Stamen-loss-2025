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
(tk.normObj.ss <- orderNorm(dat.tk$pheno$mean_ss_num))
dat.tk$pheno$norm_mean_ss <- predict(tk.normObj.ss, newdata = dat.tk$pheno$mean_ss_num)

# get relevant TK phenotypes for correlation table
dat.tk.tab <- cbind(dat.tk$pheno[c(3,6,5,7)])

# read in Royer cross data for Belm Roda

dat.si <- read.csv("./Data/belm_roda_pheno_geno/Royer final QTL data for analysis.csv")

# remove irrelevant columns
dat.si <- dat.si[,4:5]

(si.normObj.ss <- orderNorm(dat.si$mean_stamen_number))
dat.si$norm_mean_ss <- predict(si.normObj.ss, newdata = dat.si$mean_ss_num)

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

# What are parental means stamen number phenotypes quantile normalization? 
# Sweden raw mean stamen number was 5.96, Italy mean stamen number was 4.92
predict(si.normObj.ss, newdata=c(5.96, 4.92))

# Tsu mean stamen number was 1.989, Kas mean stamen number was 0.68 in Greenhouse 2014 expt
predict(tk.normObj.ss, newdata=c(0.68, 1.989))


# Export wrangled CSVs to file for QTL mapping phenotypes

write.cross(dat.tk, format="csv", filestem="./Data/tsu_kas_pheno_geno/Tsu_Kas_Rqtl_format", chr=c("1", "2", "3", "4", "5") )
write.csv(dat.si, "./Data/belm_roda_pheno_geno/wrangled_belm_roda_phenos.csv" , row.names=FALSE)

# Make histograms that indicate parental genotypes for Belm Roda

svg("./Results/Belm_Roda_pheno_dist.svg", width=12, height=7)
par(mfrow=c(2,3))
hist(dat.si$mean_stamen_number, xlab="Mean Stamen Number", main="", ylim=c(0,220), cex.lab=2)
arrows(5.96, 205, 5.96, 190) # sweden parents from Royer et al.
text(x = 5.96, y= 215, "Roda", cex=2)
arrows(4.92, 60, 4.92, 25) # italy parents from Royer et al.
text(x = 4.92, y= 72, "Belm", cex=2) # sweden parents from Royer et al. 

hist(dat.si$norm_mean_ss, xlab="Quantile-Normalized Mean Stamen Number", main="", cex.lab=2)
arrows(0.6623889, 110, 0.6623889, 95) # sweden parents from Royer et al.
text(x = 0.6623889, y= 115, "Roda", cex=2)
arrows(-1.8308334, 70, -1.8308334, 55) # italy parents from Royer et al.
text(x = -1.8308334, y= 75, "Belm", cex=2) # sweden parents from Royer et al.

hist(dat.si$IT_meandaystoflr, xlab="Days to Flower in Italy Chamber", main="", cex.lab=2)
arrows(123.59, 140, 123.59, 120) # sweden parents from Dittmar et al
text(x = 123.59, y= 147, "Roda", cex=2)
arrows(104.17, 140, 104.17, 120) # italy parents from Dittmar et al
text(x = 104.17, y= 147, "Belm", cex=2)

hist(dat.si$SW_meandaystoflr, xlab="Days to Flower in Sweden Chamber", main="", ylim=c(0,60), cex.lab=2)
arrows(138.08, 55, 138.08, 45) # sweden parents
text(x = 138.08, y = 58, "Roda", cex=2)
arrows(128.29, 55, 128.29, 49) # italy parents
text(x = 128.29, y= 58, "Belm", cex=2)

hist(dat.si$ril_ovule_num, xlab="Mean Ovule Number", main="", ylim=c(0,175), cex.lab=2)
arrows(34.25, 165, 34.25, 150) # sweden parents
text(x = 34.25, y = 172, "Roda", cex=2)
arrows(40.2, 145, 40.2, 130) # italy parents
text(x = 40.2, y= 152, "Belm", cex=2)

dev.off()

# Compute parental mean ovule number phenotypes
tk.ov.parent <- read_excel("./Data/tsu_kas_pheno_geno/TsuKas_parent_ovule_counts_MLPentry.xlsx")
View(tk.ov.parent)

k.ov.parent <- filter(tk.ov.parent, par.line=="Kas") %>% 
  group_by(pot.num) %>% summarize(mean_ovule = mean(ovule.num))
mean(k.ov.parent$mean_ovule) # Kas parental line has mean ovule number of 55.11

t.ov.parent <- filter(tk.ov.parent, par.line=="Tsu") %>% 
  group_by(pot.num) %>% summarize(mean_ovule = mean(ovule.num))
mean(t.ov.parent$mean_ovule) # Tsu parental line has mean ovule number of 53.04


# Make histograms that indicate parental genotypes for Tsu Kas

svg("./Results/Tsu_Kas_pheno_dist.svg", width=12, height=7)
par(mfrow=c(2,3))
hist(dat.tk$pheno$mean_ss_num, xlab="Mean Short Stamen Number", main="", ylim=c(0,220), cex.lab=2)
arrows(0.68, 80, 0.68, 50) 
text(x = 0.68, y= 95, "Kas", cex=2)
arrows(1.989, 130, 1.989, 100)
text(x = 1.989, y= 145, "Tsu", cex=2)

hist(dat.tk$pheno$norm_mean_ss, xlab="Quantile-Normalized Mean Stamen Number", main="", cex.lab=2)
arrows(-1.671206, 45, -1.671206, 30)
text(x = -1.671206, y=50, "Kas", cex=2)
arrows(1.720862, 45, 1.720862, 30)
text(x = 1.720862, y= 50, "Tsu", cex=2)

# Extracted flowering time and WUE parent data for Tsu Kas from Lovell et al 2013 
# Supp Fig S1 using plotdigitizer.com

# Tsu mean flowering time = 11.332 days
# Kas mean flowering time = 16.478
# Tsu mean WUE = -31.369
# Kas mean WUE = -29.883

# This is not really plausible as parental phenotypes for this set of RILs? I think
# these parents were phenotyped in a different common garden from the RILs. Leaving
# out parental means.

hist(dat.tk$pheno$ft, xlab="Days to Flower", main="", cex.lab=2)

hist(dat.tk$pheno$wue, xlab="Water Use Efficiency (% carbon-13)", main="", cex.lab=2)

hist(dat.tk$pheno$mean_ov_num, xlab="Mean Ovule Number", main="", ylim=c(0,175), cex.lab=2)
arrows(55.11, 120, 55.11, 100)
text(x = 56, y = 130, "Kas", cex=2)
arrows(53.04, 145, 53.04, 130)
text(x = 53.04, y= 155, "Tsu", cex=2)

dev.off()




