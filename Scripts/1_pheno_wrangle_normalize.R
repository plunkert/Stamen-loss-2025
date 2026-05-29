# Wrangle phenotype data, normalize, and produce histograms
# Make correlation tables for phenotypes in Tsu Kas, Belm Roda RILs

rm(list=ls())

# load required packages
require(ggplot2)
require(corrtable)
require(readxl)
require(qtl)
require(dplyr)
require(bestNormalize)
require(emmeans)
require(cowplot)

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
dat.tk <- convert2riself(dat.tk)

# read in flowering time data from Conner lab, not Lovell 2013, since parents were phenotyped
# in same conditions in Conner lab
ft_pl1 <- read_excel("./Data/tsu_kas_pheno_geno/Tsu-Kas RILs.xls", sheet="Pl1 germ-bolt")
ft_pl2 <- read_excel("./Data/tsu_kas_pheno_geno/Tsu-Kas RILs.xls", sheet="Pl2 germ-bolt")
colnames(ft_pl1) <- c("order", "line_num", "plant_date", "germ_date", "num_germ", "dff", "notes")
colnames(ft_pl2) <- colnames(ft_pl1)

# compute days to flower
ft_pl1$days_to_flower <- as.numeric(ft_pl1$dff - ft_pl1$germ_date)
ft_pl2$days_to_flower <- as.numeric(ft_pl2$dff - ft_pl2$germ_date)

# put two plantings together in one big spreadsheet
ft_pl1$planting <- rep(1,nrow(ft_pl1))
ft_pl2$planting <- rep(2,nrow(ft_pl2))
ft <- rbind(ft_pl1, ft_pl2)
ft$type <- case_when(substring(ft$line_num, 1,2)=="CS" ~ "parent",
                     .default="RIL")

# test for effect of planting date on days to flowering
m_planting <- lm(data=ft, days_to_flower ~ planting)

# Residual diagnostics to see if data meet assumptions of lm
res.vals <- resid(m_planting)
pred.vals <- fitted(m_planting) # Fitted values
# Is there a pattern in the residuals?
plot(pred.vals, res.vals, main = "Residuals vs. pred.vals values", 
     las = 1, xlab = "Predicted values", ylab = "Residuals", pch = 19)
abline(h = 0)
# well there's only two values of planting so not sure how meaningful this is. 
# More spread in positive than negative residuals

# Make qqnorm plot and see if points are along diagonal
qqnorm(res.vals, pch = 19)
qqline(res.vals, col = 'red') # yes except at tails

# did date of planting affect days to flowering?
summary(m_planting)
# it didn't! effect = -0.48 days, p=0.671

# But Jeff thinks I should compute a LSM that accounts for temporal block anyway,
# so here we go
m <- lm(data=ft, days_to_flower ~ planting + line_num)

# calculate reference grid and estimated marginal means for flowering time
ft.rg <- ref_grid(m)
ft.e <- emmeans(ft.rg, "line_num")
summary(ft.e)$emmean

# make dataframe of LSMs and raw mean flowering times

ft_mean <- filter(ft, type == "RIL") %>% group_by(line_num) %>% summarize(days_to_flower = mean(days_to_flower, na.rm=TRUE))

ft_lsm <- summary(ft.e)[,1:2] %>% merge(ft_mean, by="line_num")

colnames(ft_lsm) <- c("line_num", "ft_lsm", "ft_raw")
# how correlated?
plot(ft_lsm$ft_lsm ~ ft_lsm$ft_raw)
hist(ft_lsm$ft_lsm - ft_lsm$ft_raw)
cor(x=ft_lsm$ft_lsm, y=ft_lsm$ft_raw, method = "pearson") # very! cor = 0.9999993

# merge flowering time and ovule data
ft.ov.dat <- merge(x=ft_lsm, y = ov.dat.ids, by.x="line_num", by.y="ril")

# merge flowering time and ovule data with dat$pheno
dat.tk$pheno <- merge(x = dat.tk$pheno, y = ft.ov.dat, by = "id", all=TRUE)
dat.tk$pheno <- dat.tk$pheno[, c("id", "lsm_ss_num", "mean_ss_num", "lsm_ov_num", "mean_ov_num", "ft_lsm", "ft_raw")]

# quantile normalization of mean ss num
(tk.normObj.ss <- orderNorm(dat.tk$pheno$mean_ss_num))
dat.tk$pheno$norm_mean_ss <- predict(tk.normObj.ss, newdata = dat.tk$pheno$mean_ss_num)

# get relevant TK phenotypes for correlation table
dat.tk.tab <- cbind(dat.tk$pheno[c(3, 8, 5, 7)])

# read in Royer cross data for Belm Roda

dat.si <- read.csv("./Data/belm_roda_pheno_geno/Royer final QTL data for analysis.csv")

# remove irrelevant columns
dat.si <- dat.si[,4:5]

(si.normObj.ss <- orderNorm(dat.si$mean_stamen_number))
dat.si$norm_mean_ss <- predict(si.normObj.ss, newdata = dat.si$mean_ss_num)

# read in flowering time data for Swe It
si.ft <- read_excel("./Data/belm_roda_pheno_geno/Chamber_FloweringTime.xlsx", sheet="Flowering Time_RILs")

# merge stamen number and flowering time into one dataframe
dat.si <- merge(dat.si, si.ft, by.x="id", by.y = "RIL", all=TRUE)

# read in SW IT ovule number
si.ov <- read_excel("./Data/belm_roda_pheno_geno/SW IT OvuleNo.xlsx")
colnames(si.ov) <- c("id", "ril_ovule_num", "spacer", "parent", "line", "parent_ovule_num")

dat.si <- merge(dat.si, si.ov[,1:2], by="id", all=TRUE)

# put parental phenotypes in their own columns to use later
si_parent_ov <-si.ov[,4:6] %>% na.omit()

# remove irrelevant columns
dat.si <- dat.si[,c(1:3, 5, 7:8)]

# make days to flower a number
dat.si$IT_meandaystoflr <- as.numeric(dat.si$IT_meandaystoflr)
dat.si$SW_meandaystoflr <- as.numeric(dat.si$SW_meandaystoflr)

# subtract 4 from stamen number since long stamens were invariant
dat.si$short_stamen_number <- dat.si$mean_stamen_number - 4

si_cor_table <- correlation_matrix(
  dat.si[,c(7, 6, 4)],
  type = "pearson",
  digits = 3,
  decimal.mark = ".",
  use = "lower",
  show_significance = TRUE,
  replace_diagonal = TRUE,
  replacement = ""
) 


tk_cor_table <- correlation_matrix(
  dat.tk.tab[,c(1,3,4)],
  type = "pearson",
  digits = 3,
  decimal.mark = ".",
  use = "lower",
  show_significance = TRUE,
  replace_diagonal = TRUE,
  replacement = ""
)

write.csv(rbind(tk_cor_table, rep("-", 3), si_cor_table), file="./Results/cor_table_both_RILs.csv")

# how correlated are LSMs and raw short stamen number?
cor(dat.tk$pheno$lsm_ss_num, dat.tk$pheno$mean_ss_num, use="complete.obs")
cor(dat.tk$pheno$lsm_ov_num, dat.tk$pheno$mean_ov_num, use="complete.obs")

# read in short stamen number (from Conner lab onedrive, just to see if I should 
# rerun analyses using these instead of Anne's published total stamen count)
si_ss <- read.csv("./Data/belm_roda_pheno_geno/RIL mean ShortStamenNo.csv")
si_ss_ts <- merge(x=si_ss, y=dat.si, by.x="RIL", by.y="id")
si_ss_ts$diff <- si_ss_ts$mean_stamen_number - 4 - si_ss_ts$MeanShortStamen.
View(si_ss_ts)
# the difference between short stamen number and total stamen number - 4 is 0 
# (on scale of 1e-16, which is due to rounding, not an observed difference)
# So I can proceed with total stamen number knowing that long stamens are invariant

# What are parental means stamen number phenotypes quantile normalization? 
# Sweden raw mean stamen number was 5.96, Italy mean stamen number was 4.92
predict(si.normObj.ss, newdata=c(5.96, 4.92))

# Tsu mean stamen number was 1.989, Kas mean stamen number was 0.68 in Greenhouse 2014 expt
predict(tk.normObj.ss, newdata=c(0.68, 1.989))

# Export wrangled CSVs to file for QTL mapping phenotypes

write.cross(dat.tk, format="csv", filestem="./Data/tsu_kas_pheno_geno/Tsu_Kas_Rqtl_format", chr=c("1", "2", "3", "4", "5") )
write.csv(dat.si, "./Data/belm_roda_pheno_geno/wrangled_belm_roda_phenos.csv" , row.names=FALSE)

# Statistical test for difference in parental mean phenotypes. Get means and CIs to
# overlay on histograms
ft_parents <- ft[grepl(pattern="/", x=ft$line_num),]
ft_parents$accession <- case_when(grepl("CS903", ft_parents$line_num) ~ "Tsu-1",
                                  grepl("CS1640", ft_parents$line_num) ~ "Kas-1") %>% as.factor()
t.test(data=ft_parents, days_to_flower ~ accession)
ft_tk_m <- lm(data=ft_parents, days_to_flower ~ accession)
ft_tk_emm <- emmeans(ft_tk_m, ~ accession) %>% as.data.frame()

ov_si_m <- lm(data=si_parent_ov, parent_ovule_num ~ parent)
summary(ov_si_m)
ov_si_emm <- emmeans(ov_si_m, ~ parent) %>% as.data.frame()

ov_tk_parents <- read_excel("./Data/tsu_kas_pheno_geno/TsuKas_RIL_ovule_counts.xlsx") %>%
  filter(grepl(pattern="/", plant.id))

ov_tk_parents$accession <- case_when(grepl("1640", ov_tk_parents$plant.id) ~ "Kas-1", 
                                     grepl("903", ov_tk_parents$plant.id) ~ "Tsu-1")

ov_tk_m <- lm(data=ov_tk_parents, MeanOvuleNo ~ accession)
summary(ov_tk_m)
ov_tk_emm <- emmeans(ov_tk_m, ~ accession) %>% as.data.frame()

tk_ss_raw <- read_excel("./Data/tsu_kas_pheno_geno/Tsu-Kas RILs.xls", sheet="Stamen counts")
colnames(tk_ss_raw) <- c("tube_num", "plant_id", "long_n1", "short_n1", "long_n2", "short_n2", 
                         "long_n3", "short_n3", "long_n4", "short_n4")
tk_ss_parents <- tk_ss_raw[grepl(pattern="/", x=tk_ss_raw$plant_id),]

tk_ss_parents <- tk_ss_parents %>% mutate(short_n1 = as.numeric(short_n1), short_n2 = as.numeric(short_n2),
                                          short_n3 = as.numeric(short_n3), short_n4 = as.numeric(short_n4))
tk_ss_parents$mean_short <- rowMeans(x=tk_ss_parents[,c(4,6,8,10)], na.rm=TRUE) # compute mean ss# for each tube

tk_parent_ss_means <- tk_ss_parents %>% group_by(plant_id) %>% summarize(plant_mean_ss = mean(mean_short)) # compute mean ss# each plant
tk_parent_ss_means$parent <- case_when(grepl("1640", tk_parent_ss_means$plant_id) ~ "Kas-1", 
                                       grepl("903", tk_parent_ss_means$plant_id) ~ "Tsu-1")

ss_tk_m <- lm(data=tk_parent_ss_means, plant_mean_ss ~ parent)
summary(ss_tk_m)
ss_tk_emm <- emmeans(ss_tk_m, ~ parent) %>% as.data.frame()

ss_si_dat <- read.csv("./Data/belm_roda_pheno_geno/SW_IT_RILparents_forVariance.csv")
ss_si_dat$plant <- paste(ss_si_dat$Pop, ss_si_dat$Line, ss_si_dat$Rep, sep="_")
ss_si_means <- ss_si_dat %>% group_by(plant) %>% summarize(mean=mean(ShortStamen))
ss_si_means$parent <- substr(ss_si_means$plant, 0,2)

ss_si_m <- lm(data=ss_si_means, mean ~ parent)
summary(ss_si_m)
ss_si_emm <- emmeans(ss_si_m, ~ parent) %>% as.data.frame()

# Make histograms indicating RIL phenotype distribution and parental genotypes with CIs
tk_ft_hist <- ggplot(dat.tk$pheno, aes(x = ft_raw)) +
  geom_histogram(color = "black", fill = "#81D4FA", binwidth=10) +
  labs(x = "Days to Flowering", y = "") +
  ylim(c(0,180))+
  theme_minimal()+theme(axis.title = element_text(size = 18), axis.text= element_text(size = 14))+
  annotate("pointrange", x = ft_tk_emm[2,2], y = 100, xmin = ft_tk_emm[2,5], xmax = ft_tk_emm[2,6]) +
  annotate("text", x = ft_tk_emm[2,2], y = 110, label = "Tsu-1", size=6)+
  annotate("pointrange", x = ft_tk_emm[1,2], y = 160, xmin = ft_tk_emm[1,5], xmax = ft_tk_emm[1,6])+
  annotate("text", x = ft_tk_emm[1,2], y = 170, label = "Kas-1", size=6)

tk_ov_hist <- ggplot(dat.tk$pheno, aes(x = mean_ov_num)) +
  geom_histogram(color = "black", fill = "#81D4FA", binwidth=5) +
  labs(x = "Ovule Number", y = "") +
  ylim(c(0,180))+ xlim(c(0,65))+
  theme_minimal()+theme(axis.title = element_text(size = 18), axis.text= element_text(size = 14))+
  annotate("pointrange", x = ov_tk_emm[2,2], y = 140, xmin = ov_tk_emm[2,5], xmax = ov_tk_emm[2,6]) +
  annotate("text", x = ov_tk_emm[2,2], y = 150, label = "Tsu-1", size=6)+
  annotate("pointrange", x = ov_tk_emm[1,2], y = 165, xmin = ov_tk_emm[1,5], xmax = ov_tk_emm[1,6])+
  annotate("text", x = ov_tk_emm[1,2], y = 175, label = "Kas-1", size=6)

tk_ss_hist <- ggplot(dat.tk$pheno, aes(x = mean_ss_num)) +
  geom_histogram(color = "black", fill = "#81D4FA", binwidth=0.2) +
  labs(x = "Short Stamen Number", y = "") +
  ylim(c(0,110))+
  theme_minimal()+theme(axis.title = element_text(size = 18), axis.text= element_text(size = 14))+
  annotate("pointrange", x = ss_tk_emm[2,2], y = 90, xmin = ss_tk_emm[2,5], xmax = ss_tk_emm[2,6]) +
  annotate("text", x = ss_tk_emm[2,2], y = 100, label = "Tsu-1", size=6)+
  annotate("pointrange", x = ss_tk_emm[1,2], y = 90, xmin = ss_tk_emm[1,5], xmax = ss_tk_emm[1,6])+
  annotate("text", x = ss_tk_emm[1,2], y = 100, label = "Kas-1", size=6)


# Dittmar 2014 didn't report raw data for parents, so extracted 95% confidence intervals
# from their Figure 2A plot
# 95% CI for Italy parent in Italy enviro: 102.531-105.823
# 95% CI for Sweden parent in Italy enviro: 123.250-123.885
si_ft_hist <- ggplot(dat.si, aes(x = IT_meandaystoflr)) +
  geom_histogram(color = "black", fill = "#01579B", binwidth=2) +
  labs(x = "Days to Flowering (Italy Chamber)", y = "") +
  ylim(c(0,180))+
  theme_minimal()+theme(axis.title = element_text(size = 18), axis.text= element_text(size = 14))+
  annotate("pointrange", x = 123.59, y = 140, xmin = 123.25, xmax = 123.885) +
  annotate("text", x = 124.5, y = 148, label = "Roda-47", size=6)+
  annotate("pointrange", x = 104.17, y = 140, xmin = 102.531, xmax = 105.823)+
  annotate("text", x = 104.17, y = 148, label = "Belm-12", size=6)

si_ov_hist <- ggplot(dat.si, aes(x = ril_ovule_num)) +
  geom_histogram(color = "black", fill = "#01579B", binwidth=5) +
  labs(x = "Ovule Number", y = "") +
  ylim(c(0,180))+ xlim(c(0,65))+
  theme_minimal()+theme(axis.title = element_text(size = 18), axis.text= element_text(size = 14))+
  annotate("pointrange", x = ov_si_emm[1,2], y = 152, xmin = ov_si_emm[1,5], xmax = ov_si_emm[1,6]) +
  annotate("text", x = ov_si_emm[1,2]+5, y = 162, label = "Belm-12", size=6)+
  annotate("pointrange", x = ov_si_emm[2,2], y = 152, xmin = ov_si_emm[2,5], xmax = ov_si_emm[2,6])+
  annotate("text", x = ov_si_emm[2,2]-5, y = 162, label = "Roda-47", size=6)

si_ss_hist <- ggplot(dat.si, aes(x = short_stamen_number)) +
  geom_histogram(color = "black", fill = "#01579B", binwidth=0.2) +
  labs(x = "Short Stamen Number", y = "") +
  ylim(c(0,220))+
  theme_minimal()+theme(axis.title = element_text(size = 18), axis.text= element_text(size = 14))+
  annotate("pointrange", x = ss_si_emm[1,2], y = 196, xmin = ss_si_emm[1,5], xmax = ss_si_emm[1,6]) +
  annotate("text", x = ss_si_emm[1,2], y = 209, label = "Belm-12", size=6)+
  annotate("pointrange", x = ss_si_emm[2,2], y = 196, xmin = ss_si_emm[2,5], xmax = ss_si_emm[2,6])+
  annotate("text", x = ss_si_emm[2,2], y = 209, label = "Roda-47", size=6)

# Make multipanel figure with all the histograms!!!

histograms <- plot_grid(tk_ss_hist, si_ss_hist, tk_ov_hist, si_ov_hist, ncol = 2,
                        tk_ft_hist, si_ft_hist, labels = "AUTO", label_size=18, align="hv")

save_plot("./Results/all_histograms.svg", plot=histograms, base_width = 8, base_height = 10)

# Make supplemental figure with just quantile normalized short stamen number histograms
tk_qn_ss_hist <- ggplot(dat.tk$pheno, aes(x = norm_mean_ss)) +
  geom_histogram(color = "black", fill = "#81D4FA", binwidth=0.2) +
  labs(x = "Quantile-Normalized Short Stamen Number", y = "") +
  ylim(c(0,60))+
  theme_minimal()+theme(axis.title = element_text(size = 18), axis.text= element_text(size = 12))+
  annotate("point", x = -1.671206, y = 40) +
  annotate("text", x = -1.671206, y =  50, label = "Kas-1", size=6)+
  annotate("point", x = 1.720862, y = 40)+
  annotate("text", x = 1.720862, y = 50, label = "Tsu-1", size=6)

si_qn_ss_hist <- ggplot(dat.si, aes(x = norm_mean_ss)) +
  geom_histogram(color = "black", fill = "#01579B", binwidth=0.2) +
  labs(x = "Quantile-Normalized Short Stamen Number", y = "") +
  ylim(c(0,175))+
  theme_minimal()+theme(axis.title = element_text(size = 18), axis.text= element_text(size = 12))+
  annotate("point", x = -1.8317004, y = 140) +
  annotate("text", x = -1.8317004, y = 150, label = "Belm-12", size=6)+
  annotate("point", x = 0.6579079, y = 140)+
  annotate("text", x = 0.6579079, y = 150, label = "Roda-47", size=6)

QN_hist <- plot_grid(tk_qn_ss_hist, si_qn_ss_hist, labels = "AUTO")
save_plot("./Results/QN_ss_histograms.svg", plot=QN_hist, base_width = 8, base_height = 5)

# let's check whether marker order on chr5 matches the Col-0 reference genome physical positions

tk.map <- read.csv("./Data/tsu_kas_pheno_geno/Tsu_Kas_Rqtl_format.csv", header=FALSE)[1:3,-c(1:9)] %>% t() %>% as.data.frame()
si.map <- read.csv("./Data/belm_roda_pheno_geno/IT SW RIL genotypes.csv", header=FALSE)[1:3,] %>% t() %>% as.data.frame()

colnames(tk.map) <- c("marker", "chr", "cM")
colnames(si.map) <- c("marker", "chr", "cM")

# if marker name includes a C in Tsu Kas or an X in Belm Roda, remove it
tk.map$marker <- gsub("C","",tk.map$marker)
si.map$marker <- gsub("X","",si.map$marker)

# Now most of the markers are named as chromNum-bpPosition or chromNum_bpPosition
# The markers whose names don't match that will be coerced to NA
tk.map$bp <- substr(tk.map$marker, 3, nchar(tk.map$marker)) %>% as.numeric()
si.map$bp <- substr(si.map$marker, 3, nchar(si.map$marker)) %>% as.numeric()


# Subset to get chromosomes of interest
tk.map.chr5 <- tk.map[which(tk.map$chr==5),]
si.map.chr5 <- si.map[which(si.map$chr==5),]
tk.map.chr1 <- tk.map[which(tk.map$chr==1),]
si.map.chr1 <- si.map[which(si.map$chr==1),]

# Let's plot genetic map against Col-0 physical map for chromosomes 1 and 5

svg("./Results/bp_vs_cM_chr1_chr5.svg", width=12, height=12, pointsize=18)

par(mfrow=c(2,2))


plot(x=tk.map.chr1$cM, y=tk.map.chr1$bp, xlab="Tsu-1 x Kas-1 Genetic Map (cM)",
     ylab="Col-0 Physical Map (bp)")
abline(v=80, col="red")
abline(v=108, col="red")
plot(x=si.map.chr1$cM, y=si.map.chr1$bp, xlab="Belm-12 x Roda-47 Genetic Map (cM)",
     ylab="Col-0 Physical Map (bp)")
abline(v=74.4, col="red")

# plot chr5
plot(x=tk.map.chr5$cM, y=tk.map.chr5$bp, xlab="Tsu-1 x Kas-1 Genetic Map (cM)",
     ylab="Col-0 Physical Map (bp)")
abline(v=2.8, col="red")
abline(v=9.8, col="red")
plot(x=si.map.chr5$cM, y=si.map.chr5$bp, xlab="Belm-12 x Roda-47 Genetic Map (cM)",
     ylab="Col-0 Physical Map (bp)")
abline(v=7.7, col="red")
#abline(v=13.6, col="red")

dev.off()
