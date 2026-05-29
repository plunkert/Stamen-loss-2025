# Fit multi-QTL models

library(qtl)
library(dplyr)
library(ggplot2)
library(RCurl)

setwd("~/Documents/GitHub/Stamen-loss-2025/")

# Read in crosses. This part is repetitive with scanone scantwo script
dat.si <- read.cross("csvs", genotypes=c("a","b"),
                     genfile="./Data/belm_roda_pheno_geno/IT SW RIL genotypes.csv",
                     phefile= "./Data/belm_roda_pheno_geno/wrangled_belm_roda_phenos.csv",
                     na.strings=c("NA","-"))

dat.tk <-read.cross("csv", "./Data/tsu_kas_pheno_geno/", file = "Tsu_Kas_Rqtl_format.csv",
                    genotypes=c("AA", "BB"), na.strings=c("NA","-"))

# tell Rqtl it's RILs from a two-parent cross
class(dat.tk)[1] <-"riself"
class(dat.si)[1] <-"riself"

# jitter for overlapping markers
dat.si.jitter<-jittermap(dat.si)
dat.tk.jitter<-jittermap(dat.tk)

# imputation of genotype data
dat.si.calc <- calc.genoprob(dat.si.jitter, step=1, 
                             map.function="kosambi", stepwidth = "fixed")
dat.tk.calc <- calc.genoprob(dat.tk.jitter, step=1, 
                             map.function="kosambi", stepwidth = "fixed")


# Read in scantwo permutations

load("./Results/si_mean_stamen_scantwo_perm.RData")
load("./Results/si_norm_stamen_scantwo_perm.RData")
load("./Results/si_ovule_scantwo_perm.RData")

load("./Results/tk_mean_short_stamen_scantwo_perm.RData")
load("./Results/tk_norm_short_stamen_scantwo_perm.RData")
load("./Results/tk_ovule_scantwo_perm.RData")

# get penalties at 0.05 and 0.01 significance levels
pen.si.scantwo.0.01 <- calc.penalties(out.si.scan2.perm, alpha = 0.01)
pen.si.norm.scantwo.0.01 <- calc.penalties(out.si.norm.scan2.perm, alpha = 0.01)
pen.si.ov.scantwo.0.01 <- calc.penalties(out.si.ov.scan2.perm, alpha = 0.01)

pen.tk.scantwo.0.01 <- calc.penalties(out.tk.scan2.perm, alpha = 0.01)
pen.tk.norm.scantwo.0.01 <- calc.penalties(out.tk.norm.scan2.perm, alpha = 0.01)
pen.tk.ov.scantwo.0.01 <- calc.penalties(out.tk.ov.scan2.perm, alpha = 0.01)

pen.si.scantwo.0.05 <- calc.penalties(out.si.scan2.perm, alpha = 0.05)
pen.si.norm.scantwo.0.05 <- calc.penalties(out.si.norm.scan2.perm, alpha = 0.05)
pen.si.ov.scantwo.0.05 <- calc.penalties(out.si.ov.scan2.perm, alpha = 0.05)

pen.tk.scantwo.0.05 <- calc.penalties(out.tk.scan2.perm, alpha = 0.05)
pen.tk.norm.scantwo.0.05 <- calc.penalties(out.tk.norm.scan2.perm, alpha = 0.05)
pen.tk.ov.scantwo.0.05 <- calc.penalties(out.tk.ov.scan2.perm, alpha = 0.05)

# Fit stepwise multi-QTL models

step.out.si.01 <-stepwiseqtl(dat.si.calc, pheno.col="mean_stamen_number",
                             model="normal", method="hk", penalties=pen.si.scantwo.0.01, max.qtl=16,
                             verbose=F, refine.locations=T,additive.only=F, keeplodprofile=T,keeptrace=T)

step.out.norm.si.01 <-stepwiseqtl(dat.si.calc, pheno.col="norm_mean_ss",
                                  model="normal", method="hk", penalties=pen.si.norm.scantwo.0.01, max.qtl=16,
                                  verbose=F, refine.locations=T,additive.only=F, keeplodprofile=T,keeptrace=T)


step.out.ov.si.01 <- stepwiseqtl(dat.si.calc, pheno.col="ril_ovule_num",
                                 model="normal", method="hk", penalties=pen.si.ov.scantwo.0.01, max.qtl=16,
                                 verbose=F, refine.locations=T,additive.only=F, keeplodprofile=T,keeptrace=T)

step.out.tk.01 <-stepwiseqtl(dat.tk.calc, pheno.col="mean_ss_num",
                                  model="normal", method="hk", penalties=pen.tk.scantwo.0.01, max.qtl=16,
                                  verbose=F, refine.locations=T,additive.only=F, keeplodprofile=T,keeptrace=T)

step.out.norm.tk.01 <-stepwiseqtl(dat.tk.calc, pheno.col="norm_mean_ss",
                             model="normal", method="hk", penalties=pen.tk.norm.scantwo.0.01, max.qtl=16,
                             verbose=F, refine.locations=T,additive.only=F, keeplodprofile=T,keeptrace=T)


step.out.ov.tk.01 <- stepwiseqtl(dat.tk.calc, pheno.col="mean_ov_num",
                                 model="normal", method="hk", penalties=pen.tk.ov.scantwo.0.01, max.qtl=16,
                                 verbose=F, refine.locations=T,additive.only=F, keeplodprofile=T,keeptrace=T)



# Explore variation in max.qtl parameter with significance threshold = 0.05
pdf("./Results/vary_max_qtl_tk_mean_ss_0.01.pdf")
par(mfrow = c(4, 3))
for (i in 5:16){
  sw <- stepwiseqtl(dat.tk.calc, pheno.col="mean_ss_num", max.qtl=i,
                    method="hk", penalties = pen.tk.scantwo.0.01, keeplodprofile=TRUE, keeptrace=TRUE)
  plotLodProfile(sw, col=c("black", "red"))
  title(paste("alpha = 0.01, max.qtl = ", i))
}
dev.off()

# With significance threshold = 0.05
pdf("./Results/vary_max_qtl_tk_mean_ss_05.pdf")
par(mfrow = c(4, 3))
for (i in 5:16){
  sw <- stepwiseqtl(dat.tk.calc, pheno.col="mean_ss_num", max.qtl=i,
                    method="hk", penalties = pen.tk.scantwo.0.05, keeplodprofile=TRUE, keeptrace=TRUE)
  plotLodProfile(sw, col=c("black", "red"))
  title(paste("alpha = 0.05, max.qtl = ", i))
}
dev.off()

# At alpha=0.05, the QTL identified keep changing as max.qtl is increased, whereas 
# at alpha=0.05, the QTL identified stabilize once max.qtl=6. Since multi-QTL model 
# at alpha=0.05 is robust to changes in max.qtl, I'll use this significance threshold
# moving forward.

#81D4FA is sweden italy color
#01579B is tsu kas color

# Plot LOD profiles with significance threshold
svg("./Results/stepwise_lod_ss.svg", width=8, height=10)
par(mfrow=c(2,1))
plotLodProfile(step.out.tk.01, showallchr=TRUE, ylim=c(0,30), ylab="LOD", qtl.labels=FALSE,
               col=c('black', 'darkgray', 'black', 'black', 'darkgray', 'black'))
abline(h=pen.tk.scantwo.0.01[1], col="firebrick3", lty='dashed')

plotLodProfile(step.out.si.01, showallchr=TRUE, ylim=c(0,55), ylab="LOD", qtl.labels=FALSE)
abline(h=pen.si.scantwo.0.01[1], col="firebrick3", lty='dashed')

dev.off()

svg("./Results/stepwise_lod_qn_ss.svg", width=8, height=10)
par(mfrow=c(2,1))

plotLodProfile(step.out.norm.tk.01, showallchr=TRUE, ylim=c(0,40), ylab="LOD", qtl.labels=FALSE,
               col=c('darkgray', 'black', 'darkblue', 'black', 'black', 'black'))
abline(h=pen.tk.norm.scantwo.0.01[1], col="firebrick3", lty='dashed')

plotLodProfile(step.out.norm.si.01, showallchr=TRUE,, ylim=c(0,55), ylab="LOD", qtl.labels=FALSE)
abline(h=pen.si.norm.scantwo.0.01[1], col="firebrick3", lty='dashed')

dev.off()

svg("./Results/stepwise_lod_ovule.svg", width=8, height=10)
par(mfrow=c(2,1))
plotLodProfile(step.out.ov.tk.01, showallchr=TRUE, ylab="LOD", qtl.labels=FALSE,
               col=c('black', 'black', 'darkgray', 'black', 'black'))
abline(h=pen.tk.ov.scantwo.0.01[1], col="firebrick3", lty='dashed')

plotLodProfile(step.out.ov.si.01, showallchr=TRUE, ylab="LOD", qtl.labels=FALSE)
abline(h=pen.si.ov.scantwo.0.01[1], col="firebrick3", lty='dashed')
dev.off()

# What % variance explained by step.out.si.01 and step.out.tk.01?
fit.si <- fitqtl(dat.si.calc,
            pheno.col="mean_stamen_number",
            qtl=step.out.si.01,
            formula=formula(step.out.si.01),
            get.ests=T,dropone=T,
            method="hk")

summary(fit.si)

fit.tk <- fitqtl(dat.tk.calc,
                 pheno.col="mean_ss_num",
                 qtl=step.out.tk.01,
                 formula=formula(step.out.tk.01),
                 get.ests=T,dropone=T,
                 method="hk")
summary(fit.tk)
# Using John Lovell's stepwiseStats function to extract stats from stepwise model

u <- "https://raw.githubusercontent.com/jtlovell/r-QTL_functions/master/stepwiseStats.R"
script <- getURL(u, ssl.verifypeer = FALSE)
eval(parse(text=script))

si_table <- stepwiseStats(cross=dat.si.calc, model.in=step.out.si.01, phe="mean_stamen_number")

si_norm_table <- stepwiseStats(cross=dat.si.calc, model.in=step.out.norm.si.01, phe="norm_mean_ss")

si_ov_table <- stepwiseStats(cross=dat.si.calc, model.in=step.out.ov.si.01, phe="ril_ovule_num")

tk_table <- stepwiseStats(cross=dat.tk.calc, model.in=step.out.tk.01, phe="mean_ss_num")

tk_norm_table <- stepwiseStats(cross=dat.tk.calc, model.in=step.out.norm.tk.01, phe="norm_mean_ss")

tk_ov_table <- stepwiseStats(cross=dat.tk.calc, model.in=step.out.ov.tk.01, phe="mean_ov_num")


# Add percent of parental divergence explained as 100 * 2*effect.estimate/(Sweden - Italy)
si_table$perc_div <- 100*2*si_table$effect.estimate/(5.96-4.92)

si_norm_table$perc_div <- 100*2*si_norm_table$effect.estimate/(0.6623889 + 1.8308334)

# Add percent of parental divergence explained as 100 * 2*effect.estimate/(Tsu - Kas)
tk_table$perc_div <- 100*2*tk_table$effect.estimate/(1.989-0.68)

tk_norm_table$perc_div <- 100*2*tk_norm_table$effect.estimate/(1.720862+1.671206)

# Export tables to CSV
si_table[,c(1,3,4,17,18,7,12,8,19)] %>% mutate_if(is.numeric, round, digits = 2) %>%
  write.csv("./Results/si_mean_stamen_QTL_table.csv", row.names=FALSE)

tk_table[,c(1,3,4,17,18,7,12,8,19)] %>% mutate_if(is.numeric, round, digits = 2) %>%
  write.csv("./Results/tk_mean_stamen_QTL_table.csv", row.names=FALSE)

si_norm_table[,c(1,3,4,17,18,7,12,8,19)] %>% mutate_if(is.numeric, round, digits = 2) %>%
  write.csv("./Results/si_norm_stamen_QTL_table.csv", row.names=FALSE)

tk_norm_table[,c(1,3,4,17,18,7,12,8,19)] %>% mutate_if(is.numeric, round, digits = 2) %>%
  write.csv("./Results/tk_norm_stamen_QTL_table.csv", row.names=FALSE)


si_ov_table[,c(1,3,4,17,18,7,12,8)] %>% mutate_if(is.numeric, round, digits = 2) %>%
  write.csv("./Results/si_mean_ovule_QTL_table.csv", row.names=FALSE)

tk_ov_table[,c(1,3,4,17,18,7,12,8)] %>% mutate_if(is.numeric, round, digits = 2) %>%
  write.csv("./Results/tk_mean_ovule_QTL_table.csv", row.names=FALSE)


# Investigate effects of potential chr 5 linked QTL on stamen loss

# now check Sweden Italy

# need to check on local LOD peak that corresponds to shoulder of chr5 QTL
plotLodProfile(step.out.si.16qtl.01, chr="5", ylim=c(0,60), qtl.labels = FALSE, 
               show.marker.names=TRUE, main="Belm-12 x Roda-47 RILs")
# marker X5_5889105 looks centered on the shoulder
# also looked at actual LOD values for each marker, didn't just do it visually

si_pheno <- read.csv("./Data/belm_roda_pheno_geno/Royer final QTL data for analysis.csv")
si_geno <- read.csv("./Data/belm_roda_pheno_geno/IT SW RIL genotypes.csv")
si_geno <- si_geno[-c(1,2),]
si_geno$id <- as.numeric(si_geno$id)
si <- merge(x=si_pheno, y=si_geno, by="id")

si[si=="-"] <- NA

find.marker(dat.si.calc, chr="5", pos=7.652009)

# model containing only the genotype at the peak of 5@7.7
m_si_peak <- lm(data=si, mean_stamen_number ~ X5_2985496)
summary(m_si_peak) # significant at p=0.0457, b geno increases pheno by 0.21712 short stamens

# extract marker at shoulder of QTL5 trace
lodprof <- attr(step.out.si.16qtl.01, "lodprofile")
plot(as.numeric(lodprof$`5@7.7`$lod) ~ as.numeric(lodprof$`5@7.7`$pos))
max(lodprof$`5@7.7`[which(lodprof$`5@7.7`$pos>12),"lod"])

# looked for a local max lod, which occurs at 13.6 cM, marker X5_5889105
plot(as.numeric(lodprof$`5@7.7`$lod) ~ as.numeric(lodprof$`5@7.7`$pos))
max(lodprof$`5@7.7`[which(lodprof$`5@7.7`$pos>12),"lod"])

# include interaction of genotype at the peak of 5@7.7 and genotype at shoulder
si_int <- lm(data=si, mean_stamen_number ~ X5_2985496 * X5_5889105)
summary(si_int)
si_add <- lm(data=si, mean_stamen_number ~ X5_2985496 + X5_5889105)
summary(si_add)
AIC(si_add, si_int)
# like Tsu Kas, the two "loci" each additively contribute to short stamen number, 
# but not their interaction

# does additive model of two chr5 QTL explain data better than just detected QTL 5@7.7
si_qtl7.7 <- lm(data=si, mean_stamen_number ~ X5_2985496)
AIC(si_add, si_qtl7.7) # yes it does!
nobs(si_add)
nobs(si_qtl7.7) 
# there's only one observation that's in si_qtl7.7 but not si_add, so I'm not 
# concerned about the AIC warning that models were fitted to different numbers of observations

# create table describing models for MS
int_row <- cbind(parameter = c("5@7.7", "5@13.6", "Interaction"), 
                 effect = si_int$coefficients[2:4], 
                 pval = summary(si_int)$coefficients[2:4,4])

add_row <- cbind(parameter = c("5@7.7", "5@13.6"),
                 effect = si_add$coefficients[2:3],
                 pval = summary(si_add)$coefficients[2:3,4])

single_qtl_row <- cbind(parameter=c("5@7.7"),
                        effect=si_qtl7.7$coefficients[2],
                        pval = summary(si_add)$coefficients[2,4])

all_mods <- rbind(int_row, add_row, single_qtl_row)
aic_column <- c(AIC(si_int), NA, NA, AIC(si_add), NA, AIC(si_qtl7.7))
model_names <- c("5@7.7 + 5@13.6 + 5@7.7:5@13.6", NA, NA, "5@7.7 + 5@13.6", NA, "5@7.7")
write.csv(row.names=FALSE, cbind(model_names, all_mods, aic_column), "./Results/si_chr5_model_comparisons.csv")

# check with Tukey's
si_chr5 <- si[, c("mean_stamen_number","X5_2985496", "X5_5889105")] %>% na.omit()

si_chr5$geno_both_markers <- as.factor(paste(si_chr5$X5_2985496, si_chr5$X5_5889105, sep="_"))

TukeyHSD(aov(data=si_chr5, mean_stamen_number ~ geno_both_markers))
# so we don't need 2 loci to explain tukey test results, but that is quite conservative

# plot recombinant genotypes at potentially linked QTL on chr 1, chr 5 in Tsu Kas
plotPXG(dat.tk.calc, pheno.col="mean_ss_num", marker=c("5-1213621","5_2799359"), ylim=c(0,2))

set.seed(3)

plotPXG(dat.tk.calc, marker=find.marker(dat.tk.calc, chr="5", pos=c(2.8,9.8)), pheno.col="mean_ss_num", 
        main="Tsu-1 x Kas-1 RILs", ylab="mean short stamen number", ylim=c(0,2.5), infer=FALSE)

