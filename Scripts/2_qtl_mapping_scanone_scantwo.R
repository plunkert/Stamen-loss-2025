# This script maps short stamen number in Tsu Kas and Sweden Italy RILs together
# with comparable analysis parameters and plots the results.

library(qtl)
library(dplyr)
library(snow)

setwd("~/Documents/GitHub/Stamen-loss-2025/") # change this

# load cross data into R
dat.si <- read.cross("csvs", genotypes=c("a","b"),
                     genfile="./Data/belm_roda_pheno_geno/IT SW RIL genotypes.csv",
                     phefile= "./Data/belm_roda_pheno_geno/wrangled_belm_roda_phenos.csv",
                     na.strings=c("NA","-"))

dat.tk <-read.cross("csv", "./Data/tsu_kas_pheno_geno/", file = "Tsu_Kas_Rqtl_format.csv",
                 genotypes=c("AA", "AB", "BB"), na.strings=c("NA","-"))

# tell Rqtl it's RILs from a two-parent cross
class(dat.tk)[1] <-"riself"
class(dat.si)[1] <-"riself"

plot(dat.si, alternate.chrid = TRUE, pheno = "mean_stamen_number")
plot(dat.tk, alternate.chrid = TRUE, pheno = "mean_ss_num")

# jitter for overlapping markers
dat.si.jitter<-jittermap(dat.si)
dat.tk.jitter<-jittermap(dat.tk)

# imputation of genotype data
dat.si.calc <- calc.genoprob(dat.si.jitter, step=1, 
                             map.function="kosambi", stepwidth = "fixed")
dat.tk.calc <- calc.genoprob(dat.tk.jitter, step=1, 
                             map.function="kosambi", stepwidth = "fixed")

# one-dimensional QTl analysis with scanone
out.si.scan1 <- scanone(dat.si.calc, method = "hk", pheno.col= c("mean_stamen_number", "norm_mean_ss", "ril_ovule_num"))

out.tk.scan1 <- scanone(dat.tk.calc, method="hk", pheno.col=c("mean_ss_num", "norm_mean_ss", "mean_ov_num"))


#scanone permutations
out.si.scan1.perm <- scanone(dat.si.calc, method="hk",n.perm=10000, 
                             pheno.col=c("mean_stamen_number", "norm_mean_ss", "ril_ovule_num"))
out.tk.scan1.perm <- scanone(dat.tk.calc, method="hk", n.perm=10000,
                             pheno.col=c("mean_ss_num", "norm_mean_ss", "mean_ov_num"))

# save vectors of permutation thresholds. Stamen number, then normalized stamen number, 
# then ovule number
scan1.si.penalties0.05 <- summary(out.si.scan1.perm, alpha = 0.05) # lod ovule 2.64
scan1.tk.penalties0.05 <- summary(out.tk.scan1.perm, alpha = 0.05)

scan1.si.penalties0.01 <- summary(out.si.scan1.perm, alpha = 0.01) # lod ovule 3.41
scan1.tk.penalties0.01 <- summary(out.tk.scan1.perm, alpha = 0.01)



# plot scanone!
plot(out.si.scan1, lodcolumn = 1)
abline(h=scan1.si.penalties0.01[1])
plot(out.si.scan1, lodcolumn = 2)
abline(h=scan1.si.penalties0.01[2])
plot(out.si.scan1, lodcolumn = 3)
abline(h=scan1.si.penalties0.01[3])

plot(out.tk.scan1, lodcolumn = 1)
abline(h=scan1.tk.penalties0.01[1])
plot(out.tk.scan1, lodcolumn = 2)
abline(h=scan1.tk.penalties0.01[2])
plot(out.tk.scan1, lodcolumn = 3)
abline(h=scan1.tk.penalties0.01[3])

# two-dimensional scantwo analysis

out.si.scan2.perm <- scantwopermhk(dat.si.calc,
        pheno.col="mean_stamen_number", n.perm=10000,
        addcovar=NULL, weights=NULL,
        assumeCondIndep=FALSE, verbose=TRUE, batchsize = 100)

save(out.si.scan2.perm, file="./Results/si_mean_stamen_scantwo_perm.R")

out.si.norm.scan2.perm <- scantwopermhk(dat.si.calc,
                                   pheno.col="norm_mean_ss", n.perm=10000,
                                   addcovar=NULL, weights=NULL,
                                   assumeCondIndep=FALSE, verbose=TRUE, batchsize = 100)

save(out.si.norm.scan2.perm, file="./Results/si_norm_stamen_scantwo_perm.R")

out.si.ov.scan2.perm <- scantwopermhk(dat.si.calc,
                                        pheno.col="ril_ovule_num", n.perm=10000,
                                        addcovar=NULL, weights=NULL,
                                        assumeCondIndep=FALSE, verbose=TRUE, batchsize = 100)

save(out.si.ov.scan2.perm, file="./Results/si_ovule_stamen_scantwo_perm.R")

out.tk.scan2.perm <- scantwopermhk(dat.tk.calc,
                                      pheno.col="mean_ss_num", n.perm=10000,
                                      addcovar=NULL, weights=NULL,
                                      assumeCondIndep=FALSE, verbose=TRUE, batchsize = 100)

save(out.tk.scan2.perm, file="./Results/tk_mean_short_stamen_scantwo_perm.R")

out.tk.norm.scan2.perm <- scantwopermhk(dat.tk.calc,
                                      pheno.col="norm_mean_ss", n.perm=10000,
                                      addcovar=NULL, weights=NULL,
                                      assumeCondIndep=FALSE, verbose=TRUE, batchsize = 100)

save(out.tk.norm.scan2.perm, file="./Results/tk_norm_short_stamen_scantwo_perm.R")

out.tk.ov.scan2.perm <- scantwopermhk(dat.tk.calc,
                                      pheno.col="mean_ov_num", n.perm=10000,
                                      addcovar=NULL, weights=NULL,
                                      assumeCondIndep=FALSE, verbose=TRUE, batchsize = 100)

save(out.tk.ov.scan2.perm, file="./Results/tk_ovule_scantwo_perm.R")



