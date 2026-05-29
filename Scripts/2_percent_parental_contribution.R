# This script approximate percent contribution of each parent in the RILs and tests
# whether percent parental contribution affects phenotypes.

require(ppcor)

setwd("~/Documents/Github/Stamen-loss-2025/") # Change this

# read in phenotypes and genotypes
tk <- read.csv("./Data/tsu_kas_pheno_geno/Tsu_Kas_Rqtl_format.csv")[c(-1,-2),]
tk.geno <- tk[,c(1, 9:ncol(tk))]
tk.pheno <- tk[,1:8]
tk.pheno <- sapply(tk.pheno, as.numeric)

si.geno <- read.csv("./Data/belm_roda_pheno_geno/IT SW RIL genotypes.csv")
si.geno <- si.geno[3:nrow(si.geno),]
si.pheno <- read.csv("./Data/belm_roda_pheno_geno/wrangled_belm_roda_phenos.csv")

# calculate percent parental genome

# define function that returns number of occurrences in a vector matching a string
num_occurrences <- function(search_string, vector){
  length(which(vector==search_string))
}

# calculate number of occurrences of "a" genotype and "b" genotype
si.geno$num_a <- "a" %>% apply(X=si.geno, MARGIN=1, FUN=num_occurrences)
si.geno$num_b <- "b" %>% apply(X=si.geno, MARGIN=1, FUN=num_occurrences)

tk.geno$num_AA <- "AA" %>% apply(X=tk.geno, MARGIN=1, FUN=num_occurrences)
tk.geno$num_BB <- "BB" %>% apply(X=tk.geno, MARGIN=1, FUN=num_occurrences)


# calculate sum of known genotypes, excluding missing genotypes
si.geno$known_genos <- si.geno$num_a + si.geno$num_b
tk.geno$known_genos <- tk.geno$num_AA + tk.geno$num_BB

# calculate proportion each parent
si.geno$prop_swe <- si.geno$num_b / si.geno$known_genos
si.geno$prop_itl <- si.geno$num_a / si.geno$known_genos

tk.geno$prop_tsu <- tk.geno$num_BB / tk.geno$known_genos
tk.geno$prop_kas <- tk.geno$num_AA / tk.geno$known_genos

# view distribution of parental genome contribution
hist(si.geno$prop_swe)
hist(si.geno$prop_itl)
hist(tk.geno$prop_kas)
hist(tk.geno$prop_tsu)

svg("./Results/parental_genome_contribution_dist.svg", width=12, height=7)
par(mfrow=c(1,2))

hist(tk.geno$prop_tsu, main="",
     xlab="Proportion Tsu Genotype in Tsu x Kas RILs", cex.lab=1.6, xlim=c(0,1), col="#81D4FA")

hist(si.geno$prop_swe, main="",
     xlab="Proportion Roda Genotype in Belm x Roda RILs", cex.lab=1.6, xlim=c(0,1), col="#01579B")

dev.off()

# range of parental genome contribution from stamen retaining parents
max(tk.geno$prop_tsu)
min(tk.geno$prop_tsu) # 12.7%-81% Tsu

max(si.geno$prop_swe)
min(si.geno$prop_swe) # 0-1% Itl. Are some of them selfs? How did that happen?

# merge percent contribution with phenotypes
tk.perc <- as.data.frame(cbind(tk.geno$id, tk.geno$prop_tsu, tk.geno$prop_kas))
colnames(tk.perc) <- c("id", "prop_tsu", "prop_kas")
tk.pheno.perc <- merge(x=tk.pheno, y=tk.perc, by="id")

si.perc <- as.data.frame(cbind(si.geno$id, si.geno$prop_swe, si.geno$prop_itl))
colnames(si.perc) <- c("id", "prop_swe", "prop_itl")
si.pheno.perc <- merge(x=si.pheno, y=si.perc, by="id")

# Does parental genome contribution explain phenotype?
lm(data=si.pheno.perc, mean_stamen_number ~ prop_swe) %>% summary()
# percent genome contribution explains 18% variation in short stamen number in TK RILs

lm(data=tk.pheno.perc, mean_ss_num ~ prop_tsu) %>% summary() 
# percent genome contribution explains 15% variation in short stamen number in TK RILs


# test for correlation between flowering time and stamen number given percent genome contribution

# start with tests involving flowering time in Sweden chamber
si.pheno.perc.nona <- na.omit(si.pheno.perc)

pcor.test(x=si.pheno.perc.nona$mean_stamen_number, 
          y=si.pheno.perc.nona$SW_meandaystoflr, 
          z=si.pheno.perc.nona$prop_swe, method="pearson")
pcor.test(x=si.pheno.perc.nona$IT_meandaystoflr, 
          y=si.pheno.perc.nona$SW_meandaystoflr, 
          z=si.pheno.perc.nona$prop_swe, method="pearson")
pcor.test(x=si.pheno.perc.nona$ril_ovule_num, 
          y=si.pheno.perc.nona$SW_meandaystoflr, 
          z=si.pheno.perc.nona$prop_swe, method="pearson")

#remove days to flowering in Sweden chamber due to high rate of missingness
si.pheno.perc <- si.pheno.perc[,-5]
si.pheno.perc <-na.omit(si.pheno.perc)
pcor.test(x=si.pheno.perc$mean_stamen_number, 
          y=si.pheno.perc$IT_meandaystoflr, 
          z=si.pheno.perc$prop_swe, method="pearson")
cor.test(x=si.pheno.perc$mean_stamen_number, 
         y=si.pheno.perc$IT_meandaystoflr, method="pearson")

pcor.test(x=si.pheno.perc$mean_stamen_number, 
          y=si.pheno.perc$ril_ovule_num , 
          z=si.pheno.perc$prop_swe, method="pearson")
cor.test(x=si.pheno.perc$mean_stamen_number, 
         y=si.pheno.perc$ril_ovule_num, method="pearson")

pcor.test(x=si.pheno.perc$IT_meandaystoflr, 
          y=si.pheno.perc$ril_ovule_num , 
          z=si.pheno.perc$prop_swe, method="pearson")
cor.test(x=si.pheno.perc$IT_meandaystoflr, 
         y=si.pheno.perc$ril_ovule_num, method="pearson")

tk.pheno.perc <-na.omit(tk.pheno.perc)
pcor.test(x=tk.pheno.perc$mean_ss_num, 
          y=tk.pheno.perc$ft_raw, 
          z=tk.pheno.perc$prop_tsu, method="pearson")

tk.pheno.perc <-na.omit(tk.pheno.perc)
pcor.test(x=tk.pheno.perc$wue, 
          y=tk.pheno.perc$ft_raw, 
          z=tk.pheno.perc$prop_tsu, method="pearson")
cor.test(x=tk.pheno.perc$wue, 
         y=tk.pheno.perc$ft_raw, method="pearson")


pcor.test(x=tk.pheno.perc$ft_raw, 
          y=tk.pheno.perc$mean_ss_num, 
          z=tk.pheno.perc$prop_tsu, method="pearson")

pcor.test(x=tk.pheno.perc$wue, 
          y=tk.pheno.perc$mean_ss_num, 
          z=tk.pheno.perc$prop_tsu, method="pearson")

pcor.test(x=tk.pheno.perc$mean_ov_num, 
          y=tk.pheno.perc$ft_raw, 
          z=tk.pheno.perc$prop_tsu, method="pearson")

pcor.test(x=tk.pheno.perc$mean_ov_num, 
          y=tk.pheno.perc$wue, 
          z=tk.pheno.perc$prop_tsu, method="pearson")
