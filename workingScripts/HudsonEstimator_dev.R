# Fst Calculation Development Script

HudsonFst_slow <- function(n1, n2, p1, p2){

  numerator <- ((p1-p2)**2 - (p1*(1-p1))/(n1-1) - (p2*(1-p2))/(n2-1))

  denominator <- (p1*(1-p2) + p2*(1-p1))

  return(numerator/denominator)
}

# no assignment should make for faster calculation on large data sets.
HudsonFst <- function(n1, n2, p1, p2){
  return( ((p1-p2)**2 - (p1*(1-p1))/(n1-1) - (p2*(1-p2))/(n2-1)) / (p1*(1-p2) + p2*(1-p1)) )
}

# checking distribution rules in R

# (5 - 2)**2 = 3**2 = 9

(5-2)**2 # 9

# all distribution rules look good so far.

# TESTING: Hudson_slow and Hudson fast ... make sure we get the same results.

library(GWASpops.pheno2geno)
tml <- testMasterList
View(tml)

perAllele_1 <- tml[[2]][[1]]
perAllele_1
attributes(perAllele_1)

populations <- Populations
african_SampleCount <- as.numeric(populations$Sample_Count[2])
american_SampleCount <- as.numeric(populations$Sample_Count[10])


perAllele_1$population['1000GENOMES:phase_3:AFR']

AFR_freq <- perAllele_1[(perAllele_1$population == '1000GENOMES:phase_3:AFR') & (perAllele_1$allele != 'C'), ]$frequency
AMR_freq <- perAllele_1[(perAllele_1$population == '1000GENOMES:phase_3:AMR') & (perAllele_1$allele != 'C'), ]$frequency



hSlow <- HudsonFst_slow(african_SampleCount,american_SampleCount, AFR_freq, AMR_freq)
hSlow # 0.02688189

hFast <- HudsonFst(african_SampleCount,american_SampleCount, AFR_freq, AMR_freq)
hFast # 0.02688189

# Hudson estimations look equal and good. Will proceed with the Hudson fast.
