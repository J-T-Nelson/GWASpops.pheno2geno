#Fst using the Hudson estimator from Hudson et al. 1992 and Bhatia et al. 2013
#Fst_hud = (((p1 - p2)**2) - ((p1*(1-p1))/(n1-1)) - ((p2*(1-p2))/(n2-1))) / ((p1*(1-p2)) + (p2*(1-p1)))

#loading packages
library(ggplot2)
library(scales)
library(dplyr)
library(tidyr)

#set directory
setwd("/Users/philipbaldassari/Desktop/zim-cos_downsampled")

#reading files
ChrX <- read.csv("downsampled_allele_freq_ChrX.csv")
head(ChrX)
Chr2L <- read.csv("downsampled_allele_freq_Chr2L.csv")
Chr2R <- read.csv("downsampled_allele_freq_Chr2R.csv")
Chr3L <- read.csv("downsampled_allele_freq_Chr3L.csv")
Chr3R <- read.csv("downsampled_allele_freq_Chr3R.csv")


#pairwise Fst ZS_RAL_ZI
ChrX_ZS_RAL_ZI_Fst <- ChrX %>%
  filter(n_ZS == 4 & n_RAL == 100 & n_ZI == 100) %>%
  mutate(ZS.vs.RAL_Fst = ((((maf_ZS - maf_RAL)**2) - ((maf_ZS*(1-maf_ZS))/(n_ZS-1)) - ((maf_RAL*(1-maf_RAL))/(n_RAL-1))) / ((maf_ZS*(1-maf_RAL)) + (maf_RAL*(1-maf_ZS))))) %>%
  mutate(ZS.vs.ZI_Fst = ((((maf_ZS - maf_ZI)**2) - ((maf_ZS*(1-maf_ZS))/(n_ZS-1)) - ((maf_ZI*(1-maf_ZI))/(n_ZI-1))) / ((maf_ZS*(1-maf_ZI)) + (maf_ZI*(1-maf_ZS))))) %>%
  mutate(ZS.vs.ZS_Fst = ((((maf_ZS - maf_ZS)**2) - ((maf_ZS*(1-maf_ZS))/(n_ZS-1)) - ((maf_ZS*(1-maf_ZS))/(n_ZS-1))) / ((maf_ZS*(1-maf_ZS)) + (maf_ZS*(1-maf_ZS))))) %>%
  select(Chrom, Site, ZS.vs.RAL_Fst, ZS.vs.ZI_Fst, ZS.vs.ZS_Fst)

ChrX_ZS_RAL_ZI_Fst[is.na(ChrX_ZS_RAL_ZI_Fst)] = 0

ChrX_ZS_RAL_ZI_Fst <- ChrX_ZS_RAL_ZI_Fst %>%
  mutate(Avg_Fst_ZS.vs.RAL_ZI = (ZS.vs.RAL_Fst + ZS.vs.ZI_Fst) / 2)



Chr2L_ZS_RAL_ZI_Fst <- Chr2L %>%
  filter(n_ZS == 4 & n_RAL == 100 & n_ZI == 100) %>%
  mutate(ZS.vs.RAL_Fst = ((((maf_ZS - maf_RAL)**2) - ((maf_ZS*(1-maf_ZS))/(n_ZS-1)) - ((maf_RAL*(1-maf_RAL))/(n_RAL-1))) / ((maf_ZS*(1-maf_RAL)) + (maf_RAL*(1-maf_ZS))))) %>%
  mutate(ZS.vs.ZI_Fst = ((((maf_ZS - maf_ZI)**2) - ((maf_ZS*(1-maf_ZS))/(n_ZS-1)) - ((maf_ZI*(1-maf_ZI))/(n_ZI-1))) / ((maf_ZS*(1-maf_ZI)) + (maf_ZI*(1-maf_ZS))))) %>%
  mutate(ZS.vs.ZS_Fst = ((((maf_ZS - maf_ZS)**2) - ((maf_ZS*(1-maf_ZS))/(n_ZS-1)) - ((maf_ZS*(1-maf_ZS))/(n_ZS-1))) / ((maf_ZS*(1-maf_ZS)) + (maf_ZS*(1-maf_ZS))))) %>%
  select(Chrom, Site, ZS.vs.RAL_Fst, ZS.vs.ZI_Fst, ZS.vs.ZS_Fst)

Chr2L_ZS_RAL_ZI_Fst[is.na(Chr2L_ZS_RAL_ZI_Fst)] = 0

Chr2L_ZS_RAL_ZI_Fst <- Chr2L_ZS_RAL_ZI_Fst %>%
  mutate(Avg_Fst_ZS.vs.RAL_ZI = (ZS.vs.RAL_Fst + ZS.vs.ZI_Fst) / 2)



Chr2R_ZS_RAL_ZI_Fst <- Chr2R %>%
  filter(n_ZS == 4 & n_RAL == 100 & n_ZI == 100) %>%
  mutate(ZS.vs.RAL_Fst = ((((maf_ZS - maf_RAL)**2) - ((maf_ZS*(1-maf_ZS))/(n_ZS-1)) - ((maf_RAL*(1-maf_RAL))/(n_RAL-1))) / ((maf_ZS*(1-maf_RAL)) + (maf_RAL*(1-maf_ZS))))) %>%
  mutate(ZS.vs.ZI_Fst = ((((maf_ZS - maf_ZI)**2) - ((maf_ZS*(1-maf_ZS))/(n_ZS-1)) - ((maf_ZI*(1-maf_ZI))/(n_ZI-1))) / ((maf_ZS*(1-maf_ZI)) + (maf_ZI*(1-maf_ZS))))) %>%
  mutate(ZS.vs.ZS_Fst = ((((maf_ZS - maf_ZS)**2) - ((maf_ZS*(1-maf_ZS))/(n_ZS-1)) - ((maf_ZS*(1-maf_ZS))/(n_ZS-1))) / ((maf_ZS*(1-maf_ZS)) + (maf_ZS*(1-maf_ZS))))) %>%
  select(Chrom, Site, ZS.vs.RAL_Fst, ZS.vs.ZI_Fst, ZS.vs.ZS_Fst)

Chr2R_ZS_RAL_ZI_Fst[is.na(Chr2R_ZS_RAL_ZI_Fst)] = 0

Chr2R_ZS_RAL_ZI_Fst <- Chr2R_ZS_RAL_ZI_Fst %>%
  mutate(Avg_Fst_ZS.vs.RAL_ZI = (ZS.vs.RAL_Fst + ZS.vs.ZI_Fst) / 2)



Chr3L_ZS_RAL_ZI_Fst <- Chr3L %>%
  filter(n_ZS == 4 & n_RAL == 100 & n_ZI == 100) %>%
  mutate(ZS.vs.RAL_Fst = ((((maf_ZS - maf_RAL)**2) - ((maf_ZS*(1-maf_ZS))/(n_ZS-1)) - ((maf_RAL*(1-maf_RAL))/(n_RAL-1))) / ((maf_ZS*(1-maf_RAL)) + (maf_RAL*(1-maf_ZS))))) %>%
  mutate(ZS.vs.ZI_Fst = ((((maf_ZS - maf_ZI)**2) - ((maf_ZS*(1-maf_ZS))/(n_ZS-1)) - ((maf_ZI*(1-maf_ZI))/(n_ZI-1))) / ((maf_ZS*(1-maf_ZI)) + (maf_ZI*(1-maf_ZS))))) %>%
  mutate(ZS.vs.ZS_Fst = ((((maf_ZS - maf_ZS)**2) - ((maf_ZS*(1-maf_ZS))/(n_ZS-1)) - ((maf_ZS*(1-maf_ZS))/(n_ZS-1))) / ((maf_ZS*(1-maf_ZS)) + (maf_ZS*(1-maf_ZS))))) %>%
  select(Chrom, Site, ZS.vs.RAL_Fst, ZS.vs.ZI_Fst, ZS.vs.ZS_Fst)

Chr3L_ZS_RAL_ZI_Fst[is.na(Chr3L_ZS_RAL_ZI_Fst)] = 0

Chr3L_ZS_RAL_ZI_Fst <- Chr3L_ZS_RAL_ZI_Fst %>%
  mutate(Avg_Fst_ZS.vs.RAL_ZI = (ZS.vs.RAL_Fst + ZS.vs.ZI_Fst) / 2)



Chr3R_ZS_RAL_ZI_Fst <- Chr3R %>%
  filter(n_ZS == 4 & n_RAL == 100 & n_ZI == 100) %>%
  mutate(ZS.vs.RAL_Fst = ((((maf_ZS - maf_RAL)**2) - ((maf_ZS*(1-maf_ZS))/(n_ZS-1)) - ((maf_RAL*(1-maf_RAL))/(n_RAL-1))) / ((maf_ZS*(1-maf_RAL)) + (maf_RAL*(1-maf_ZS))))) %>%
  mutate(ZS.vs.ZI_Fst = ((((maf_ZS - maf_ZI)**2) - ((maf_ZS*(1-maf_ZS))/(n_ZS-1)) - ((maf_ZI*(1-maf_ZI))/(n_ZI-1))) / ((maf_ZS*(1-maf_ZI)) + (maf_ZI*(1-maf_ZS))))) %>%
  mutate(ZS.vs.ZS_Fst = ((((maf_ZS - maf_ZS)**2) - ((maf_ZS*(1-maf_ZS))/(n_ZS-1)) - ((maf_ZS*(1-maf_ZS))/(n_ZS-1))) / ((maf_ZS*(1-maf_ZS)) + (maf_ZS*(1-maf_ZS))))) %>%
  select(Chrom, Site, ZS.vs.RAL_Fst, ZS.vs.ZI_Fst, ZS.vs.ZS_Fst)

Chr3R_ZS_RAL_ZI_Fst[is.na(Chr3R_ZS_RAL_ZI_Fst)] = 0

Chr3R_ZS_RAL_ZI_Fst <- Chr3R_ZS_RAL_ZI_Fst %>%
  mutate(Avg_Fst_ZS.vs.RAL_ZI = (ZS.vs.RAL_Fst + ZS.vs.ZI_Fst) / 2)



#export to csv files
write.csv(ChrX_ZS_RAL_ZI_Fst, "ChrX_ZS_RAL_ZI_Fst.csv", row.names = FALSE)
write.csv(Chr2L_ZS_RAL_ZI_Fst, "Chr2L_ZS_RAL_ZI_Fst.csv", row.names = FALSE)
write.csv(Chr2R_ZS_RAL_ZI_Fst, "Chr2R_ZS_RAL_ZI_Fst.csv", row.names = FALSE)
write.csv(Chr3L_ZS_RAL_ZI_Fst, "Chr3L_ZS_RAL_ZI_Fst.csv", row.names = FALSE)
write.csv(Chr3R_ZS_RAL_ZI_Fst, "Chr3R_ZS_RAL_ZI_Fst.csv", row.names = FALSE)



##############################################################################################################################################################################################################################




#pairwise Fst ZS_RAL_ZI_FR_SAfr
ChrX_ZS_RAL_ZI_FR_SAfr_Fst <- ChrX %>%
  filter(n_ZS == 4 & n_RAL == 100 & n_ZI == 100 & n_FR == 20 & n_SAfr == 20) %>%
  mutate(ZS.vs.RAL_Fst = ((((maf_ZS - maf_RAL)**2) - ((maf_ZS*(1-maf_ZS))/(n_ZS-1)) - ((maf_RAL*(1-maf_RAL))/(n_RAL-1))) / ((maf_ZS*(1-maf_RAL)) + (maf_RAL*(1-maf_ZS))))) %>%
  mutate(ZS.vs.ZI_Fst = ((((maf_ZS - maf_ZI)**2) - ((maf_ZS*(1-maf_ZS))/(n_ZS-1)) - ((maf_ZI*(1-maf_ZI))/(n_ZI-1))) / ((maf_ZS*(1-maf_ZI)) + (maf_ZI*(1-maf_ZS))))) %>%
  mutate(ZS.vs.FR_Fst = ((((maf_ZS - maf_FR)**2) - ((maf_ZS*(1-maf_ZS))/(n_ZS-1)) - ((maf_FR*(1-maf_FR))/(n_FR-1))) / ((maf_ZS*(1-maf_FR)) + (maf_FR*(1-maf_ZS))))) %>%
  mutate(ZS.vs.SAfr_Fst = ((((maf_ZS - maf_SAfr)**2) - ((maf_ZS*(1-maf_ZS))/(n_ZS-1)) - ((maf_SAfr*(1-maf_SAfr))/(n_SAfr-1))) / ((maf_ZS*(1-maf_SAfr)) + (maf_SAfr*(1-maf_ZS))))) %>%
  select(Chrom, Site, ZS.vs.RAL_Fst, ZS.vs.ZI_Fst, ZS.vs.FR_Fst, ZS.vs.SAfr_Fst)

ChrX_ZS_RAL_ZI_FR_SAfr_Fst[is.na(ChrX_ZS_RAL_ZI_FR_SAfr_Fst)] = 0

ChrX_ZS_RAL_ZI_FR_SAfr_Fst <- ChrX_ZS_RAL_ZI_FR_SAfr_Fst %>%
  mutate(Avg_Fst_ZS.vs.RAL_ZI_FR_SAfr = (ZS.vs.RAL_Fst + ZS.vs.ZI_Fst + ZS.vs.FR_Fst + ZS.vs.SAfr_Fst) / 4)



Chr2L_ZS_RAL_ZI_FR_SAfr_Fst <- Chr2L %>%
  filter(n_ZS == 4 & n_RAL == 100 & n_ZI == 100 & n_FR == 20 & n_SAfr == 20) %>%
  mutate(ZS.vs.RAL_Fst = ((((maf_ZS - maf_RAL)**2) - ((maf_ZS*(1-maf_ZS))/(n_ZS-1)) - ((maf_RAL*(1-maf_RAL))/(n_RAL-1))) / ((maf_ZS*(1-maf_RAL)) + (maf_RAL*(1-maf_ZS))))) %>%
  mutate(ZS.vs.ZI_Fst = ((((maf_ZS - maf_ZI)**2) - ((maf_ZS*(1-maf_ZS))/(n_ZS-1)) - ((maf_ZI*(1-maf_ZI))/(n_ZI-1))) / ((maf_ZS*(1-maf_ZI)) + (maf_ZI*(1-maf_ZS))))) %>%
  mutate(ZS.vs.FR_Fst = ((((maf_ZS - maf_FR)**2) - ((maf_ZS*(1-maf_ZS))/(n_ZS-1)) - ((maf_FR*(1-maf_FR))/(n_FR-1))) / ((maf_ZS*(1-maf_FR)) + (maf_FR*(1-maf_ZS))))) %>%
  mutate(ZS.vs.SAfr_Fst = ((((maf_ZS - maf_SAfr)**2) - ((maf_ZS*(1-maf_ZS))/(n_ZS-1)) - ((maf_SAfr*(1-maf_SAfr))/(n_SAfr-1))) / ((maf_ZS*(1-maf_SAfr)) + (maf_SAfr*(1-maf_ZS))))) %>%
  select(Chrom, Site, ZS.vs.RAL_Fst, ZS.vs.ZI_Fst, ZS.vs.FR_Fst, ZS.vs.SAfr_Fst)

Chr2L_ZS_RAL_ZI_FR_SAfr_Fst[is.na(Chr2L_ZS_RAL_ZI_FR_SAfr_Fst)] = 0

Chr2L_ZS_RAL_ZI_FR_SAfr_Fst <- Chr2L_ZS_RAL_ZI_FR_SAfr_Fst %>%
  mutate(Avg_Fst_ZS.vs.RAL_ZI_FR_SAfr = (ZS.vs.RAL_Fst + ZS.vs.ZI_Fst + ZS.vs.FR_Fst + ZS.vs.SAfr_Fst) / 4)



Chr2R_ZS_RAL_ZI_FR_SAfr_Fst <- Chr2R %>%
  filter(n_ZS == 4 & n_RAL == 100 & n_ZI == 100 & n_FR == 20 & n_SAfr == 20) %>%
  mutate(ZS.vs.RAL_Fst = ((((maf_ZS - maf_RAL)**2) - ((maf_ZS*(1-maf_ZS))/(n_ZS-1)) - ((maf_RAL*(1-maf_RAL))/(n_RAL-1))) / ((maf_ZS*(1-maf_RAL)) + (maf_RAL*(1-maf_ZS))))) %>%
  mutate(ZS.vs.ZI_Fst = ((((maf_ZS - maf_ZI)**2) - ((maf_ZS*(1-maf_ZS))/(n_ZS-1)) - ((maf_ZI*(1-maf_ZI))/(n_ZI-1))) / ((maf_ZS*(1-maf_ZI)) + (maf_ZI*(1-maf_ZS))))) %>%
  mutate(ZS.vs.FR_Fst = ((((maf_ZS - maf_FR)**2) - ((maf_ZS*(1-maf_ZS))/(n_ZS-1)) - ((maf_FR*(1-maf_FR))/(n_FR-1))) / ((maf_ZS*(1-maf_FR)) + (maf_FR*(1-maf_ZS))))) %>%
  mutate(ZS.vs.SAfr_Fst = ((((maf_ZS - maf_SAfr)**2) - ((maf_ZS*(1-maf_ZS))/(n_ZS-1)) - ((maf_SAfr*(1-maf_SAfr))/(n_SAfr-1))) / ((maf_ZS*(1-maf_SAfr)) + (maf_SAfr*(1-maf_ZS))))) %>%
  select(Chrom, Site, ZS.vs.RAL_Fst, ZS.vs.ZI_Fst, ZS.vs.FR_Fst, ZS.vs.SAfr_Fst)

Chr2R_ZS_RAL_ZI_FR_SAfr_Fst[is.na(Chr2R_ZS_RAL_ZI_FR_SAfr_Fst)] = 0

Chr2R_ZS_RAL_ZI_FR_SAfr_Fst <- Chr2R_ZS_RAL_ZI_FR_SAfr_Fst %>%
  mutate(Avg_Fst_ZS.vs.RAL_ZI_FR_SAfr = (ZS.vs.RAL_Fst + ZS.vs.ZI_Fst + ZS.vs.FR_Fst + ZS.vs.SAfr_Fst) / 4)



Chr3L_ZS_RAL_ZI_FR_SAfr_Fst <- Chr3L %>%
  filter(n_ZS == 4 & n_RAL == 100 & n_ZI == 100 & n_FR == 20 & n_SAfr == 20) %>%
  mutate(ZS.vs.RAL_Fst = ((((maf_ZS - maf_RAL)**2) - ((maf_ZS*(1-maf_ZS))/(n_ZS-1)) - ((maf_RAL*(1-maf_RAL))/(n_RAL-1))) / ((maf_ZS*(1-maf_RAL)) + (maf_RAL*(1-maf_ZS))))) %>%
  mutate(ZS.vs.ZI_Fst = ((((maf_ZS - maf_ZI)**2) - ((maf_ZS*(1-maf_ZS))/(n_ZS-1)) - ((maf_ZI*(1-maf_ZI))/(n_ZI-1))) / ((maf_ZS*(1-maf_ZI)) + (maf_ZI*(1-maf_ZS))))) %>%
  mutate(ZS.vs.FR_Fst = ((((maf_ZS - maf_FR)**2) - ((maf_ZS*(1-maf_ZS))/(n_ZS-1)) - ((maf_FR*(1-maf_FR))/(n_FR-1))) / ((maf_ZS*(1-maf_FR)) + (maf_FR*(1-maf_ZS))))) %>%
  mutate(ZS.vs.SAfr_Fst = ((((maf_ZS - maf_SAfr)**2) - ((maf_ZS*(1-maf_ZS))/(n_ZS-1)) - ((maf_SAfr*(1-maf_SAfr))/(n_SAfr-1))) / ((maf_ZS*(1-maf_SAfr)) + (maf_SAfr*(1-maf_ZS))))) %>%
  select(Chrom, Site, ZS.vs.RAL_Fst, ZS.vs.ZI_Fst, ZS.vs.FR_Fst, ZS.vs.SAfr_Fst)

Chr3L_ZS_RAL_ZI_FR_SAfr_Fst[is.na(Chr3L_ZS_RAL_ZI_FR_SAfr_Fst)] = 0

Chr3L_ZS_RAL_ZI_FR_SAfr_Fst <- Chr3L_ZS_RAL_ZI_FR_SAfr_Fst %>%
  mutate(Avg_Fst_ZS.vs.RAL_ZI_FR_SAfr = (ZS.vs.RAL_Fst + ZS.vs.ZI_Fst + ZS.vs.FR_Fst + ZS.vs.SAfr_Fst) / 4)



Chr3R_ZS_RAL_ZI_FR_SAfr_Fst <- Chr3R %>%
  filter(n_ZS == 4 & n_RAL == 100 & n_ZI == 100 & n_FR == 20 & n_SAfr == 20) %>%
  mutate(ZS.vs.RAL_Fst = ((((maf_ZS - maf_RAL)**2) - ((maf_ZS*(1-maf_ZS))/(n_ZS-1)) - ((maf_RAL*(1-maf_RAL))/(n_RAL-1))) / ((maf_ZS*(1-maf_RAL)) + (maf_RAL*(1-maf_ZS))))) %>%
  mutate(ZS.vs.ZI_Fst = ((((maf_ZS - maf_ZI)**2) - ((maf_ZS*(1-maf_ZS))/(n_ZS-1)) - ((maf_ZI*(1-maf_ZI))/(n_ZI-1))) / ((maf_ZS*(1-maf_ZI)) + (maf_ZI*(1-maf_ZS))))) %>%
  mutate(ZS.vs.FR_Fst = ((((maf_ZS - maf_FR)**2) - ((maf_ZS*(1-maf_ZS))/(n_ZS-1)) - ((maf_FR*(1-maf_FR))/(n_FR-1))) / ((maf_ZS*(1-maf_FR)) + (maf_FR*(1-maf_ZS))))) %>%
  mutate(ZS.vs.SAfr_Fst = ((((maf_ZS - maf_SAfr)**2) - ((maf_ZS*(1-maf_ZS))/(n_ZS-1)) - ((maf_SAfr*(1-maf_SAfr))/(n_SAfr-1))) / ((maf_ZS*(1-maf_SAfr)) + (maf_SAfr*(1-maf_ZS))))) %>%
  select(Chrom, Site, ZS.vs.RAL_Fst, ZS.vs.ZI_Fst, ZS.vs.FR_Fst, ZS.vs.SAfr_Fst)

Chr3R_ZS_RAL_ZI_FR_SAfr_Fst[is.na(Chr3R_ZS_RAL_ZI_FR_SAfr_Fst)] = 0

Chr3R_ZS_RAL_ZI_FR_SAfr_Fst <- Chr3R_ZS_RAL_ZI_FR_SAfr_Fst %>%
  mutate(Avg_Fst_ZS.vs.RAL_ZI_FR_SAfr = (ZS.vs.RAL_Fst + ZS.vs.ZI_Fst + ZS.vs.FR_Fst + ZS.vs.SAfr_Fst) / 4)



#export to csv files
write.csv(ChrX_ZS_RAL_ZI_FR_SAfr_Fst, "ChrX_ZS_RAL_ZI_FR_SAfr_Fst.csv", row.names = FALSE)
write.csv(Chr2L_ZS_RAL_ZI_FR_SAfr_Fst, "Chr2L_ZS_RAL_ZI_FR_SAfr_Fst.csv", row.names = FALSE)
write.csv(Chr2R_ZS_RAL_ZI_FR_SAfr_Fst, "Chr2R_ZS_RAL_ZI_FR_SAfr_Fst.csv", row.names = FALSE)
write.csv(Chr3L_ZS_RAL_ZI_FR_SAfr_Fst, "Chr3L_ZS_RAL_ZI_FR_SAfr_Fst.csv", row.names = FALSE)
write.csv(Chr3R_ZS_RAL_ZI_FR_SAfr_Fst, "Chr3R_ZS_RAL_ZI_FR_SAfr_Fst.csv", row.names = FALSE)



##############################################################################################################################################################################################################################




#pairwise Fst Zim_RAL_ZI
ChrX_Zim_RAL_ZI_Fst <- ChrX %>%
  filter(n_zim >= 8 & n_RAL == 100 & n_ZI == 100) %>%
  mutate(Zim.vs.RAL_Fst = ((((maf_zim - maf_RAL)**2) - ((maf_zim*(1-maf_zim))/(n_zim-1)) - ((maf_RAL*(1-maf_RAL))/(n_RAL-1))) / ((maf_zim*(1-maf_RAL)) + (maf_RAL*(1-maf_zim))))) %>%
  mutate(Zim.vs.ZI_Fst = ((((maf_zim - maf_ZI)**2) - ((maf_zim*(1-maf_zim))/(n_zim-1)) - ((maf_ZI*(1-maf_ZI))/(n_ZI-1))) / ((maf_zim*(1-maf_ZI)) + (maf_ZI*(1-maf_zim))))) %>%
  select(Chrom, Site, Zim.vs.RAL_Fst, Zim.vs.ZI_Fst)

ChrX_Zim_RAL_ZI_Fst[is.na(ChrX_Zim_RAL_ZI_Fst)] = 0

ChrX_Zim_RAL_ZI_Fst <- ChrX_Zim_RAL_ZI_Fst %>%
  mutate(Avg_Fst_Zim.vs.RAL_ZI = (Zim.vs.RAL_Fst + Zim.vs.ZI_Fst) / 2)



Chr2L_Zim_RAL_ZI_Fst <- Chr2L %>%
  filter(n_zim >= 8 & n_RAL == 100 & n_ZI == 100) %>%
  mutate(Zim.vs.RAL_Fst = ((((maf_zim - maf_RAL)**2) - ((maf_zim*(1-maf_zim))/(n_zim-1)) - ((maf_RAL*(1-maf_RAL))/(n_RAL-1))) / ((maf_zim*(1-maf_RAL)) + (maf_RAL*(1-maf_zim))))) %>%
  mutate(Zim.vs.ZI_Fst = ((((maf_zim - maf_ZI)**2) - ((maf_zim*(1-maf_zim))/(n_zim-1)) - ((maf_ZI*(1-maf_ZI))/(n_ZI-1))) / ((maf_zim*(1-maf_ZI)) + (maf_ZI*(1-maf_zim))))) %>%
  select(Chrom, Site, Zim.vs.RAL_Fst, Zim.vs.ZI_Fst)

Chr2L_Zim_RAL_ZI_Fst[is.na(Chr2L_Zim_RAL_ZI_Fst)] = 0

Chr2L_Zim_RAL_ZI_Fst <- Chr2L_Zim_RAL_ZI_Fst %>%
  mutate(Avg_Fst_Zim.vs.RAL_ZI = (Zim.vs.RAL_Fst + Zim.vs.ZI_Fst) / 2)



Chr2R_Zim_RAL_ZI_Fst <- Chr2R %>%
  filter(n_zim >= 8 & n_RAL == 100 & n_ZI == 100) %>%
  mutate(Zim.vs.RAL_Fst = ((((maf_zim - maf_RAL)**2) - ((maf_zim*(1-maf_zim))/(n_zim-1)) - ((maf_RAL*(1-maf_RAL))/(n_RAL-1))) / ((maf_zim*(1-maf_RAL)) + (maf_RAL*(1-maf_zim))))) %>%
  mutate(Zim.vs.ZI_Fst = ((((maf_zim - maf_ZI)**2) - ((maf_zim*(1-maf_zim))/(n_zim-1)) - ((maf_ZI*(1-maf_ZI))/(n_ZI-1))) / ((maf_zim*(1-maf_ZI)) + (maf_ZI*(1-maf_zim))))) %>%
  select(Chrom, Site, Zim.vs.RAL_Fst, Zim.vs.ZI_Fst)

Chr2R_Zim_RAL_ZI_Fst[is.na(Chr2R_Zim_RAL_ZI_Fst)] = 0

Chr2R_Zim_RAL_ZI_Fst <- Chr2R_Zim_RAL_ZI_Fst %>%
  mutate(Avg_Fst_Zim.vs.RAL_ZI = (Zim.vs.RAL_Fst + Zim.vs.ZI_Fst) / 2)



Chr3L_Zim_RAL_ZI_Fst <- Chr3L %>%
  filter(n_zim >= 8 & n_RAL == 100 & n_ZI == 100) %>%
  mutate(Zim.vs.RAL_Fst = ((((maf_zim - maf_RAL)**2) - ((maf_zim*(1-maf_zim))/(n_zim-1)) - ((maf_RAL*(1-maf_RAL))/(n_RAL-1))) / ((maf_zim*(1-maf_RAL)) + (maf_RAL*(1-maf_zim))))) %>%
  mutate(Zim.vs.ZI_Fst = ((((maf_zim - maf_ZI)**2) - ((maf_zim*(1-maf_zim))/(n_zim-1)) - ((maf_ZI*(1-maf_ZI))/(n_ZI-1))) / ((maf_zim*(1-maf_ZI)) + (maf_ZI*(1-maf_zim))))) %>%
  select(Chrom, Site, Zim.vs.RAL_Fst, Zim.vs.ZI_Fst)

Chr3L_Zim_RAL_ZI_Fst[is.na(Chr3L_Zim_RAL_ZI_Fst)] = 0

Chr3L_Zim_RAL_ZI_Fst <- Chr3L_Zim_RAL_ZI_Fst %>%
  mutate(Avg_Fst_Zim.vs.RAL_ZI = (Zim.vs.RAL_Fst + Zim.vs.ZI_Fst) / 2)



Chr3R_Zim_RAL_ZI_Fst <- Chr3R %>%
  filter(n_zim >= 8 & n_RAL == 100 & n_ZI == 100) %>%
  mutate(Zim.vs.RAL_Fst = ((((maf_zim - maf_RAL)**2) - ((maf_zim*(1-maf_zim))/(n_zim-1)) - ((maf_RAL*(1-maf_RAL))/(n_RAL-1))) / ((maf_zim*(1-maf_RAL)) + (maf_RAL*(1-maf_zim))))) %>%
  mutate(Zim.vs.ZI_Fst = ((((maf_zim - maf_ZI)**2) - ((maf_zim*(1-maf_zim))/(n_zim-1)) - ((maf_ZI*(1-maf_ZI))/(n_ZI-1))) / ((maf_zim*(1-maf_ZI)) + (maf_ZI*(1-maf_zim))))) %>%
  select(Chrom, Site, Zim.vs.RAL_Fst, Zim.vs.ZI_Fst)

Chr3R_Zim_RAL_ZI_Fst[is.na(Chr3R_Zim_RAL_ZI_Fst)] = 0

Chr3R_Zim_RAL_ZI_Fst <- Chr3R_Zim_RAL_ZI_Fst %>%
  mutate(Avg_Fst_Zim.vs.RAL_ZI = (Zim.vs.RAL_Fst + Zim.vs.ZI_Fst) / 2)



#export to csv files
write.csv(ChrX_Zim_RAL_ZI_Fst, "ChrX_Zim_RAL_ZI_Fst.csv", row.names = FALSE)
write.csv(Chr2L_Zim_RAL_ZI_Fst, "Chr2L_Zim_RAL_ZI_Fst.csv", row.names = FALSE)
write.csv(Chr2R_Zim_RAL_ZI_Fst, "Chr2R_Zim_RAL_ZI_Fst.csv", row.names = FALSE)
write.csv(Chr3L_Zim_RAL_ZI_Fst, "Chr3L_Zim_RAL_ZI_Fst.csv", row.names = FALSE)
write.csv(Chr3R_Zim_RAL_ZI_Fst, "Chr3R_Zim_RAL_ZI_Fst.csv", row.names = FALSE)




##############################################################################################################################################################################################################################




#pairwise Fst ZH_RAL_ZI
ChrX_ZH_RAL_ZI_Fst <- ChrX %>%
  filter(n_ZH == 3 & n_RAL == 100 & n_ZI == 100) %>%
  mutate(ZH.vs.RAL_Fst = ((((maf_ZH - maf_RAL)**2) - ((maf_ZH*(1-maf_ZH))/(n_ZH-1)) - ((maf_RAL*(1-maf_RAL))/(n_RAL-1))) / ((maf_ZH*(1-maf_RAL)) + (maf_RAL*(1-maf_ZH))))) %>%
  mutate(ZH.vs.ZI_Fst = ((((maf_ZH - maf_ZI)**2) - ((maf_ZH*(1-maf_ZH))/(n_ZH-1)) - ((maf_ZI*(1-maf_ZI))/(n_ZI-1))) / ((maf_ZH*(1-maf_ZI)) + (maf_ZI*(1-maf_ZH))))) %>%
  select(Chrom, Site, ZH.vs.RAL_Fst, ZH.vs.ZI_Fst)

ChrX_ZH_RAL_ZI_Fst[is.na(ChrX_ZH_RAL_ZI_Fst)] = 0

ChrX_ZH_RAL_ZI_Fst <- ChrX_ZH_RAL_ZI_Fst %>%
  mutate(Avg_Fst_ZH.vs.RAL_ZI = (ZH.vs.RAL_Fst + ZH.vs.ZI_Fst) / 2)



Chr2L_ZH_RAL_ZI_Fst <- Chr2L %>%
  filter(n_ZH == 3 & n_RAL == 100 & n_ZI == 100) %>%
  mutate(ZH.vs.RAL_Fst = ((((maf_ZH - maf_RAL)**2) - ((maf_ZH*(1-maf_ZH))/(n_ZH-1)) - ((maf_RAL*(1-maf_RAL))/(n_RAL-1))) / ((maf_ZH*(1-maf_RAL)) + (maf_RAL*(1-maf_ZH))))) %>%
  mutate(ZH.vs.ZI_Fst = ((((maf_ZH - maf_ZI)**2) - ((maf_ZH*(1-maf_ZH))/(n_ZH-1)) - ((maf_ZI*(1-maf_ZI))/(n_ZI-1))) / ((maf_ZH*(1-maf_ZI)) + (maf_ZI*(1-maf_ZH))))) %>%
  select(Chrom, Site, ZH.vs.RAL_Fst, ZH.vs.ZI_Fst)

Chr2L_ZH_RAL_ZI_Fst[is.na(Chr2L_ZH_RAL_ZI_Fst)] = 0

Chr2L_ZH_RAL_ZI_Fst <- Chr2L_ZH_RAL_ZI_Fst %>%
  mutate(Avg_Fst_ZH.vs.RAL_ZI = (ZH.vs.RAL_Fst + ZH.vs.ZI_Fst) / 2)



Chr2R_ZH_RAL_ZI_Fst <- Chr2R %>%
  filter(n_ZH == 3 & n_RAL == 100 & n_ZI == 100) %>%
  mutate(ZH.vs.RAL_Fst = ((((maf_ZH - maf_RAL)**2) - ((maf_ZH*(1-maf_ZH))/(n_ZH-1)) - ((maf_RAL*(1-maf_RAL))/(n_RAL-1))) / ((maf_ZH*(1-maf_RAL)) + (maf_RAL*(1-maf_ZH))))) %>%
  mutate(ZH.vs.ZI_Fst = ((((maf_ZH - maf_ZI)**2) - ((maf_ZH*(1-maf_ZH))/(n_ZH-1)) - ((maf_ZI*(1-maf_ZI))/(n_ZI-1))) / ((maf_ZH*(1-maf_ZI)) + (maf_ZI*(1-maf_ZH))))) %>%
  select(Chrom, Site, ZH.vs.RAL_Fst, ZH.vs.ZI_Fst)

Chr2R_ZH_RAL_ZI_Fst[is.na(Chr2R_ZH_RAL_ZI_Fst)] = 0

Chr2R_ZH_RAL_ZI_Fst <- Chr2R_ZH_RAL_ZI_Fst %>%
  mutate(Avg_Fst_ZH.vs.RAL_ZI = (ZH.vs.RAL_Fst + ZH.vs.ZI_Fst) / 2)



Chr3L_ZH_RAL_ZI_Fst <- Chr3L %>%
  filter(n_ZH == 3 & n_RAL == 100 & n_ZI == 100) %>%
  mutate(ZH.vs.RAL_Fst = ((((maf_ZH - maf_RAL)**2) - ((maf_ZH*(1-maf_ZH))/(n_ZH-1)) - ((maf_RAL*(1-maf_RAL))/(n_RAL-1))) / ((maf_ZH*(1-maf_RAL)) + (maf_RAL*(1-maf_ZH))))) %>%
  mutate(ZH.vs.ZI_Fst = ((((maf_ZH - maf_ZI)**2) - ((maf_ZH*(1-maf_ZH))/(n_ZH-1)) - ((maf_ZI*(1-maf_ZI))/(n_ZI-1))) / ((maf_ZH*(1-maf_ZI)) + (maf_ZI*(1-maf_ZH))))) %>%
  select(Chrom, Site, ZH.vs.RAL_Fst, ZH.vs.ZI_Fst)

Chr3L_ZH_RAL_ZI_Fst[is.na(Chr3L_ZH_RAL_ZI_Fst)] = 0

Chr3L_ZH_RAL_ZI_Fst <- Chr3L_ZH_RAL_ZI_Fst %>%
  mutate(Avg_Fst_ZH.vs.RAL_ZI = (ZH.vs.RAL_Fst + ZH.vs.ZI_Fst) / 2)



Chr3R_ZH_RAL_ZI_Fst <- Chr3R %>%
  filter(n_ZH == 3 & n_RAL == 100 & n_ZI == 100) %>%
  mutate(ZH.vs.RAL_Fst = ((((maf_ZH - maf_RAL)**2) - ((maf_ZH*(1-maf_ZH))/(n_ZH-1)) - ((maf_RAL*(1-maf_RAL))/(n_RAL-1))) / ((maf_ZH*(1-maf_RAL)) + (maf_RAL*(1-maf_ZH))))) %>%
  mutate(ZH.vs.ZI_Fst = ((((maf_ZH - maf_ZI)**2) - ((maf_ZH*(1-maf_ZH))/(n_ZH-1)) - ((maf_ZI*(1-maf_ZI))/(n_ZI-1))) / ((maf_ZH*(1-maf_ZI)) + (maf_ZI*(1-maf_ZH))))) %>%
  select(Chrom, Site, ZH.vs.RAL_Fst, ZH.vs.ZI_Fst)

Chr3R_ZH_RAL_ZI_Fst[is.na(Chr3R_ZH_RAL_ZI_Fst)] = 0

Chr3R_ZH_RAL_ZI_Fst <- Chr3R_ZH_RAL_ZI_Fst %>%
  mutate(Avg_Fst_ZH.vs.RAL_ZI = (ZH.vs.RAL_Fst + ZH.vs.ZI_Fst) / 2)



#export to csv files
write.csv(ChrX_ZH_RAL_ZI_Fst, "ChrX_ZH_RAL_ZI_Fst.csv", row.names = FALSE)
write.csv(Chr2L_ZH_RAL_ZI_Fst, "Chr2L_ZH_RAL_ZI_Fst.csv", row.names = FALSE)
write.csv(Chr2R_ZH_RAL_ZI_Fst, "Chr2R_ZH_RAL_ZI_Fst.csv", row.names = FALSE)
write.csv(Chr3L_ZH_RAL_ZI_Fst, "Chr3L_ZH_RAL_ZI_Fst.csv", row.names = FALSE)
write.csv(Chr3R_ZH_RAL_ZI_Fst, "Chr3R_ZH_RAL_ZI_Fst.csv", row.names = FALSE)




##############################################################################################################################################################################################################################





#pairwise Fst ZW_RAL_ZI
ChrX_ZW_RAL_ZI_Fst <- ChrX %>%
  filter(n_ZW == 4 & n_RAL == 100 & n_ZI == 100) %>%
  mutate(ZW.vs.RAL_Fst = ((((maf_ZW - maf_RAL)**2) - ((maf_ZW*(1-maf_ZW))/(n_ZW-1)) - ((maf_RAL*(1-maf_RAL))/(n_RAL-1))) / ((maf_ZW*(1-maf_RAL)) + (maf_RAL*(1-maf_ZW))))) %>%
  mutate(ZW.vs.ZI_Fst = ((((maf_ZW - maf_ZI)**2) - ((maf_ZW*(1-maf_ZW))/(n_ZW-1)) - ((maf_ZI*(1-maf_ZI))/(n_ZI-1))) / ((maf_ZW*(1-maf_ZI)) + (maf_ZI*(1-maf_ZW))))) %>%
  select(Chrom, Site, ZW.vs.RAL_Fst, ZW.vs.ZI_Fst)

ChrX_ZW_RAL_ZI_Fst[is.na(ChrX_ZW_RAL_ZI_Fst)] = 0

ChrX_ZW_RAL_ZI_Fst <- ChrX_ZW_RAL_ZI_Fst %>%
  mutate(Avg_Fst_ZW.vs.RAL_ZI = (ZW.vs.RAL_Fst + ZW.vs.ZI_Fst) / 2)



Chr2L_ZW_RAL_ZI_Fst <- Chr2L %>%
  filter(n_ZW == 4 & n_RAL == 100 & n_ZI == 100) %>%
  mutate(ZW.vs.RAL_Fst = ((((maf_ZW - maf_RAL)**2) - ((maf_ZW*(1-maf_ZW))/(n_ZW-1)) - ((maf_RAL*(1-maf_RAL))/(n_RAL-1))) / ((maf_ZW*(1-maf_RAL)) + (maf_RAL*(1-maf_ZW))))) %>%
  mutate(ZW.vs.ZI_Fst = ((((maf_ZW - maf_ZI)**2) - ((maf_ZW*(1-maf_ZW))/(n_ZW-1)) - ((maf_ZI*(1-maf_ZI))/(n_ZI-1))) / ((maf_ZW*(1-maf_ZI)) + (maf_ZI*(1-maf_ZW))))) %>%
  select(Chrom, Site, ZW.vs.RAL_Fst, ZW.vs.ZI_Fst)

Chr2L_ZW_RAL_ZI_Fst[is.na(Chr2L_ZW_RAL_ZI_Fst)] = 0

Chr2L_ZW_RAL_ZI_Fst <- Chr2L_ZW_RAL_ZI_Fst %>%
  mutate(Avg_Fst_ZW.vs.RAL_ZI = (ZW.vs.RAL_Fst + ZW.vs.ZI_Fst) / 2)



Chr2R_ZW_RAL_ZI_Fst <- Chr2R %>%
  filter(n_ZW == 4 & n_RAL == 100 & n_ZI == 100) %>%
  mutate(ZW.vs.RAL_Fst = ((((maf_ZW - maf_RAL)**2) - ((maf_ZW*(1-maf_ZW))/(n_ZW-1)) - ((maf_RAL*(1-maf_RAL))/(n_RAL-1))) / ((maf_ZW*(1-maf_RAL)) + (maf_RAL*(1-maf_ZW))))) %>%
  mutate(ZW.vs.ZI_Fst = ((((maf_ZW - maf_ZI)**2) - ((maf_ZW*(1-maf_ZW))/(n_ZW-1)) - ((maf_ZI*(1-maf_ZI))/(n_ZI-1))) / ((maf_ZW*(1-maf_ZI)) + (maf_ZI*(1-maf_ZW))))) %>%
  select(Chrom, Site, ZW.vs.RAL_Fst, ZW.vs.ZI_Fst)

Chr2R_ZW_RAL_ZI_Fst[is.na(Chr2R_ZW_RAL_ZI_Fst)] = 0

Chr2R_ZW_RAL_ZI_Fst <- Chr2R_ZW_RAL_ZI_Fst %>%
  mutate(Avg_Fst_ZW.vs.RAL_ZI = (ZW.vs.RAL_Fst + ZW.vs.ZI_Fst) / 2)



Chr3L_ZW_RAL_ZI_Fst <- Chr3L %>%
  filter(n_ZW == 4 & n_RAL == 100 & n_ZI == 100) %>%
  mutate(ZW.vs.RAL_Fst = ((((maf_ZW - maf_RAL)**2) - ((maf_ZW*(1-maf_ZW))/(n_ZW-1)) - ((maf_RAL*(1-maf_RAL))/(n_RAL-1))) / ((maf_ZW*(1-maf_RAL)) + (maf_RAL*(1-maf_ZW))))) %>%
  mutate(ZW.vs.ZI_Fst = ((((maf_ZW - maf_ZI)**2) - ((maf_ZW*(1-maf_ZW))/(n_ZW-1)) - ((maf_ZI*(1-maf_ZI))/(n_ZI-1))) / ((maf_ZW*(1-maf_ZI)) + (maf_ZI*(1-maf_ZW))))) %>%
  select(Chrom, Site, ZW.vs.RAL_Fst, ZW.vs.ZI_Fst)

Chr3L_ZW_RAL_ZI_Fst[is.na(Chr3L_ZW_RAL_ZI_Fst)] = 0

Chr3L_ZW_RAL_ZI_Fst <- Chr3L_ZW_RAL_ZI_Fst %>%
  mutate(Avg_Fst_ZW.vs.RAL_ZI = (ZW.vs.RAL_Fst + ZW.vs.ZI_Fst) / 2)



Chr3R_ZW_RAL_ZI_Fst <- Chr3R %>%
  filter(n_ZW == 4 & n_RAL == 100 & n_ZI == 100) %>%
  mutate(ZW.vs.RAL_Fst = ((((maf_ZW - maf_RAL)**2) - ((maf_ZW*(1-maf_ZW))/(n_ZW-1)) - ((maf_RAL*(1-maf_RAL))/(n_RAL-1))) / ((maf_ZW*(1-maf_RAL)) + (maf_RAL*(1-maf_ZW))))) %>%
  mutate(ZW.vs.ZI_Fst = ((((maf_ZW - maf_ZI)**2) - ((maf_ZW*(1-maf_ZW))/(n_ZW-1)) - ((maf_ZI*(1-maf_ZI))/(n_ZI-1))) / ((maf_ZW*(1-maf_ZI)) + (maf_ZI*(1-maf_ZW))))) %>%
  select(Chrom, Site, ZW.vs.RAL_Fst, ZW.vs.ZI_Fst)

Chr3R_ZW_RAL_ZI_Fst[is.na(Chr3R_ZW_RAL_ZI_Fst)] = 0

Chr3R_ZW_RAL_ZI_Fst <- Chr3R_ZW_RAL_ZI_Fst %>%
  mutate(Avg_Fst_ZW.vs.RAL_ZI = (ZW.vs.RAL_Fst + ZW.vs.ZI_Fst) / 2)



#export to csv files
write.csv(ChrX_ZW_RAL_ZI_Fst, "ChrX_ZW_RAL_ZI_Fst.csv", row.names = FALSE)
write.csv(Chr2L_ZW_RAL_ZI_Fst, "Chr2L_ZW_RAL_ZI_Fst.csv", row.names = FALSE)
write.csv(Chr2R_ZW_RAL_ZI_Fst, "Chr2R_ZW_RAL_ZI_Fst.csv", row.names = FALSE)
write.csv(Chr3L_ZW_RAL_ZI_Fst, "Chr3L_ZW_RAL_ZI_Fst.csv", row.names = FALSE)
write.csv(Chr3R_ZW_RAL_ZI_Fst, "Chr3R_ZW_RAL_ZI_Fst.csv", row.names = FALSE)




##############################################################################################################################################################################################################################





#pairwise Fst ZS_ZH_ZW
ChrX_ZS_ZH_ZW_Fst <- ChrX %>%
  filter(n_ZS == 4 & n_ZH == 3 & n_ZW == 4) %>%
  mutate(ZS.vs.ZH_Fst = ((((maf_ZS - maf_ZH)**2) - ((maf_ZS*(1-maf_ZS))/(n_ZS-1)) - ((maf_ZH*(1-maf_ZH))/(n_ZH-1))) / ((maf_ZS*(1-maf_ZH)) + (maf_ZH*(1-maf_ZS))))) %>%
  mutate(ZS.vs.ZW_Fst = ((((maf_ZS - maf_ZW)**2) - ((maf_ZS*(1-maf_ZS))/(n_ZS-1)) - ((maf_ZW*(1-maf_ZW))/(n_ZW-1))) / ((maf_ZS*(1-maf_ZW)) + (maf_ZW*(1-maf_ZS))))) %>%
  select(Chrom, Site, ZS.vs.ZH_Fst, ZS.vs.ZW_Fst)

ChrX_ZS_ZH_ZW_Fst[is.na(ChrX_ZS_ZH_ZW_Fst)] = 0

ChrX_ZS_ZH_ZW_Fst <- ChrX_ZS_ZH_ZW_Fst %>%
  mutate(Avg_Fst_ZS.vs.ZH_ZW = (ZS.vs.ZH_Fst + ZS.vs.ZW_Fst) / 2)



Chr2L_ZS_ZH_ZW_Fst <- Chr2L %>%
  filter(n_ZS == 4 & n_ZH == 3 & n_ZW == 4) %>%
  mutate(ZS.vs.ZH_Fst = ((((maf_ZS - maf_ZH)**2) - ((maf_ZS*(1-maf_ZS))/(n_ZS-1)) - ((maf_ZH*(1-maf_ZH))/(n_ZH-1))) / ((maf_ZS*(1-maf_ZH)) + (maf_ZH*(1-maf_ZS))))) %>%
  mutate(ZS.vs.ZW_Fst = ((((maf_ZS - maf_ZW)**2) - ((maf_ZS*(1-maf_ZS))/(n_ZS-1)) - ((maf_ZW*(1-maf_ZW))/(n_ZW-1))) / ((maf_ZS*(1-maf_ZW)) + (maf_ZW*(1-maf_ZS))))) %>%
  select(Chrom, Site, ZS.vs.ZH_Fst, ZS.vs.ZW_Fst)

Chr2L_ZS_ZH_ZW_Fst[is.na(Chr2L_ZS_ZH_ZW_Fst)] = 0

Chr2L_ZS_ZH_ZW_Fst <- Chr2L_ZS_ZH_ZW_Fst %>%
  mutate(Avg_Fst_ZS.vs.ZH_ZW = (ZS.vs.ZH_Fst + ZS.vs.ZW_Fst) / 2)



Chr2R_ZS_ZH_ZW_Fst <- Chr2R %>%
  filter(n_ZS == 4 & n_ZH == 3 & n_ZW == 4) %>%
  mutate(ZS.vs.ZH_Fst = ((((maf_ZS - maf_ZH)**2) - ((maf_ZS*(1-maf_ZS))/(n_ZS-1)) - ((maf_ZH*(1-maf_ZH))/(n_ZH-1))) / ((maf_ZS*(1-maf_ZH)) + (maf_ZH*(1-maf_ZS))))) %>%
  mutate(ZS.vs.ZW_Fst = ((((maf_ZS - maf_ZW)**2) - ((maf_ZS*(1-maf_ZS))/(n_ZS-1)) - ((maf_ZW*(1-maf_ZW))/(n_ZW-1))) / ((maf_ZS*(1-maf_ZW)) + (maf_ZW*(1-maf_ZS))))) %>%
  select(Chrom, Site, ZS.vs.ZH_Fst, ZS.vs.ZW_Fst)

Chr2R_ZS_ZH_ZW_Fst[is.na(Chr2R_ZS_ZH_ZW_Fst)] = 0

Chr2R_ZS_ZH_ZW_Fst <- Chr2R_ZS_ZH_ZW_Fst %>%
  mutate(Avg_Fst_ZS.vs.ZH_ZW = (ZS.vs.ZH_Fst + ZS.vs.ZW_Fst) / 2)



Chr3L_ZS_ZH_ZW_Fst <- Chr3L %>%
  filter(n_ZS == 4 & n_ZH == 3 & n_ZW == 4) %>%
  mutate(ZS.vs.ZH_Fst = ((((maf_ZS - maf_ZH)**2) - ((maf_ZS*(1-maf_ZS))/(n_ZS-1)) - ((maf_ZH*(1-maf_ZH))/(n_ZH-1))) / ((maf_ZS*(1-maf_ZH)) + (maf_ZH*(1-maf_ZS))))) %>%
  mutate(ZS.vs.ZW_Fst = ((((maf_ZS - maf_ZW)**2) - ((maf_ZS*(1-maf_ZS))/(n_ZS-1)) - ((maf_ZW*(1-maf_ZW))/(n_ZW-1))) / ((maf_ZS*(1-maf_ZW)) + (maf_ZW*(1-maf_ZS))))) %>%
  select(Chrom, Site, ZS.vs.ZH_Fst, ZS.vs.ZW_Fst)

Chr3L_ZS_ZH_ZW_Fst[is.na(Chr3L_ZS_ZH_ZW_Fst)] = 0

Chr3L_ZS_ZH_ZW_Fst <- Chr3L_ZS_ZH_ZW_Fst %>%
  mutate(Avg_Fst_ZS.vs.ZH_ZW = (ZS.vs.ZH_Fst + ZS.vs.ZW_Fst) / 2)



Chr3R_ZS_ZH_ZW_Fst <- Chr3R %>%
  filter(n_ZS == 4 & n_ZH == 3 & n_ZW == 4) %>%
  mutate(ZS.vs.ZH_Fst = ((((maf_ZS - maf_ZH)**2) - ((maf_ZS*(1-maf_ZS))/(n_ZS-1)) - ((maf_ZH*(1-maf_ZH))/(n_ZH-1))) / ((maf_ZS*(1-maf_ZH)) + (maf_ZH*(1-maf_ZS))))) %>%
  mutate(ZS.vs.ZW_Fst = ((((maf_ZS - maf_ZW)**2) - ((maf_ZS*(1-maf_ZS))/(n_ZS-1)) - ((maf_ZW*(1-maf_ZW))/(n_ZW-1))) / ((maf_ZS*(1-maf_ZW)) + (maf_ZW*(1-maf_ZS))))) %>%
  select(Chrom, Site, ZS.vs.ZH_Fst, ZS.vs.ZW_Fst)

Chr3R_ZS_ZH_ZW_Fst[is.na(Chr3R_ZS_ZH_ZW_Fst)] = 0

Chr3R_ZS_ZH_ZW_Fst <- Chr3R_ZS_ZH_ZW_Fst %>%
  mutate(Avg_Fst_ZS.vs.ZH_ZW = (ZS.vs.ZH_Fst + ZS.vs.ZW_Fst) / 2)



#export to csv files
write.csv(ChrX_ZS_ZH_ZW_Fst, "ChrX_ZS_ZH_ZW_Fst.csv", row.names = FALSE)
write.csv(Chr2L_ZS_ZH_ZW_Fst, "Chr2L_ZS_ZH_ZW_Fst.csv", row.names = FALSE)
write.csv(Chr2R_ZS_ZH_ZW_Fst, "Chr2R_ZS_ZH_ZW_Fst.csv", row.names = FALSE)
write.csv(Chr3L_ZS_ZH_ZW_Fst, "Chr3L_ZS_ZH_ZW_Fst.csv", row.names = FALSE)
write.csv(Chr3R_ZS_ZH_ZW_Fst, "Chr3R_ZS_ZH_ZW_Fst.csv", row.names = FALSE)













