# Data integration: Taking data structures from fsst_GWAS_annotaion_lists and combining them together to create super data structures.
#
#
# Biggest problem with this plan is some DS (specifically the population data structures may not be small enough to smash into single DS while holding them in memory... but we will have to see.)



#Strategy is to do this piece wise, 1 sub process per data structure. Right now the only one that REALLY matters is the Fst one, so we will get that one fixed up first.

getwd()

fileList <- list.files('./workingData/fst_GWAS_annotation_lists/')
setwd("./workingData/fst_GWAS_annotation_lists/")

allFst <- lapply(fileList, function(i){
  load(i)
  return(masterList[[4]])
})

master_fst <- data.frame()

for(i in 1:length(allFst)){
  master_fst <- rbind.data.frame(master_fst, allFst[[i]])
}
# Success!


# We want to correct all negative values within it and make them 0. as That is the appropriate value for meaningful Fst.

na_count <- sum(is.na(master_fst))
na_count # 14,244,420

total_fst_slots <- 159559*496
na_count/total_fst_slots # 0.1799 ... so ~ 18% of slots are NA ... not terrible all things considered.

master_fst[master_fst < 0 ] <- 0

# Converted all neg values to 0. Now we need to store the new object in a new folder.. and consider what steps are next for the analysis we wish to do
full_fst <- master_fst

setwd("..")
setwd("./full_data_for_analysis/")
save(full_fst, file = "full_fst.rds")

# first DS saved, lets go for the next, larger 2. Gonna clear up some ram from Chrome first.

rm(allFst, allFst_df, fstDF_test, fstDF2_test, full_fst, gwasData, masterList, pops)

setwd("./workingData/fst_GWAS_annotation_lists/")


full_popAlleleFreq_perAllele <- lapply(fileList, function(i){
  load(i)
  return(masterList[[2]])
})
# 2 GB size... not bad. May actually be able to get them all smashed together in memory at once if its ever needed.

######************************************************************************************************************
fTest <- purrr::flatten(full_popAlleleFreq_perAllele) # 213,710 elements ... number of SNPs I succesfully got population data on here
###########################################
full_popAlleleFreq_perAllele <- fTest
save(full_popAlleleFreq_perAllele, file = "full_popAlleleFreq_PerAllele.rds")


# 3rd DS ### 1GB here...

full_popAlleleFreq_perPopulation <- lapply(fileList, function(i){
  load(i)
  return(masterList[[3]])
})


full_popAlleleFreq_perPopulation <- purrr::flatten(full_popAlleleFreq_perPopulation)
perPopNames <- names(full_popAlleleFreq_perPopulation)
uni_perPopNames <- unique(perPopNames) # 73 unique names.. good.
1679/73 # 23

# binding all tables with the same names together

testFilter <- full_popAlleleFreq_perPopulation[ perPopNames %in% uni_perPopNames[1] ]
test_rbind <- rbind.data.frame(testFilter) # same error... different number of rows... which is very strange.. shouldn't matter so long as cols are consistent, and they are
test_rbind <- cbind.data.frame(testFilter) # same error... doesn't matter that I am trying to do by col now.. weird

test_rbind <- dplyr::bind_rows(testFilter) # works... ok

full_perPopulation_popAlleleFreq <- lapply(uni_perPopNames, function(x){
  print(x)
  tempList <- full_popAlleleFreq_perPopulation[perPopNames %in% x] # filtering down to only single population
  fullSquash <- dplyr::bind_rows(tempList)
  return(fullSquash)
})

names(full_perPopulation_popAlleleFreq) <- uni_perPopNames # setting names now that the proper number of data frames exists.

# saving
save(full_perPopulation_popAlleleFreq, file = "full_perPopulation_popAlleleFreq.rds")

rm(full_popAlleleFreq_perPopulation)
rm(test_rbind, testFilter)


# 4th DS, the master Table  ## Less than 1 GB here.



full_SNP_Annotations_GWASc_Ensembl <- lapply(fileList, function(i){
  load(i)
  return(masterList[[1]])
})


test <- dplyr::bind_rows(full_SNP_Annotations_GWASc_Ensembl) # many repeated variants here... but unique in other points being reported. Thus I will not cut anything out of here until I understand it makes sense to in the future.

full_SNP_Annotations_GWASc_Ensembl <- test

#saving final DS
save(full_SNP_Annotations_GWASc_Ensembl, file = "full_SNP_Annotations_GWASc_Ensembl.rds")


# constructing the complete SUPER data structure.

full_fst <- master_fst
full_masterList <- list(full_SNP_Annotations_GWASc_Ensembl, full_popAlleleFreq_perAllele, full_perPopulation_popAlleleFreq, full_fst) # 4GB total Which is workable.

save(full_masterList, file = "full_masterList.rds")


# all data is setup and saved... Now onto the basic analyses.

rm(full_masterList)


# Data integration complete.































