# novel analyses after recalculation of Fst and update of SNP annotation table with parental terms for categorical splitting

# 3-21-23


# Split out disease parent terms
# perform analyses:
#   1. Fst sum on all pop-pairs for all SNPs and just disease SNPs, Generate visualization for top 1,5,10% ; generate and save table for top 1,5,10%
#   2. Fst Sum on all pops against 1000Genomes:All population, for all SNPs and just disease SNPs, Generate visualization for top 1,5,10% ; generate and save table for top 1,5,10%
#   3. Fst sum for all pops split up into 32 groups :: for all SNPs and just disease SNPs, Generate visualization for top 1,5,10% ; generate and save table for top 1,5,10%
#   4. Brainstorm on further reportings / findings... figure out how to share info more meaningfully... figure out how to represent the intersection data of which SNPs are most significant to which populations ... Consider what sorts of reports best address your specific research questions


unique(EFO_term_mapping$`Parent term`)
# [1] "Cancer"                           "NR"                               "Other disease"                    "Digestive system disorder"
# [5] "Other trait"                      "Other measurement"                "Cardiovascular disease"           "Neurological disorder"
# [9] "Immune system disorder"           "Cardiovascular measurement"       "Hematological measurement"        "Biological process"
# [13] "Body measurement"                 "Lipid or lipoprotein measurement" "Metabolic disorder"               "Inflammatory measurement"
# [17] "Liver enzyme measurement"

diseaseParentTerms <- c('Cancer', 'Other disease', 'Digestive system disorder', 'Cardiovascular disease', 'Neurological disorder', 'Immune system disorder', 'Metabolic disorder')

SNPannoTable_disease <- full_SNP_Anno_withParentalTerms[full_SNP_Anno_withParentalTerms$`Parent term` %in% diseaseParentTerms ,]
# 77221 obs

uniDiseaseSNPs <- unique(SNPannoTable_disease$VariantID)
length(uniDiseaseSNPs) # 47,363 unique SNPS
numDiseaseSNP_withFst <- sum(uniDiseaseSNPs %in% rownames(full_fst_fixed)) # 36,097  .... reason for loss of values is unclear, however synonyms may be the cause, or once source of loss.
length(unique(rownames(full_fst_fixed))) # 159271 .. all row names unique.

36097 / 159271 # 0.2266389 => ~ 23% of our Fst values calculated come from disease associated traits.

# SANITY CHECK ...
numTotalSNP_withFst <- sum(unique(full_SNP_Anno_withParentalTerms$VariantID) %in% rownames(full_fst_fixed))
numTotalSNP_withFst # 158173 ... hm some are missing. weird... but for the most part we see everything. Wonder if the creation of full_SNP_Annotations_GWASc_Ensembl was done improperly in any way... we really shouldn't see any values missing even if its only 1.1k or so

uniqueSNP_asso <- unique(asso$SNPS)
uniqueSNP_postEnsemblDataGrab <- unique(full_SNP_Annotations_GWASc_Ensembl$VariantID)
# we expect to see substantial loss here... as many SNPs were not actually grabbed due to failed data transfer that was never solved.
length(uniqueSNP_asso) # 249,792
length(uniqueSNP_postEnsemblDataGrab) # 209,455
249792 - 209455 # 40337 # CRITICAL VALUE ... # SNPS EITHER NOT RETURNED FROM ENSEMBL OR NOT CONTAINED WITHIN.. ALL WERE REQUESTED.
209455/249792 # 0.8385176



# Begin analyses: Generate sums
diseaseFst_fixed <- full_fst_fixed[rownames(full_fst_fixed) %in% uniDiseaseSNPs,]

# Each pop against ALL pop
#popsVSall is the analysis of each population against the 1000kGenomes:ALL population
allFst_popsVSall <- full_fst_fixed[, grep("ALL", colnames(full_fst_fixed) )]

diseaseFst_popsVSall <- diseaseFst_fixed[, grep("ALL", colnames(diseaseFst_fixed) )]

# each pop on its own
allFst_allPops_list <- lapply(thousGenPops$Population_Abbreviation, function(i){ full_fst_fixed[, grep(i, colnames(full_fst_fixed))] })
diseaseFst_allPops_list <- lapply(thousGenPops$Population_Abbreviation, function(i){ diseaseFst_fixed[, grep(i, colnames(diseaseFst_fixed))] })


# Summing Fst for analysis on all subsets:
allSNPfst_sum <- sapply(1:nrow(full_fst_fixed), function(i){ return(sum(full_fst_fixed[i, ])) })
allDiseaseSNPfst_sum <- sapply(1:nrow(diseaseFst_fixed), function(i){ return(sum(diseaseFst_fixed[i, ])) })

allFst_popsVsall_sum <- sapply(1:nrow(allFst_popsVSall), function(i){ return(sum(allFst_popsVSall[i, ])) })
diseaseFst_popsVsall_sum <- sapply(1:nrow(diseaseFst_popsVSall), function(i){ return(sum(diseaseFst_popsVSall[i, ])) })


allFst_allPops_list_sum <- lapply(allFst_allPops_list, function(x){
  sapply(1:nrow(x), function(i){
    return(sum(x[i, ]))
    })
})
diseaseFst_allPops_list_sum <- lapply(diseaseFst_allPops_list, function(x){
  sapply(1:nrow(x), function(i){
    return(sum(x[i, ]))
  })
})

# Saving summed Fst vals:

getwd()
setwd("./workingData/fst_sums")
save(allSNPfst_sum, file="allSNPfst_sum.rds")
save(allDiseaseSNPfst_sum, file="allDiseaseSNPfst_sum.rds")
save(allFst_popsVsall_sum, file="allFst_popsVsall_sum.rds")
save(diseaseFst_popsVsall_sum, file="diseaseFst_popsVsall_sum.rds")
save(allFst_allPops_list_sum, file="allFst_allPops_list_sum.rds")
save(diseaseFst_allPops_list_sum, file="diseaseFst_allPops_list_sum.rds")
setwd("..")

# Generating plots and Tables:


##### QUESTION: Are there genes involved in human phenotypes/disease harboring SNPs which are population-specific?

# to answer this we will get the top 10% SNPs for each popVS_ALLpop, and for each pop itself

# I think we may not need to actually sum across the cols of the ALLpopVSothers DF... as each column itself should capture the top SNPs per population against the metapopulation.. which in itself may be interesting.

# 3-22 ... Start with generation of top 1,5,10% tables. Graph them for visualization and understanding. Merge interesting features back into those tables. ... figure out whats next


allDiseaseSNPfst_sum[1:10]
allFst_popsVsall_sum[1:10]
allSNPfst_sum[1:10] # not seeing names on values of these vectors... need to add that into them before doing anything else with them
diseaseFst_popsVsall_sum[1:10]

rownames(diseaseFst_fixed)[1:10]

names(allDiseaseSNPfst_sum) <- rownames(diseaseFst_fixed)
names(diseaseFst_popsVsall_sum) <- rownames(diseaseFst_fixed)
names(allSNPfst_sum) <- rownames(full_fst_fixed)
names(allFst_popsVsall_sum) <- rownames(full_fst_fixed)

# apply names to list sums
for(i in 1:length(allFst_allPops_list_sum)){
  names(allFst_allPops_list_sum[[i]]) <- rownames(allFst_allPops_list[[i]])
}


for(i in 1:length(diseaseFst_allPops_list_sum)){
  names(diseaseFst_allPops_list_sum[[i]]) <- rownames(diseaseFst_allPops_list[[i]])
}

names(allFst_allPops_list_sum) <- thousGenPops$Population_Abbreviation
names(diseaseFst_allPops_list_sum) <- thousGenPops$Population_Abbreviation

# All names have been applied now... Going to resave all data

# Generate top tables:

?quantile
quantile(x <- rnorm(1001)) # Extremes & Quartiles by default
quantile(x,  probs = c(0.1, 0.5, 1, 2, 5, 10, 50, NA)/100)

# get quantiles for filtering
allSNP_quantiles <- quantile(allSNPfst_sum, probs=c(.99,.95,.90))
allSNP_popsVSall_quantiles <- quantile(allFst_popsVsall_sum, probs=c(.99,.95,.90))
allDiseaseSNPfst_quantiles <- quantile(allDiseaseSNPfst_sum, probs=c(.99,.95,.90))
diseaseFst_popsVsall_quantiles <- quantile(diseaseFst_popsVsall_sum, probs=c(.99,.95,.90))

# top 1%, 5%, 10% tables for all subsets
allSNPfst_topTables <- lapply(allSNP_quantiles, function(x){ return(allSNPfst_sum[ allSNPfst_sum > x])}) # theyre named properly too!

allSNPfst_popsVSall_topTables <- lapply(allSNP_quantiles, function(x){
                                        return(allSNPfst_sum[ allSNPfst_sum > x])})

allDiseaseSNPfst_topTables <- lapply(allDiseaseSNPfst_quantiles, function(x){
                                        return(allDiseaseSNPfst_sum[ allDiseaseSNPfst_sum > x])})

diseaseFst_popsVsalltopTables <- lapply(diseaseFst_popsVsall_quantiles, function(x){
                                        return(diseaseFst_popsVsall_sum[ diseaseFst_popsVsall_sum > x])})


# generate list of quantiles for lists of Fst by population:

allFst_allPops_list_quantiles <- lapply(allFst_allPops_list_sum, function(x){ return(quantile(x, probs=c(.99,.95,.90)))})
diseaseFst_allPops_list_quantiles <- lapply(diseaseFst_allPops_list_sum, function(x){ return(quantile(x, probs=c(.99,.95,.90)))})


# top tables for list of per pop subsets

allFst_allPops_list_quantiles[[1]][1]
length(allFst_allPops_list_quantiles[[1]])


head(allFst_allPops_list_sum[[1]])

names(allFst_allPops_list_sum[1])


allFst_allPops_list_topTables <- list()
tempTableList <- list()
for(i in 1:length(allFst_allPops_list_quantiles)){

  tempTableList <- list("99%" = allFst_allPops_list_sum[[i]][allFst_allPops_list_sum[[i]] > allFst_allPops_list_quantiles[[i]][1]],
                        "95%" = allFst_allPops_list_sum[[i]][allFst_allPops_list_sum[[i]] > allFst_allPops_list_quantiles[[i]][2]],
                        "90%" = allFst_allPops_list_sum[[i]][allFst_allPops_list_sum[[i]] > allFst_allPops_list_quantiles[[i]][3]])

  allFst_allPops_list_topTables[[names(allFst_allPops_list_sum[i])]] <- tempTableList
}


diseaseFst_allPops_list_topTables <- list()
tempTableList <- list()
for(i in 1:length(diseaseFst_allPops_list_quantiles)){

  tempTableList <- list("99%" = diseaseFst_allPops_list_sum[[i]][diseaseFst_allPops_list_sum[[i]] > diseaseFst_allPops_list_quantiles[[i]][1]],
                        "95%" = diseaseFst_allPops_list_sum[[i]][diseaseFst_allPops_list_sum[[i]] > diseaseFst_allPops_list_quantiles[[i]][2]],
                        "90%" = diseaseFst_allPops_list_sum[[i]][diseaseFst_allPops_list_sum[[i]] > diseaseFst_allPops_list_quantiles[[i]][3]])

  diseaseFst_allPops_list_topTables[[names(diseaseFst_allPops_list_sum[i])]] <- tempTableList
}


# Establish important data to asssociate with top tables - (associated by 'VariantID')
colnames(full_SNP_Annotations_GWASc_Ensembl) # honestly looking at the cols reveals that it may just be best to do this sort of association and reporting of genes or other interesting links once it comes up in an actual report for the project.


# save all top tables


allSNPfst_topTables[[1]][1:10]
allSNPfst_popsVSall_topTables[[1]][1:2]

allDiseaseSNPfst_topTables[[1]][1:10]
diseaseFst_popsVsalltopTables[[1]][1:10]

allFst_allPops_list_topTables
diseaseFst_allPops_list_topTables

# everything looks good with the tables.. all named^^^

getwd()
setwd('..')
setwd('./fst_sum_topTables/')

save(allSNPfst_topTables, file = "allSNPfst_topTables.rds")
save(allSNPfst_popsVSall_topTables, file = "allSNPfst_popsVSall_topTables.rds")
save(allDiseaseSNPfst_topTables, file = "allDiseaseSNPfst_topTables.rds")
save(diseaseFst_popsVsalltopTables, file = "diseaseFst_popsVsalltopTables.rds")
save(allFst_allPops_list_topTables, file = "allFst_allPops_list_topTables.rds")
save(diseaseFst_allPops_list_topTables, file = "diseaseFst_allPops_list_topTables.rds")


# top 1, 5 10% tables all saved for all variations... Should do a comparison of the ALL population table when processed in the list vs when processed individually before proceeding to analyses to ensure things are consistent...
#
# MOVE INTO VISUAL ANALYSIS AND REPORTING ON PRIMARY QUESTIONS NEXT SESSION!! (3-22-23)

# replicating previous analyses:

load("./workingData/fst_sums/allSNPfst_sum.rds") # its a numeric object! (maybe why I couldn't find it at first)


q99 <- quantile(allSNPfst_sum, 0.99)
q95 <- quantile(allSNPfst_sum, 0.95)
q90 <- quantile(allSNPfst_sum, 0.90)

# pch = point character, cex = character expand
plot(allSNPfst_sum, main = "All Fst Sums, All Population Pairs", ylab = "Sum of Variant Fst", xlab = "Index of Sample", pch = 20, cex = .7)

abline(h = q99, col = 'green', lwd = 2)
abline(h = q95, col = 'blue', lwd = 2)
abline(h = q90, col = 'red', lwd = 2)

text(-4000, y = q99, "1%", pos = 3, offset = .5, col = 'darkgreen')
text(-4000, y = q95, "5%", pos = 3, offset = .5, col = 'darkblue')
text(-4000, y = q90, "10%", pos = 3, offset = .5, col = 'darkred')
# Beautiful

dens_allSNPfst <- density(allSNPfst_sum)
plot(dens_allSNPfst, main = "Density of All SNPs All Pop-Pairs")

summary(allSNPfst_sum)
sd(allSNPfst_sum)

abline(v = q99, col = 'green', lwd = 2)
abline(v = q95, col = 'blue', lwd = 2)
abline(v = q90, col = 'red', lwd = 2)
sumSts <- c("Mean: 27.44", "Median: 20.88", "Standard Deviation: 24.25", "Min: 0.00", "Max: 197.02")
legend("topright", inset = 0.02, legend = sumSts, bg = "white", cex = .80, text.font = 4, box.lwd = 0)

# Now lets do the same for ALL disease SNPs only



q99 <- quantile(allDiseaseSNPfst_sum, 0.99)
q95 <- quantile(allDiseaseSNPfst_sum, 0.95)
q90 <- quantile(allDiseaseSNPfst_sum, 0.90)

# pch = point character, cex = character expand
plot(allDiseaseSNPfst_sum, main = "All Disease Fst Sums w/ All Population Pairs", ylab = "Sum of Variant Fst", xlab = "Index of Sample", pch = 20, cex = .7)

abline(h = q99, col = 'green', lwd = 2)
abline(h = q95, col = 'blue', lwd = 2)
abline(h = q90, col = 'red', lwd = 2)

text(-4000, y = q99, "1%", pos = 3, offset = .5, col = 'darkgreen')
text(-4000, y = q95, "5%", pos = 3, offset = .5, col = 'darkblue')
text(-4000, y = q90, "10%", pos = 3, offset = .5, col = 'darkred')
# Beautiful

dens_allDiseaseSNPfst <- density(allDiseaseSNPfst_sum)
plot(dens_allDiseaseSNPfst, main = "Density of All Disease SNPs w/ All Pop-Pairs")

summary(allDiseaseSNPfst_sum)
sd(allDiseaseSNPfst_sum)

abline(v = q99, col = 'green', lwd = 2)
abline(v = q95, col = 'blue', lwd = 2)
abline(v = q90, col = 'red', lwd = 2)
sumSts <- c("Mean: 28.08", "Median: 21.68", "Standard Deviation: 23.96", "Min: 0.00", "Max: 165.72")
legend("topright", inset = 0.02, legend = sumSts, bg = "white", cex = .80, text.font = 4, box.lwd = 0)
# they look practically identical. Which is great! Meaning the slice of disease data does have some standout SNPs in terms of Fst

# lets take a look at a few poppulations against the meta population to get a sense of whether or not any disease-SNPs stick out without an Fst Sum being taken.. basically the question here is: are we seeing any SNPs which obviously stand out for any populations... (maybe if we set like a threshold value we can quickly scan the whole DF)

diseaseFst_popsVSall
dfpvInspect <- dplyr::copy(diseaseFst_popsVSall) #not a dplyr func.. there is copy for data.table though
dfpvInspect <- diseaseFst_popsVSall


?copy
?deepcopy #nothing
?clone #nothgin

# gonna have to figure out how shallow copies work in R...

dfpvInspect[dfpvInspect < .8] <- NA
class(dfpvInspect)
sum(is.na.data.frame(dfpvInspect)) # 1,119,007
sum(is.na.data.frame(diseaseFst_popsVSall)) # 0

# ok so somehow R knows to make a deep copy when modifying the shallow copy.. or has some way to work around the usual constraints of a shallow copy.

dim(dfpvInspect) # 36097    31
36097 * 31 # 1119007 ... means no values above .8 against the total population

quantile(as.matrix(diseaseFst_popsVSall), probs = c(.999, .99, .95, .90, .75, .5, .25))
?quantile

summary_diseaseFst_popsVSall <- summary(as.matrix(diseaseFst_popsVSall))
summary_diseaseFst_popsVSall # we can see just a handfull of values breach .5 but there are some real standouts per population.






