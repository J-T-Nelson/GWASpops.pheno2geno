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

names(diseaseFst_allPops_list_quantiles) <- names(diseaseFst_allPops_list_sum)

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

L_diseaseFst_popsVSall <- split(diseaseFst_popsVSall, colnames(diseaseFst_popsVSall)) # split it row-wise

L_diseaseFst_popsVSall <- lapply(diseaseFst_popsVSall, as.list) # this does it.
## lapply(1:length(L_diseaseFst_popsVSall), function(x){ names(L_diseaseFst_popsVSall[[x]]) <- rownames(diseaseFst_popsVSall)})  # lapply just not a full replacement to 'for()' loops.
##
for(i in 1:length(L_diseaseFst_popsVSall)) {names(L_diseaseFst_popsVSall[[i]]) <- rownames(diseaseFst_popsVSall)} # worked


L_diseaseFst_popsVSall_topSnps <- lapply(L_diseaseFst_popsVSall, function(x){ return(x[x > .4])})


n <- names(L_diseaseFst_popsVSall_topSnps)
n <- gsub('1000GENOMES:phase_3:ALL', '', n)
n <- gsub('-X-', '', n)


k <- sapply(L_diseaseFst_popsVSall_topSnps, length)
topSnps_popsVSall <- cbind(n,k)
tsp <- merge(topSnps_popsVSall, pops, by.x = 'n', by.y = 'Population_Abbreviation')
topSnps_popsVSall <- tsp


# naming hard...
names(topSnps_popsVSall$n) <- 'population_abbreviation' # wrong.

colnames(topSnps_popsVSall)[colnames(topSnps_popsVSall)== 'n'] <- 'population_abbreviation'
colnames(topSnps_popsVSall)[colnames(topSnps_popsVSall)== 'k'] <- 'numSnps_>_0.4'


colnames(topSnps_popsVSall) <- c("population_abbreviation", 'numSnps_>_0.4', "Sample_Count","Pop_Ancestry","PopAnces_Graph_Labels") # effective for multiple renames .. grab names.. reenter what you want, rename everything

topSnps_popsVSall <- topSnps_popsVSall[order(as.numeric(topSnps_popsVSall$`numSnps_>_0.4`), decreasing = T), ]

sort(topSnps_popsVSall$`numSnps_>_0.4`)

# lowering threshold to see how many SNPs we capture.. lets try .35
L_diseaseFst_popsVSall_topSnps_35 <- lapply(L_diseaseFst_popsVSall, function(x){ return(x[x > .35])})

# considering how to find degree of intersection for each SNP and unique SNPs for each population...
#
# could find uniques by intersecting each pop with each other, and then compiling a novel vector containing all SNPs found in an intersection,
# then by taking those out of each of our vectors we would reveal all SNPs unqiue to a single population
#
# The issue of getting the degree of intersection per snp still exists however... I suppose the way to find this would then be to take this list of intersection positive SNPs and to look for each one across all populations and to count the number of times it comes up.. could additonally mark which populations we see it in within a separate vector.. thus capturing the information we want.


# testing some Gpt4 code for the job :

n <- names(L_diseaseFst_popsVSall_topSnps)
n[1]
n <- gsub("1000GENOMES:phase_3:ALL",'',n)
n <- gsub("-X-",'',n)
n[1:10]

names(L_diseaseFst_popsVSall_topSnps) <- n
names(L_diseaseFst_popsVSall_topSnps) # good


# Generate all possible pairs of indices
index_pairs <- combn(length(L_diseaseFst_popsVSall_topSnps), 2, simplify = FALSE)
?combn
index_pairs

intersection_results <- purrr::map(index_pairs, ~dplyr::intersect(L_diseaseFst_popsVSall_topSnps[[.x[1]]],
                                                                  L_diseaseFst_popsVSall_topSnps[[.x[2]]]))


n1 <- names(L_diseaseFst_popsVSall_topSnps[[30]])
n2 <- names(L_diseaseFst_popsVSall_topSnps[[31]])
n1_2intersect <- dplyr::intersect(n1,n2)
n1_2intersect
# looking at the intersection of names to verify we are getting real matches.. as its currently the fst values being compared. Real matches indeed.

#naming the intersections:
pop_n <- gsub("1000GENOMES:phase_3:",'',names(L_diseaseFst_popsVSall_topSnps))

intersect_names <- purrr::map(index_pairs, ~paste0(pop_n[.x[1]], "_", pop_n[.x[2]]))
?purrr::map
intersect_names[1:10]

names(intersection_results) <- intersect_names
names(intersection_results) # good

# Removing intersections with value of 0, seeking multiple intersections:


intResultsFiltered <- intersection_results[ sapply(intersection_results, length) > 0 ]

sum(sapply(intersection_results, length)) # 48 ... so 48 SNPs which have SOME intersection across all pop-pairs

# checking for higher degree of intersection
allIntersections <- purrr::flatten(intResultsFiltered)
length(allIntersections) # 48
length(unique(allIntersections)) # 22 .. so we see several multi-intersectional SNPs.. Lets count up how many times we see them
allIntersections <- names(allIntersections) # grabbing just the SNP names now
allIntTable <- table(allIntersections)
allIntTable
#rs10873298 rs10873299 rs11114149 rs11692588  rs12075 rs12297948 rs12479436 rs13003464   rs143384  rs2058619   rs228768  rs2643826 rs35085068
# 1          1          3          3         10          3          3          1          1          1          1          6          1
# rs3741353  rs4553272  rs4565870  rs4842266 rs56335113 rs61826828  rs7575465  rs8181996   rs943451
# 3          1          1          3          1          1          1          1          1

# ^^ and above we have identified the degree of intersecting SNPs.. one is totally dominant at 10 shared groups.. while the next highest is at 6. Others are mostly 3 or 1.. and the majority of SNPs generally are uniquely high for one population against the metapopulation... lets see that count and enumerate those?


sum(sapply(L_diseaseFst_popsVSall_topSnps, length)) # 1207 for SNPs greater than .4

sum(sapply(L_diseaseFst_popsVSall_topSnps_35, length)) # 2746 for SNPs greater than .35


# 3-26 start:
#
uniIntersectingVariants <- unique(allIntersections)
uniIntersectingVariants

names(L_diseaseFst_popsVSall_topSnps) <- pop_n

# check for each term in each pop, if pop contains, paste pop name to vector.

intVariantPopulations <- list()
for(i in 1:length(uniIntersectingVariants)){
  temp <- character()
  for(j in 1:length(L_diseaseFst_popsVSall_topSnps)){
    if(uniIntersectingVariants[i] %in% names(L_diseaseFst_popsVSall_topSnps[[j]])){
      temp <- paste0(temp,"_", names(L_diseaseFst_popsVSall_topSnps)[j])
    }

  }
  intVariantPopulations[[ uniIntersectingVariants[i] ]] <- temp
}

uniIntersectingVariants[1]

names(L_diseaseFst_popsVSall_topSnps[[2]])

# refining names for readability
names(diseaseFst_allPops_list_quantiles) <- gsub("1000GENOMES:phase_3:",'',names(diseaseFst_allPops_list_quantiles))

# extract top 1% of each, then compare to top SNPs against the metapopulation

disFst_ap_top1percent <- lapply(diseaseFst_allPops_list_quantiles, function(x){ x[1]}) # looks like my initial grab of the data wasn't any good actually.. lets just grab the top 1% now then.

# Grabbing top .5% of values as top 1% is 361 values per population which may be too inclusive for comparison to top values found in comparing to metapopulation
disFst_ap_top1percent <- lapply(diseaseFst_allPops_list_sum, function(x){
  tempQuant <- quantile(x, .995);
  return(x[x>tempQuant])})

names(disFst_ap_top1percent) <- gsub("1000GENOMES:phase_3:",'',names(disFst_ap_top1percent))

disFst_ap_top1percent <- disFst_ap_top1percent[-1]

comp_top1_above.4 <- lapply(1:length(L_diseaseFst_popsVSall_topSnps), function(x){
    names(L_diseaseFst_popsVSall_topSnps[[x]]) %in% names(disFst_ap_top1percent[[names(L_diseaseFst_popsVSall_topSnps[x])]])
  })

names(comp_top1_above.4) <- names(disFst_ap_top1percent)

# revealing any values not in top .5% compared to above .4 vs all
sapply(comp_top1_above.4, sum)
sapply(comp_top1_above.4, length)

# naming truth vector
for(i in 1:length(comp_top1_above.4)){
  names(comp_top1_above.4[[i]]) <- names(L_diseaseFst_popsVSall_topSnps[[i]])
}

#get variant IDs of interest. Against big table, grab rows and cols of interest. Then add to it with the Fst and further data on SNPs of interest... 2 tables for the two sets of SNPs

names(comp_top1_above.4)
above.4_variants <- purrr::flatten(lapply(L_diseaseFst_popsVSall_topSnps, function(x){names(x)}))
above.4_variants <- as.character(above.4_variants)

top.5percent_variants <- as.character(purrr::flatten(lapply(disFst_ap_top1percent, function(x){names(x)})))

length(above.4_variants) # 1207
length(top.5percent_variants) # 5611
length(unique(above.4_variants)) # 477
length(unique(top.5percent_variants)) # 1876
# not nearly as many unique as I would have thought... hm.... specifically in the above.4_variants I wouldn't expect to see that many non unique names.. maybe our earlier intersection comparison was inadequate.

# going to redo the intersection analysis with names instead of values to be sure
above.4_names <- lapply(L_diseaseFst_popsVSall_topSnps, function(x){names(x)})

intersection_results <- purrr::map(index_pairs, ~dplyr::intersect(above.4_names[[.x[1]]],
                                                                  above.4_names[[.x[2]]]))


sum(sapply(purrr::flatten(intersection_results), length)) # 1677 intersections total.
sum(sapply( unique(purrr::flatten(intersection_results)) , length)) # 277 unique intersecting variants.

# lets get the new list of intersecting values and the degree of intersection now:

uniIntersectingVariants <- unique(as.character(purrr::flatten(intersection_results)))
length(uniIntersectingVariants) # 277
uniIntersectingVariants[1:20]


variantIntersectDegreeTable <- table(as.character(purrr::flatten(intersection_results)))
variantIntersectDegreeTable[1:10] # looks good.
variantIntersectDegreeTable
max(variantIntersectDegreeTable) # 45
summary(variantIntersectDegreeTable)
# rs10007784  rs10035291   rs1009840   rs1012621  rs10245867    rs103294  rs10551501  rs10714879  rs10761659  rs10772040  rs10858374  rs10873298
# 6          10           6           6           3           1           1           3           1          15           1          10
# rs10873299  rs10884064  rs10929757  rs10935182  rs10935185  rs11038927  rs11114149  rs11150602  rs11249906  rs11250135  rs11263955   rs1129038
# 10           1           1           1           1          15           6           1          10           3           1           3
# rs1151988  rs11590283  rs11611680  rs11673591  rs11675358  rs11680058  rs11692588 rs116945310  rs11814448  rs11845134  rs12040273   rs1205863
# 6           3           6          15           3          10           3          10          15           1           6          15
# rs12068879  rs12074934     rs12075  rs12171500  rs12263288  rs12297948  rs12425451  rs12479436  rs12548184  rs12621647  rs12913832  rs12959570
# 3          10          45           1           6           6          21           3           1           1           3           3
#######............. ...

sum(!(unique(above.4_variants) %in% uniIntersectingVariants)) # 200 variants which are unique to a given population.

# at this point all the groundwork is laid out to make sense of the data when reporting on it.. I don't think I will reasonably get a report done today without sacrificing time spent elsewhere.. thus I am going to just do 45 mins of report writing and focus on finishing it tomorrow morning.


singlePop_above.4_vars <- above.4_variants[!(unique(above.4_variants) %in% uniIntersectingVariants)]
length(unique(singlePop_above.4_vars))


above.4_names[1]












































