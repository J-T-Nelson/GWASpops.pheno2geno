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




# Data digest format / function development -------------------------------

#  Need to get my current data structures efficiently queried or transformed for processing via HudsonFst()
#
#  I think using the perAllele tables makes sense, as I can essentially create the pairwise combinations as rows which can than easily be iterated over by using an apply(HudsonFst()) sort of call. Hit ever row to generate value for a HudsonValue row essentially
#
#  At some point I need to be able to see PER allele PER pairwise population combination Fst values I think. So, it would make sense to just digest current data frames into a dataFrame which has all pairs represented as their own rows, then iterating rowise with the HudsonFst() func to create a new column that is then added to the DF by a wrapper func (wrapper to the perAlleleFstTransform() func)

# what should each row look like? i.e., what should my cols be?
#     Pair, n1, p1, n2, p2, Fst_Hudson <- hudson col should only be there after calculated.. transform won't add it right away.
#
# EX: 1000GENOMES:phase_3:AFR-X-1000GENOMES:phase_3:AMR | 555 | .3 | 241 | .45 | .4145


thousGenPops <- populations[grep("1000GENOMES",populations$Population_Abbreviation)]
# Num pairs = 32 Choose 2
numPairs <- choose(32,2)
numPairs # 496

# so the goal is to digest a list of DFs and to produce a list of DFs. Each DF should contain the 'Ancestral_Allele' and 'VariantID' attributes still
# Each DF will be made, assigned attributes, the preDigest frame will have irrelevant values removed, then will be digested for relevant values, the value

AA <- attr(perAllele_1, 'Ancestral_Allele')

perAllele_1$population %in% thousGenPops$Population_Abbreviation
perAllele_POI_t <- perAllele_1[perAllele_1$population %in% thousGenPops$Population_Abbreviation, ]
perAllele_POI <- perAllele_1[perAllele_1$population %in% thousGenPops$Population_Abbreviation & perAllele_1$allele != AA, ]


n1ex <- thousGenPops[thousGenPops$Population_Abbreviation == perAllele_POI$population[1]]$Sample_Count # not sure how efficient this is.. but I may need to find a less computationally intense lookup method to get the right number in the right place. ...
# 1. could just put them into alleleDF off the bat then I just grab based on index... Not sure if more efficient.
# 2. Also could sort the populations DF to match the indices of the alleleDF s.t. the i and j values can be used to reference by position without comparison based lookup.

# making a dummy row
dummyRow <- c(paste0(perAllele_POI$population[1], "-X-",perAllele_POI$population[2]),
              thousGenPops[thousGenPops$Population_Abbreviation == perAllele_POI$population[1]]$Sample_Count,
              perAllele_POI$frequency[1],
              thousGenPops[thousGenPops$Population_Abbreviation == perAllele_POI$population[2]]$Sample_Count,
              perAllele_POI$frequency[2])
dummyRow
# [1] "1000GENOMES:phase_3:ACB-X-1000GENOMES:phase_3:AFR" "96"
# [3] "0.135416666666667"                                 "661"
# [5] "0.119515885022693"

# ^^ looks good. Now the question is how computationally expensive is this compared to just sorting both DFs and doing no comparisons when grabbing the data? I suppose I don't have a good sense of this, and will probably be best off just rolling with the ugly method first, and getting more tricky adhoc if performance is a problem. !!!!




# perAllele... will take in a single allele DF and the populations desired to be digested, and will return a single DF
#   `populations` should contain the sample # values as well as the populations names (abbreviations) so it will be a filtered DF from the `Populations` obj of my GWASp2g package

perAlleleFst_transform <- function(alleleDF, populations){

  # Extract data of interest from alleleDF
  ancestralAllele <- attr(alleleDF, "Ancestral_Allele")
  alleleDF <- alleleDF[alleleDF$population %in% populations$Population_Abbreviation & alleleDF$allele != ancestralAllele , ]

  # Digest DF in to create DF out:

  DF_rows <- nrow(alleleDF) # number for efficient pairwise iteration
  rowHolder <- list()
  for(i in 1:(DF_rows-1)){ # i correlates to a population
    for(j in (1+i):DF_rows){ # j correlates to the second population used to pair with i's population
      rName <- paste0(alleleDF$population[i], "-X-",alleleDF$population[j])
      row <- c(rName,
               populations[populations$Population_Abbreviation == alleleDF$population[i]]$Sample_Count,
               alleleDF$frequency[i],
               populations[populations$Population_Abbreviation == alleleDF$population[j]]$Sample_Count,
               alleleDF$frequency[j])

      rowHolder[[rName]] <- row
    }
  }

  #setup return DF, name cols, fix col types, assign attributes
  retDF <- as.data.frame(t(dplyr::bind_rows(rowHolder)))

  names(retDF) <- c("pop_pair", 'n1', "p1", 'n2', 'p2')
  retDF['n1'] <- as.numeric(retDF[['n1']])
  retDF['n2'] <- as.numeric(retDF[['n2']])
  retDF['p1'] <- as.numeric(retDF[['p1']])
  retDF['p2'] <- as.numeric(retDF[['p2']])

  attr(retDF, "Ancestral_Allele") <- ancestralAllele
  attr(retDF, "VariantID") <- attr(alleleDF, "VariantID")

  return(retDF)
}



tDF <- data.frame()
attr(tDF, "Ancestral_Allele") <- attr(perAllele_1, "Ancestral_Allele")
attributes(tDF)

length(perAllele_1[,1])

rh <- list()
append(rh, 'g')
append(rh, 'good') # not working for appending new elements to list
rh
r1 <- c(rh, 'g', 'gg', 5)
r2 <- c(rh, 'good', 5 , 6)
rh <- c(rh, r1)
rh <- c(rh, r2)
rh <- c(rh, dummyRow) # its appending things as individual object which is undesirable.... maybe if I convert to a DF?
rh # works well

rh <- list()
r1 <- list( 'g', 'gg', 5)
r2 <- list('good', 5 , 6)
rh <- list( r1, r2)
rh <- list(rh, r1, r2) # this is working to put each item listed as its own entry in the list without properly flat structure.. so would induce highly nested structure if called within a for loop like I am doing.

?append # appears to be for vectors alone, not for lists.



rh <- list()
rh[["ham1"]] <- dummyRow
rh[["hammm"]] <- dummyRow
rh[["hammgf"]] <- r1
rh


# thinking that a named list is ideal... then each row will carry the name with it.. we still need to name cols though.
rhdf <- dplyr::bind_rows(rh)
rhdf_col <- dplyr::bind_cols(rh)
rhdfFinal <- as.data.frame(t(rhdf)) # its a fucking character vec
names(rhdfFinal) <- c("pop_pair", 'n1', "p1", 'n2', 'p2')
rhdfFinal
str(rhdfFinal) # bunch of character vecs
rhdfFinal['n1'] <- as.numeric(rhdfFinal[['n1']])
rhdfFinal['n2'] <- as.numeric(rhdfFinal[['n2']])
rhdfFinal['p1'] <- as.numeric(rhdfFinal[['p1']])
rhdfFinal['p2'] <- as.numeric(rhdfFinal[['p2']])
str(rhdfFinal) # looks good now.




?bind_rows
# all pairwise relationships should be made by this.. we expect to see 496 for only 1kgenomes stuff.




# TEST perAlleleFst_transform()
#



hudmoPrep <- perAlleleFst_transform(perAllele_1, thousGenPops) # success. debugged it pretty quick. Had to adjust loop amounts just a tad

debug(perAlleleFst_transform)
undebug(perAlleleFst_transform)

# Now we should test using HudsonFst() on the thing to generate a new row..

HudsonFst <- function(n1, n2, p1, p2){
  return( ((p1-p2)**2 - (p1*(1-p1))/(n1-1) - (p2*(1-p2))/(n2-1)) / (p1*(1-p2) + p2*(1-p1)) )
}

# I am realizing that my current method of storing n1, p1, n2, p2, is very memory inefficient. as they are stored MANY times as I do more SNPs and truthfully they can likely be looked up adhoc when calculating HudsonFst, and we only really need HudsonFst to be associated with the pairwise comparison... We will have to consider making a more memory efficient scheme in the future possibly as I will have like... tens of thousands of redundant data entries in this current scheme, which may bog down the feasability of producing meaningful graphs. Hm.
# ^^^ I could just store the values temporarily, then calculate Fst and delete them or store them in whatever way is useful for visual exploration later on...
#
#
# Need to brainstorm how I would visually explore things tomorrow morning.




# Developing options for perAlleleFst_transform ---------------------------

perAlleleFst_transform <- function(alleleDF, populations, deleteRedundants = FALSE){

  # Extract data of interest from alleleDF
  ancestralAllele <- attr(alleleDF, "Ancestral_Allele")
  alleleDF <- alleleDF[alleleDF$population %in% populations$Population_Abbreviation & alleleDF$allele != ancestralAllele , ]

  # Digest DF in to create DF out:

  DF_rows <- nrow(alleleDF) # number for efficient pairwise iteration
  rowHolder <- list()
  for(i in 1:(DF_rows-1)){ # i correlates to a population
    for(j in (1+i):DF_rows){ # j correlates to the second population used to pair with i's population
      rName <- paste0(alleleDF$population[i], "-X-",alleleDF$population[j])
      row <- c(rName,
               populations[populations$Population_Abbreviation == alleleDF$population[i]]$Sample_Count,
               alleleDF$frequency[i],
               populations[populations$Population_Abbreviation == alleleDF$population[j]]$Sample_Count,
               alleleDF$frequency[j])

      rowHolder[[rName]] <- row
    }
  }

  #setup return DF, name cols, fix col types, assign attributes
  retDF <- as.data.frame(t(dplyr::bind_rows(rowHolder)))

  names(retDF) <- c("pop_pair", 'n1', "p1", 'n2', 'p2')
  retDF['n1'] <- as.numeric(retDF[['n1']])
  retDF['n2'] <- as.numeric(retDF[['n2']])
  retDF['p1'] <- as.numeric(retDF[['p1']])
  retDF['p2'] <- as.numeric(retDF[['p2']])

  attr(retDF, "Ancestral_Allele") <- ancestralAllele
  attr(retDF, "VariantID") <- attr(alleleDF, "VariantID")

  fstVec <- numeric(nrow(retDF))
  for(i in 1:nrow(retDF)){
    fstVec[i] <- HudsonFst(retDF[i,2],retDF[i,4],retDF[i,3],retDF[i,5])
  }
  retDF['Fst_Hudson'] <- fstVec

  if(deleteRedundants){ #removing all but population pairs and Fst value to save memory
    retDF <- retDF[,c(1,6)]
  }

  return(retDF)
}


# Test calculating Fst against new table.

fstCol <- sapply(hudmoPrep, \(x) HudsonFst(x$n1, x$n2, x$p1, x$p2))
fstCol <- apply(hudmoPrep, 1, \(x) HudsonFst(x[,2], x[,4], x[,3], x[,5])) # neither type of referencing is working properly... I wonder if there is a way to get this to go... chatGPT gave me a suggestion to use a rapper that changes the means of referencing inherently so I will try that

hudmoPrep[1,2]
hudmoPrep[1,1]
class(hudmoPrep[1,1])
rowSubset <- hudmoPrep[1,]
str(rowSubset)
rowSubset[2] # names of col and row are attached



#   Positional arguments for sake of apply usage.. ERROR... doesn't work
hFst_applyWrapper <- function(alleleDF_row) {
  return(HudsonFst(alleleDF_row[2], alleleDF_row[4], alleleDF_row[3], alleleDF_row[5]))
}



fstCol <- apply(hudmoPrep, 1, hFst_applyWrapper)
?apply

tret <- hFst_applyWrapper(hudmoPrep[1,])
class(tret)
tHudsFST <- HudsonFst(hudmoPrep[1,2],hudmoPrep[1,4],hudmoPrep[1,3],hudmoPrep[1,5])
tHudsFST # properly returning just a number.. however it is a negative number which makes me wonder if my function has been written correctly


# trying to use apply is wasting time... lets use  for loop

fstVec <- numeric(nrow(hudmoPrep))
for(i in 1:nrow(hudmoPrep)){
  fstVec[i] <- HudsonFst(hudmoPrep[i,2],hudmoPrep[i,4],hudmoPrep[i,3],hudmoPrep[i,5])
}

?append


# TEST hudson FST addition ------------------------------------------------

debug(perAlleleFst_transform)
hudsonTest2 <- perAlleleFst_transform(perAllele_1, thousGenPops)

# test dropping cols

less <- hudsonTest2[,c(1,6)]



# Develop Ancestral_Allele filling ----------------------------------------

testCase <- testAlleleList[[5]]

Avals <- sum(testCase[testCase$allele == "A", ]$frequency)
Cvals <- sum(testCase[testCase$allele == "C", ]$frequency)
Gvals <- sum(testCase[testCase$allele == "G", ]$frequency)
Tvals <- sum(testCase[testCase$allele == "T", ]$frequency)

which.max(c(Avals, Cvals, Gvals, Tvals))

switchTest <- switch(which.max(c(Avals, Cvals, Gvals, Tvals)),
                     '1' = "A", '2' = "C", '3' = "G", '4' = "T")

calc_ancestralAllele <- function(population_alleleDF){

  A_mag <- sum(population_alleleDF[population_alleleDF$allele == "A", ]$frequency)
  C_mag <- sum(population_alleleDF[population_alleleDF$allele == "C", ]$frequency)
  G_mag <- sum(population_alleleDF[population_alleleDF$allele == "G", ]$frequency)
  T_mag <- sum(population_alleleDF[population_alleleDF$allele == "T", ]$frequency)
  tempVec <- c(A_mag, C_mag, G_mag, T_mag)

  ancestralAllele <- switch(which.max(tempVec),
                            '1' = "A",
                            '2' = "C",
                            '3' = "G",
                            '4' = "T")
  return(ancestralAllele)
}

tCalc <- calc_ancestralAllele(testCase) # works

# -------------------------------------------------------------------------



# Developing: Wrapper to process list of perAllele tables  ----------------

testAlleleList <- tml[[2]]
names(testAlleleList[1])
# list of DFs in, list of DFs out.
testAlleleList[1]
testAlleleList[[1]]


hudsonFst_alleleList <- function(alleleList, populationsDF, deleteRedundants = FALSE){

  captureList <- list()

  for(i in 1:length(alleleList)){
    tableName <- names(alleleList[i])
    captureList[[tableName]] <- perAlleleFst_transform(alleleList[[i]], populationsDF, deleteRedundants)
  }

  return(captureList)
}

# TESTING hudsonFst_alleleLIst

hudTest3 <- perAlleleFst_transform(testAlleleList[[2]], thousGenPops)
hudTest4 <- perAlleleFst_transform(testAlleleList[[3]], thousGenPops) # both work, variable numbers of obs. out (rows) .. must just be related to a variable number of valid pairs... so far all outputs are odd values.. which is curious.. I wonder if there is any meaning behind that.



debug(hudsonFst_alleleList)
undebug(hudsonFst_alleleList)
tHAL <- hudsonFst_alleleList(testAlleleList, thousGenPops) # after solving the two issues noted just below this executed relatively quickly at around 10 seconds?
tHAL <- system.time(hudsonFst_alleleList(testAlleleList, thousGenPops))
tHAL
# user  system elapsed     so ... about 16 seconds to complete on ~200 variants for the 1000k Genomes.. this seems like a fast enough rate for future processing.
# 15.04    0.24   15.53

debug(perAlleleFst_transform)
undebug(perAlleleFst_transform)

# ERROR:1
tPAF <- perAlleleFst_transform(testAlleleList[[14]], thousGenPops)
# for this one, we are actually seeing that there are no populations of interest (thousand genomes populations), and thus there are some interesting issues cropping up.. will have to check for this case and discard tables which fail to have any data of interest.

# ^^ returning NULL works standalone... lest see if it works in the wrapper .. it should. FIXED


# ERROR:2
tPAF <- perAlleleFst_transform(testAlleleList[[5]], thousGenPops)
# we are seeing an issue where ancestral alleles can be NA, which totally ruins the function. Best quick solution I can imagine is to sum the various alleles frequencies and to assign the highest one as the major (ancestral) allele...

# FIXED. lets find the next bug lol












# -------------------------------------------------------------------------

head(tHAL[[1]], 2)
tHAL[[1]][1:2, c(2:6)]























