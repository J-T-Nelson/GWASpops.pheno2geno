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









































































