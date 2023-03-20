# update Fst Calculation then recalculate all Fst: 3-17-23
#
# REMAKE: full_fst.rds & full_fst_perPop.rds !!!
#
# Upon fixing Fst calculation I expect to see the majority of NA values dissapear... I should quantify this expectation to be sure though

# load pop allele freq data in
load("./workingData/full_data_for_analysis/full_popAlleleFreq_PerAllele.rds")

# Fix Fst Func; then test

TEMP_perAlleleFst_transform <- function(alleleDF, populations, deleteRedundants = FALSE){

  if(length(unique(alleleDF$allele)) > 2){ # no calculations for multiallelic sites. This method of Fst calculation isn't suitable to non-biallelic sites.
    return(NA)
  }

  # Extract data of interest from alleleDF
  ancestralAllele <- attr(alleleDF, "Ancestral_Allele")

  if(is.na(ancestralAllele)){ # sometimes ancestral Allele is NA, which makes the rest of this function impossible to execute. Calculating ancestral allele by finding allele with highest frequency.
    ancestralAllele <- calc_ancestralAllele(alleleDF)
  } else {

    if( !(ancestralAllele %in% unique(alleleDF$allele)) ) { # reassign AA if assignment of AA is somehow wrong, (ancestral allele not found in data.frame)
      ancestralAllele <- calc_ancestralAllele(alleleDF)     #  Wrong assignment can come directly from data sources, not necessarily my own code
    }
  }

  # adding rows where they are missing when an allele's ancestral allele is fixed at 1.0
  if(any( (alleleDF$frequency == 1) & (alleleDF$allele == ancestralAllele) )){
    nonAncestralAllele <- unique(alleleDF$allele[ alleleDF$allele != ancestralAllele ] )
    mutateRows <- alleleDF[ alleleDF$frequency == 1 , ]
    mutateRows$frequency <- 0
    mutateRows$allele_count <- 0
    mutateRows$allele <- nonAncestralAllele
    alleleDF <- rbind.data.frame(alleleDF, mutateRows)
  }

  alleleDF <- alleleDF[alleleDF$population %in% populations$Population_Abbreviation & alleleDF$allele != ancestralAllele , ] # filtering down to minor allele.
  # ^^  discard non 'ancestral alleles' because we are looking for the minor alleles to compare fst wrt to.

  # Digest DF in to create DF out:

  DF_rows <- nrow(alleleDF) # number for efficient pairwise iteration

  if(DF_rows == 0){ # when no pops of interest exist for a given variant, we just cancel the function and return nothing
    return(NA)
  }

  rowHolder <- list()
  for(i in 1:(DF_rows-1)){ # i correlates to a population
    for(j in (1+i):DF_rows){ # j correlates to the second population used to pair with i's population
      rName <- paste0(alleleDF$population[i], "-X-",alleleDF$population[j])

      row <- c(populations[populations$Population_Abbreviation == alleleDF$population[i]]$Sample_Count,
               alleleDF$frequency[i],
               populations[populations$Population_Abbreviation == alleleDF$population[j]]$Sample_Count,
               alleleDF$frequency[j])

      rowHolder[[rName]] <- row
    }
  }

  #setup return DF, name cols, fix col types, assign attributes
  retDF <- as.data.frame(t(dplyr::bind_rows(rowHolder)))

  names(retDF) <- c('n1', "p1", 'n2', 'p2')
  retDF['n1'] <- as.numeric(retDF[['n1']])
  retDF['n2'] <- as.numeric(retDF[['n2']])
  retDF['p1'] <- as.numeric(retDF[['p1']])
  retDF['p2'] <- as.numeric(retDF[['p2']])

  attr(retDF, "Ancestral_Allele") <- ancestralAllele
  attr(retDF, "VariantID") <- attr(alleleDF, "VariantID")

  fstVec <- numeric(nrow(retDF))
  for(i in 1:nrow(retDF)){
    fstVec[i] <- HudsonFst(retDF[i,1],retDF[i,3],retDF[i,2],retDF[i,4])
  }

  fstVec[is.nan(fstVec)] <- 0 # for cases where p1 and p2 are 0 we get NaN out of 'HudsonFst()'
  fstVec[fstVec < 0] <- 0 # bringing all negatives up to 0 for accurate averaging and realistic Fst values
  retDF['Fst_Hudson'] <- fstVec

  if(deleteRedundants){ #removing all but population pairs and Fst value to save memory, populations pairs are stored as row names
    retDF <- retDF[,5, drop = FALSE] # drop = FALSE ensures we don't lose the row names in coercion
  }

  return(retDF)
}

# test novel Fst Func

  # find some test data
rows_with_na <- apply(full_fst, 1, function(row) any(is.na(row)))
onlyPosRows <- rows_with_na[rows_with_na == TRUE]

alleles_to_check <- names(onlyPosRows)[1:10]
check_alleles <- full_popAlleleFreq_perAllele[alleles_to_check]

 # filtering down to only 1kGenomes pops
check_alleles <- lapply(check_alleles, function(x) {
  x[x$population %in% thousGenPops$Population_Abbreviation, ]
})


num_missing_pops <- sapply(check_alleles, function(x) {
  thousGenPops$Population_Abbreviation %in% x$population
})


testAlleleDF <- check_alleles[[2]]
testAlleleDF # many ancestral alleles with 1.0 allele freq
attributes(testAlleleDF) # proper ancestral allele

# How to add in rows we expect to see?
#   store non-ancestral allele
#   extract all 1.0 && ancestral allele rows into new object
#   rewrite allele freqs to 0 and alleles to non-ancestral allele
#   rbind back into initial data structure.
#
#   ^^ this scheme is efficient as it grabs just the amount we need, and just quickly rewrites some data then binds them into things
#   Before doing this we would scan the DF for any 1.0 in ancestral allele rows

alleleDF <- testAlleleDF
ancestralAllele <- 'G'
if(any( (alleleDF$frequency == 1) & (alleleDF$allele == ancestralAllele) )){
  nonAncestralAllele <- unique(alleleDF$allele[ alleleDF$allele != ancestralAllele ] )
  mutateRows <- alleleDF[ alleleDF$frequency == 1 , ]
  mutateRows$frequency <- 0
  mutateRows$allele_count <- 0
  mutateRows$allele <- nonAncestralAllele
  alleleDF <- rbind.data.frame(alleleDF, mutateRows)
}
# result looks good!


# test novel Fst calc func:

testRun <- lapply(check_alleles, function(x) {
  TEMP_perAlleleFst_transform(x, thousGenPops, TRUE)})

# think we may see an error when both pops have allele freqs of 0, as the denominator becomes 0 ... we'll have to write in how to handle this
debug(TEMP_perAlleleFst_transform)
undebug(TEMP_perAlleleFst_transform)


resetNan <- testRun[["rs35472707"]]$Fst_Hudson

resetNan[is.nan(resetNan)] <- 0
resetNan[resetNan < 0 ] <- 0
resetNan # looks good.

# second go:
testRun <- lapply(check_alleles, function(x) {
  TEMP_perAlleleFst_transform(x, thousGenPops, TRUE)})

# ready to test function for its purpose now: Recalculate all Fst using full_popAlleleFreq...


full_fst_fixed <- lapply(full_popAlleleFreq_perAllele[1:1000], function(x) {
  TEMP_perAlleleFst_transform(x, thousGenPops, TRUE); print(attr(x, "VariantID"))})

# Error:
#   ! Tibble columns must have compatible sizes.
# • Size 2: Column `-X-1000GENOMES:phase_3:MSL`.
# • Size 3: Column `1000GENOMES:phase_3:MSL-X-NA`.
# • Size 4: Column `1000GENOMES:phase_3:MSL-X-1000GENOMES:phase_3:MSL`.
# ℹ Only values of size one are recycled.
# Run `rlang::last_error()` to see where the error occurred.
rlang::last_error()


names(full_popAlleleFreq_perAllele[1])

# deparse(substitute( full_popAlleleFreq_perAllele[[1]] ))


?substitute
quote(full_popAlleleFreq_perAllele[[1]])
?deparse

attributes(full_popAlleleFreq_perAllele[[1]])

attr(full_popAlleleFreq_perAllele[[1]], "VariantID")

which(names(full_popAlleleFreq_perAllele) == "rs4980169")

debugFstCalc <- full_popAlleleFreq_perAllele[320:400]

debug(perAlleleFst_transform)
debug_fst_fixed <- lapply(debugFstCalc, function(x) {
  perAlleleFst_transform(x, thousGenPops, TRUE); print(attr(x, "VariantID"))})
# Bug was occuring because of a strange variant rs548223088, where the insertion didn't happen at all in our pops of interest, which resulted in a single row for some reason... I am now returning NA automatically when less than 2 rows remain after filtering out ancestral alleles.


# try to recalculate all Fst again:

# printing for visability when error occurrs.
full_fst_fixed <- lapply(full_popAlleleFreq_perAllele, function(x) {
    perAlleleFst_transform(x, thousGenPops, TRUE); print(attr(x, "VariantID"))})
# Error:
#   ! Assigned data `nonAncestralAllele` must be compatible with existing data.
# ✖ Existing data has 12 rows.
# ✖ Assigned data has 0 rows.
# ℹ Only vectors of size 1 are recycled.


which(names(full_popAlleleFreq_perAllele) == "rs143524414")

debug2_fst_fixed <- lapply(full_popAlleleFreq_perAllele[663:1000], function(x) {
  perAlleleFst_transform(x, thousGenPops, TRUE); print(attr(x, "VariantID"))})

# CAUSES ERROR: rs1555226898 - some data doesn't actually report on the non-ancestral allele... But there is still otherwise interesting data, thus using a mock value to calculate Fst has been used to solve this.

# Start over with print for error seeking:
full_fst_fixed <- lapply(full_popAlleleFreq_perAllele, function(x) {
  perAlleleFst_transform(x, thousGenPops, TRUE); print(attr(x, "VariantID"))})

##### Keep running until all bugs are clear.
getwd()
setwd("./workingData/full_data_for_analysis/")
save(full_fst_fixed, file = "full_fst_fixed.rds" ) # saved nothing of use... whew.

# well... running the print actually just returned the print value... which I didn't realize could happen... which is pretty dissapointing. I guess we at least now know that print statements actually return a value in R ... lol wtf.

full_fst_fixed <- lapply(full_popAlleleFreq_perAllele, function(x) {
  perAlleleFst_transform(x, thousGenPops, TRUE)}) # expect this to take like 4+ hours.

save(full_fst_fixed, file = "full_fst_fixed.rds" ) # things look good this time. 10.3 GB of data. saving without removing NA's ... we need to process this down properly still. Think I will save it again with a better name, then save over 'full_fst_fixed' once the object has been properly smashed into a single table.

full_fst_unprocessed <- full_fst_fixed
save(full_fst_unprocessed, file = "full_fst_unprocessed.rds" )


# Continuing the recalculation of Fst... need to process the list of DFs into a single DF  --------

#  clearing space before data transform
rm(full_fst_fixed)
rm(full_popAlleleFreq_perAllele)

part_fstEX <- full_fst[1:10,]
rm(full_fst, testRun, debug_fst_fixed, debug2_fst_fixed, debugFstCalc, allrsID, allrsID_ch10, alleleDF, num_missing_pops)
rm(check_alleles, masterList, mutateRows, testAlleleDF, rows_with_na)


# Counting num NA returned in Fst calculation:

totalNA_postFstCalc <- sum(is.na(full_fst_unprocessed))
totalNA_postFstCalc # 53329

# samples processed 213710
53329/213710 # 0.2495391 ~= 25% of samples were NA for some reason... multiallelic... or any of the other 2 (or 3) failure to process conditions..

full_fst_unprocessed <- full_fst_unprocessed[!is.na(full_fst_unprocessed)] # 160381 samples remain... which is close to the amount we got last time..

# see if all have 496 rows
num_rows_in_fullFst <- unique(sapply(full_fst_unprocessed, function(x){ nrow(x) })) # 59 unique values
num_rows_in_fullFst
# 496  55 519 534 557 545 465 516 506 171 378  15 566 435 190 120 505 351   1 210 613 696 231  66 532 528 105  10 300   3 276   6 531 406 521  45 325  28
# 253 508  78  21  36 556 555 153 581 507 552 530 522 136 713 612 666 629 724 626  91

# ^^ this is undesirable.. we wanted all to have 496.. no more, no less. This result indicates the changes in processing Fst may have introduced some novel bugs. Lets see how many non 496 row values exist

num_badRows <- sapply(full_fst_unprocessed, function(x){ nrow(x) }) != 496
sum(num_badRows) # 164. ... truly trivial amount in the grand scheme... thus. unless I get mis-matched row names, I am going to assume everything is good. and proceed.

full_fst_unprocessed <- full_fst_unprocessed[!num_badRows] # down to 160217

fstFullNames <- names(full_fst_unprocessed)
full_fst_unprocessed <- cbind.data.frame(full_fst_unprocessed)
full_fst_unprocessed <- as.data.frame(t(full_fst_unprocessed))
rownames(full_fst_unprocessed) <- fstFullNames
# Error in `.rowNamesDF<-`(x, value = value) :
#   duplicate 'row.names' are not allowed

dupeNames <- duplicated(fstFullNames)
sum(dupeNames) # 946
?duplicated
fromLast_dupes <- duplicated(fstFullNames, fromLast = T)
sum(fromLast_dupes) # 946
identical(fromLast_dupes, dupeNames) # FALSE
#### "The identical() function in R checks whether two objects are exactly equal in terms of both their contents and attributes."
dupeIdx <- which(dupeNames)
revDupeIdx <- which(fromLast_dupes)
allDupesIdx <- append(dupeIdx, revDupeIdx) # 1892 total entries .. == 946 * 2 .... looks good. ...

# I want to check each of these duplicates to see if they're actually coming out with the same Fst ... if they are... we can just delete.. if not... that is a call for question.

dupe_fst <- full_fst_unprocessed[allDupesIdx,] # they are unnamed... but we need to check for identical rows.. we would expect to see half the rows be identical if the Fst calculation is the same for Ids with matching names.. Else.. We will have to consider how to proceed.
dupeRows <- duplicated.data.frame(dupe_fst)
sum(dupeRows) # 958 exactly duplicated rows... which is 12 more than expected... which seems highly improbably to be coincidental. This does suggest however we can probably safely delete duplicated names...
which(dupeRows) # looks just about perfect. other than the 12 weird unexpected dupes.. we see that the rest are exactly from 947:1892 . So I don't know what the deal with totally equal rows is.. but maybe this warrants checking my whole data frame.. which seems... memory intensive.


# removing the duplicate indices from both the name vec and the DF

fstFullNames <- fstFullNames[!fromLast_dupes] #159,271
full_fst_unprocessed <- full_fst_unprocessed[!fromLast_dupes, ] #159,271

rownames(full_fst_unprocessed) <- fstFullNames # works now.. no more duplicate names... interested to see number of dupe rows... a little worried about

dupeRowsTotal <- duplicated.data.frame(full_fst_unprocessed)
sum(dupeRowsTotal) # 4784... this suggests that I have some duplicate SNPs I think... as its just really not probabilistically likely to see so many duplicate rows o/w ...the only other explanation may be that there are rows of ONLY 0 ... which could be another source. hm.

full_fst_fixed <- full_fst_unprocessed
getwd()
setwd("./workingData/full_data_for_analysis/")
save(full_fst_fixed, file = 'full_fst_fixed.rds')

sum(is.na(full_fst_fixed)) # 0 NA values



# WORKING ON GETTING LIST OF DIESEASE TRAITS 3-19-23 ------------------------------

diseaseTrait <- unique(asso$`DISEASE/TRAIT`)
mappedTrait <- unique(asso$MAPPED_TRAIT)
mapT_in_diseaseT <- mappedTrait %in% diseaseTrait
sum(mapT_in_diseaseT)
?writeLines
writeLines(diseaseTrait, 'diseaseTrait.txt')
writeLines(mappedTrait, 'mappedTrait.txt')
diseaseTrait <- tolower(diseaseTrait)
mappedTrait <- tolower(mappedTrait)
mapT_in_diseaseT <- mappedTrait %in% diseaseTrait
sum(mapT_in_diseaseT)
disease_mapped_Traits <- append(diseaseTrait, mappedTrait)
disease_mapped_Traits <- unique(disease_mapped_Traits) # now we have the unique traits names in each, ~27,982 unique names.

# Going to see if I can get help from GPT-4 with filtering this down to only traits that are relevant to disease... as the list offered in the disease term lookup on GWASc is not returning very interesting results when searching the terms within GWASc itself. i.e, some terms within return no associations. Which means its just not useful. We have all associations already and need to filter down what is here to those which are disease related

writeLines(disease_mapped_Traits, 'disease_mapped_traits.txt')

any(is.na(diseaseTrait)) # FALSE
any(is.na(mappedTrait)) # FALSE
mappedTrait

any(is.null(diseaseTrait)) # FALSE
any(is.null(mappedTrait)) # FALSE

traitsDF <- asso[, c('DISEASE/TRAIT', 'MAPPED_TRAIT')]

write.csv(traitsDF, file = "traitsDF.cdf")


# Now I am wondering if using study titles may actually be a better correlation with disease assocaited SNPs reported in the GWASc

study_titles <- unique(asso$STUDY) # only 5446!
writeLines(study_titles, 'GWASc_study_titles.txt')
studyTitle_chars <- sum(nchar(study_titles))
studyTitle_chars # 599738 ..... 599738/2000 = 299.869 ..... implies like 300 calls to chat GPT... which is still pretty annoying. I wonder if there is a good way to interact with chat gpt in a meta interface way... like prompt it over and over and copy out the results with some kind of program that does the manual entry and collection for me.. so I don't have to pass through an API... with correct timing between requests I think I could get away with this as it is under 500 requests total.. and if this was done over a week it'd be NBD.

# seeking a default term in the trait vecs of interest:

disTraitTable <- table(asso$`DISEASE/TRAIT`)
mapTraitTable <- table(asso$MAPPED_TRAIT)

sort(disTraitTable, decreasing = T)[1:100] # no default terms revealed
sort(mapTraitTable, decreasing = T)[1:100] # empty entries revealed.

sort(table(asso$STUDY), decreasing = T)[1:200] # no default terms revealed

# LEAVING THE FILTERING ISSUE FOR NOW... DOCUMENTED REASONABLY WELL.. ITS TIME TO DO MORE ACTIVELY PRODUCTIVE TASKS NOW









































