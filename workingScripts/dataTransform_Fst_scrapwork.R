# scrapwork which distracts from core actions of script is held here:


# setup env ---------------------------------------------------------------

getwd()
setwd("D:/Programming/R_projects/Kulathinal_Lab/GWASpops.pheno2geno/")

library(GWASpops.pheno2geno)
source("bootCalls2.R")
source("./R/fst_funcs.R")
source("./workingScripts/dataTransform_Fst.R")
load("./WorkingData/GwasAssocitions.rda")
thousGenPops <- Populations[grep("1000GENOMES", Populations$Population_Abbreviation)]


# -------------------------------------------------------------------------




# read data from memory,
# transform data using modified package transform func (only modify adhoc)
# calculate fst stats for that chunk
# save as object to memory

# transform_fst_save() takes in the large GWAS association table, a number of chunks, and starting chunk and transforms the lists within the chunks into a list of DFs and lists OF DFs.. 4 objects specified above. Fst will be calculated within then the object will be saved.
transform_fst_save <- function(GWAS_associations,
                               numChunks,
                               startChunk = 1,
                               Fst_populations,
                               return_DS = FALSE,
                               saveData = TRUE){

  # load data from memory and flatten for processing
  variantList <- load_n_flatten(numChunks = numChunks, startChunk = startChunk)

  # compose list then tranform into GWASpops.geno2pheno masterList format
  dataList <- list(GWAS_associations, variantList)
  masterList <- ensListTransform_mod(dataList, TRUE)

  # calculate Fst, delete redudant data vals, and discard multiallelic sites
  fstList <- hudsonFst_alleleList(masterList[[2]], Fst_populations, TRUE,TRUE)

  # make single table of Fst Value list, then bind to masterList data structure
  #
  fstList <- fill_rows(fstList) # making all sublists compatible for binding together as data.frame
  names <- names(fstList)
  fstDF <- cbind.data.frame(fstList)
  colnames(fstDF) <- names
  fstDF <- as.data.frmae(t(fstDF)) # transpose s.t. rows are alleles, cols are population-pairs

  masterList[['Fst_per_allele']] <- fstDF

  # save new data structure in memory
  if(saveData){
    fileName <- paste0("D:\\Programming\\R_projects\\Kulathinal_Lab\\GWASpops.pheno2geno\\workingData\\fst_GWAS_annotation_lists",'fullData_', numChunks, '_', startChunk)
    save(masterList, file = fileName)
  }

  # return nothing if desired, or resulting data structure if desired.
  if(return_DS){
    return(masterList)
  } else{
      return()
  }

}



# PRIMARY FUNCS PREPPING FOR transform_fst_save ---------------------------------------------



# planning inspecting .... ------------------------------------------------

# flatten chunk list into single list of rsIDs (get rid of empty entries )
# filter 'asso' (MT) according to names of list elements (variants within a single flat list)
# combine filtered 'asso' (fAsso) with flattened list of rsIDs .... [fAsso element 1 ; variant list element 2]
# pass into transform func.
#
# .... probably create modifed transform func which trucates data forcibly in order to ensure successful transforms.
#
#  After all above done, compose transform_fst_save which perfroms steps above and will reasonably create the 4 specified data structures needed for analyses.


# test flattening

testLNF <- load_n_flatten(numChunks = 3) # looks good. runs fast.

# now to filter asso and transform

fAsso <- asso[asso$SNPS %in% names(testLNF)] # down to 3523 obs

# tranform
dList <- list(fAsso, testLNF)
tTransform <- ensListTransform_mod(dList, T) # notably no failure of transformation

fassoCOPY <- fAsso
data.table::setnames(fassoCOPY, old = 'SNPS', new = 'VariantID') # cannot use normal syntax to reference names on a data.table. 'DF['colname']'
names(fassoCOPY)


##############################################################################################################################################
####### After transform we calculate Fst per SNP, we bind all snps into one DF, transpose and append to the list before saving ###############
##############################################################################################################################################
setwd("..")
getwd()
source("./R/fst_funcs.R")

fstTest <- hudsonFst_alleleList(tTransform[[2]], Populations, TRUE) # taking a long time ... might be having some issue with populations including those without listed sample_count
#      many warnings stating that "NAs introduced by coercion"

thousGenPops <- Populations[grep("1000GENOMES", Populations$Population_Abbreviation)]
fstTest <- hudsonFst_alleleList(tTransform[[2]], thousGenPops, TRUE) # started at 3:20 finish at 3:26 .. so 6 mins for 3 100 size chunks rn. w/ 234 chunks that 468 mins for just these transform w/o anything else.. not terrible 7.8 hrs for total time for calculating Hudson Fst... though we want Wrights.. which should be the same time probably?

fstTest <- hudsonFst_alleleList(tTransform[[2]][1:100], thousGenPops, TRUE)
fstTest <- do.call(cbind, fstTest) # rbind isn't the call we want.. we want the rows to be matched up, issue is there are many missing values.. so would be necessary to fill missing values with NULL or NA .. maybe dplyr has this built in?


tTransform[['fstTest']] <- fstTest

# fill rows func test: ----------------------------------------------------

filledTest <- fill_rows(fstTest) # error
debug(fill_rows)
undebug(fill_rows)


# THIRD SESH START --------------------------------------------------------
# 2-27-23 ... FIXED ERROR. wasn't what I thought it was. chat gpt wasn't quite able to get things right, but it was really my poor instructions.
#
# need to fix the issue where extra values were added to the biggest row in our test data set within perAllelFst..


fstTest_2noNA <- purrr::discard(fstTest_2, is.logical) # this worked for discarding them 99 - 76 = 23 total multiallelics discarded.

# now we move onto binding snps together, transposing then appending.. then writing code up for transform_fst_save() and testing it!

fstTest_3 <- hudsonFst_alleleList(tTransform[[2]][1:100], thousGenPops, TRUE) # should autodiscard multiallelics SUCCESS

filledTest_3 <- fill_rows_gpt(fstTest_3) # looks good.

# filledTest_3bind <- dplyr::bind_cols(filledTest_3) # currently discarding names.. need to rename cols first or store?

ns <- names(filledTest_3)

filledTest_3bind <- cbind.data.frame(filledTest_3) # this is faster than dplyr.

# checking for order sameness.. if order is preserved by cbind.data.frame, we'll have no issue extracting and re implementing names ...

# according to my inspection and chat GPT order is presevered, so simple name assignment works here.

colnames(filledTest_3bind) <- ns

filledTest_3bind_t <- as.data.frame(t(filledTest_3bind)) # looks good

tTransform[["HudsonFst_PerAllele"]] <- filledTest_3bind_t


# Sesh 3 end --------------------------------------------------------------
# ENDING NOTES:
#  successfully debugged 2 prmary issues in 'fill_rows_gpt()' and dealing with multiallelic sites within multiple funcs. Now I need to draft up transfrom_fst_save() and start running it to get data all processed... HOWEVER. I should work out how to calculate Wright's FST before actually getting all my data calculated...



GWASpops.pheno2geno::Populations


# planning inspecting .... ------------------------------------------------

load("./workingData/unprocessedChunks/chunk1-100.rds")

chuk2 <- load("./workingData/unprocessedChunks/chunk101-200.rds") # doesn't work for renamming

?load
# make small script to read in all chunks and rename them.. unless there is some option to load objects with custom name? (no option... )

# test flattening
tl <- list()

load("./workingData/unprocessedChunks/chunk1-100.rds")
tl[["ret1"]] <- retList
load("./workingData/unprocessedChunks/chunk101-200.rds")
tl[["ret2"]] <- retList

tl2 <- purrr::flatten(tl) # first layer of flattening
tl2 <- purrr::flatten(tl2) # second. both succeeded


# DETERMINE HOW TO GET NAMED VEC FROM perAlleleFst_transform()

namedVec <- perAlleleFst_transform(tTransform[[2]]$rs6921580, thousGenPops, F)

nVec <- namedVec[,6, drop = FALSE] # this is it!
##############################################################


# refactored for efficiency version: WILL BE BUGGED UNTIL FIXES ARE IMPLEMENTED

fill_rowsREFACTOR_THIS <- function(DF_list){

  largestRowIndex <- which.max(sapply(DF_list, nrow))

  # find largest row
  max_rows <- nrow(DF_list[[largestRowIndex]])

  # get names of largest row as rowNames
  rowNames <- names(DF_list[[largestRowIndex]])

  # for each DF in DF_list add missing rows with NA filled in using rowNames
  for (i in seq_along(DF_list)) {
    if (nrow(DF_list[[i]]) < max_rows) {

      # create DF with NA values, then name the rows, then bind with the DF before reassigning to the DF_list
      missing_rows <- data.frame(matrix(NA, nrow = max_rows - nrow(DF_list[[i]]), ncol = ncol(DF_list[[i]])))
      row.names(missing_rows) <- setdiff(rowNames, row.names(DF_list[[i]]))
      DF_list[[i]] <- rbind(DF_list[[i]], missing_rows)
    }
  }
  return(DF_list)
}


# THIRD SESH START --------------------------------------------------------


# 2-27-23 ... FIXED ERROR. wasn't what I thought it was. chat gpt wasn't quite able to get things right, but it was really my poor instructions.
#
# need to fix the issue where extra values were added to the biggest row in our test data set within perAllelFst..


badRow <- fstTest[which.max(sapply(fstTest, nrow))]
names(badRow) # rs9378815
search <- tTransform[[2]][['rs9378815']] # grabbing the problem data...

unique(search$allele) # "A" "C" "T" "G"
# ^^ here is our problem I think. We are seeing a multi-allelic site which means that when filtering down to calculate hudsonFst there is a problem.. So we need to determine how to deal with these sort of multi-allelic sites ....
#
# There are multiple ways of dealing with this.. for calculating fst where population differentiation is the insight.. it may make sense to just treat each non-ancestral allele as its own unique mutation/variant. In the case where the mutation is meaningful (not a synonymous mutation) we would want this sort of differentiation I think.

attributes(search)
search$allele

# explore data visually
searchSplit <- split(search, search$allele)
ssA <- searchSplit[[1]]
ssC <- searchSplit[[2]]
ssG <- searchSplit[[3]]
ssT <- searchSplit[[4]]

plot(ssA$population, ssA$frequency) # not working ... think this might not work fo non-numeric args
?plot

library(ggplot2)

ggplot(ssA, aes(x = population, y = frequency))+geom_bar(stat = 'identity', fill = 'blue') + ggtitle("A alleles")
ggplot(ssC, aes(x = population, y = frequency))+geom_bar(stat = 'identity', fill = 'green') + ggtitle("C alleles")
ggplot(ssG, aes(x = population, y = frequency))+geom_bar(stat = 'identity', fill = 'orange') + ggtitle("G alleles")
ggplot(ssT, aes(x = population, y = frequency))+geom_bar(stat = 'identity', fill = 'red') + ggtitle("T alleles")

# decided to discard multiallelic sites, as calculating fst for them is fundamentally different.. and thus I should not impose an improper means of calculating without looking into the subject more deeply.
# NOW: need to test creating of fst list by

fstTest_2 <- hudsonFst_alleleList(tTransform[[2]][1:100], thousGenPops, TRUE) # many sites were discarded.. lets see how many

numNA = 0
for(i in 1:length(fstTest_2)){ if(is.na(fstTest_2[[i]])){numNA = numNA+1} } # not going
fstTest_2[[1]]

numNA <- sum(is.na(unlist(fstTest_2))) # not working as it counts NA from within the objects...
fstTest_2[[2]]

fstTest_2noNA <- purrr::compact(fstTest_2) # not working for this case

has_single_na <- function(x) {
  sum(is.na(x)) == 1
}

fstTest_2noNA <- purrr::compact(fstTest_2, has_single_na) # still not working... chat gpt aint doing so hot

fstTest_2noNA <- purrr::discard(fstTest_2, is.logical) # this worked for discarding them 99 - 76 = 23 total multiallelics discarded.

# now we move onto binding snps together, trasposing then appending.. then writing code up for transform_fst_save() and testing it!

fstTest_3 <- hudsonFst_alleleList(tTransform[[2]][1:100], thousGenPops, TRUE) # should autodiscard multiallelics SUCCESS

filledTest_3 <- fill_rows_gpt(fstTest_3) # looks good.

# filledTest_3bind <- dplyr::bind_cols(filledTest_3) # currently discarding names.. need to rename cols first or store?

ns <- names(filledTest_3)

filledTest_3bind <- cbind.data.frame(filledTest_3) # this is faster than dplyr.

# checking for order sameness.. if order is preserved by cbind.data.frame, we'll have no issue extracting and re implementing names ...

# according to my inspection and chat GPT order is presevered, so simple name assignment works here.

colnames(filledTest_3bind) <- ns

filledTest_3bind_t <- as.data.frame(t(filledTest_3bind)) # looks good

tTransform[["HudsonFst_PerAllele"]] <- filledTest_3bind_t


# Sesh 3 end --------------------------------------------------------------
# ENDING NOTES:
#  successfully debugged 2 prmary issues in 'fill_rows_gpt()' and dealing with multiallelic sites within multiple funcs. Now I need to draft up transfrom_fst_save() and start running it to get data all processed... HOWEVER. I should work out how to calculate Wright's FST before actually getting all my data calculated...

# I should probably clean this workspace a bit at the start of next session to orient myself for clean efficient work.




# DEPRECATED.... SAVING CODE BELOW WORKSPACE NOW::

# read, rename script dev -------------------------------------------------


renameData <- function(numChunk, startChunk = 1) {

  setwd("D:\\Programming\\R_projects\\Kulathinal_Lab\\GWASpops.pheno2geno\\workingData\\unprocessedChunks")

  for(i in 1:numChunk){
    s <- startChunk + (i-1)*100
    end <- s + 99
    dataName <- paste0("chunk", s, "-", end, ".rds")
    load(dataName)
    newObj <- paste0("chunk", s, "-", end)
    assign(newObj , retList)
    save(list = newObj, file = dataName) # assign() should properly rename objects which all will have the name 'retList'
  }

  setwd("D:/Programming/R_projects/Kulathinal_Lab/GWASpops.pheno2geno/")

}


# -------------------------------------------------------------------------

names(tl2[1])

names(tl2)[1:20]
# testing save to see what names gets stored when saving objects of a list

save(tl2[1], file = "testSave.rds")
# Error in save(tl2[1], file = "testSave.rds") : object ‘tl2[1]’ not found

# so no we cannot save elements from a list like this ... lets check save options to see if we can modify the object name
?save
save(list = names(tl2), file = "testSave.rds") # getting error, because its trying to take this character vec as a list of object names... I am trying to specify the name of the object saved.

?assign # this is what I want.

assign(names(tl2[1]), tl2[1]) # it fuckin worked. That is the good shit.


# testing rename func DEPRECATED!!!!!!!!! -----------------------------------------------------
debug(renameData)
renameData(1) # works now

load("./workingData/unprocessedChunks/chunk1-100.rds")

renameData(233)
#ERROR Error in readChar(con, 5L, useBytes = TRUE) : cannot open the connection ... some seem to work... having problems with others..

getwd()
load("chunk2001-2100.rds") # still ret list
load("chunk1501-1600.rds")
load("chunk1401-1500.rds")
load("chunk1401-1500.rds")


# missing chunk 1301-1400 .. this is our issue.

grabChunks(allrsID_ch10, 1301, numCalls = 1)

renameData(220, startChunk = 1301) # restarting after grabbing missing chunk

# Don't think I need to rename data.. think it was unnecessary from the start... mistakes are ok.
#
# Renaming wrong named chunks
setwd('./unprocessedChunks/')
load("chunk1-100.rds")
load("chunk101-200.rds")
load("chunk201-300.rds")
load("chunk301-400.rds")
load("chunk401-500.rds")
load("chunk501-600.rds")
load("chunk601-700.rds")
load("chunk701-800.rds")
load("chunk801-900.rds")
load("chunk901-1000.rds")
load("chunk1001-1100.rds")
load("chunk1101-1200.rds")
load("chunk1201-1300.rds")
load("chunk1301-1400.rds")


retList <- `chunk1-100`
save(retList, file = "chunk1-100.rds")
retList <- `chunk101-200`
save(retList, file = "chunk101-200.rds")
retList <- `chunk201-300`
save(retList, file = "chunk201-300.rds")
retList <- `chunk301-400`
save(retList, file = "chunk301-400.rds")
retList <- `chunk401-500`
save(retList, file = "chunk401-500.rds")
retList <- `chunk501-600`
save(retList, file = "chunk501-600.rds")
retList <- `chunk601-700`
save(retList, file = "chunk601-700.rds")
retList <- `chunk701-800`
save(retList, file = "chunk701-800.rds")
retList <- `chunk801-900`
save(retList, file = "chunk801-900.rds")
retList <- `chunk901-1000`
save(retList, file = "chunk901-1000.rds")
retList <- `chunk1001-1100`
save(retList, file = "chunk1001-1100.rds")
retList <- `chunk1101-1200`
save(retList, file = "chunk1101-1200.rds")
retList <- `chunk1201-1300`
save(retList, file = "chunk1201-1300.rds")


rm(`chunk1-100`,
   `chunk101-200`,
   `chunk201-300`,
   `chunk301-400`,
   `chunk401-500`,
   `chunk501-600`,
   `chunk601-700`,
   `chunk701-800`,
   `chunk801-900`,
   `chunk901-1000`,
   `chunk1001-1100`,
   `chunk1101-1200`,
   `chunk1201-1300`,
   `chunk1301-1400`)

# EVERYTHING NAMED retList again.
#
#

# -------------------------------------------------------------------------


# Testing transform_fst_save() --------------------------------------------

testRun_1 <- transform_fst_save(asso, 1, 1, thousGenPops, TRUE, TRUE) # ran successfully, seeing far too many population pairs forming, woundering if an environment reset may help.

# POST ENV RESET TESTING:

testRun_2 <- transform_fst_save(asso, 1, 1, thousGenPops, TRUE, TRUE)
# seeing same issue... getting 528 pop-pairs.. which is 30 extra.. gotta figure out what is going wrong.
debug(transform_fst_save)
undebug(transform_fst_save)

# found an entry with an ancestral allele that didn't match to alleles within the data frame.. Which breaks our system.. reset based on calc_ancestralAllele() ... if its possible to have an ancestral allele which is different from ANY found within a population (which I presume is in fact possible), then maybe I should change this code somehow.. but this is my fix for now.
testRun_3 <- transform_fst_save(asso, 1, 1, thousGenPops, TRUE, TRUE) # again too many pop-pairs, 524 this time though..
# looks like a deletion as the allele is going from TT to '-' ... coming for inadequate calc_ancestralAllele


# REWRITTEN:
calc_ancestralAllele <- function(population_alleleDF){

  alleleSet <- unique(population_alleleDF$allele)

  highestSum <- 0
  AA <- ''
  for(allele_character in alleleSet){
    tempSum <- sum(population_alleleDF[population_alleleDF$allele == allele_character, ]$frequency)
    if(tempSum > highestSum){
      highestSum <- tempSum
      AA <- allele_character
    }
  }
  return(AA)
}


# SESH END 3-1-23 : Need to rewrite calc_ancestralallele to handle any characters such as 'TT' and '-', with proper assignment. Then we can get back to testing... once we get a good result we'll start running against more data and inspecting results to ensure everything looks good on say 15% of the data.

# SESH 3-2 start:


testRun_4 <- transform_fst_save(asso, 1, 1, thousGenPops, TRUE, TRUE) # want to see if calc_ancestralAllele works now... SUCCESS! we are seeing the right number of rows.. but I would also like to check that the proper ancestral alleles are being assigned.

# data isn't currently being saved despite the fact it SHOULD be getting saved... will have to work on that. Think I figured out the saving issue


testMany_1 <- transform_fst_save(asso,10, 2,thousGenPops, TRUE, TRUE) # this threw an error:
# Error in readChar(con, 5L, useBytes = TRUE) : cannot open the connection
# In addition: Warning message:
#   In readChar(con, 5L, useBytes = TRUE) :
#
#  ^^ Error was likely due to incorrect referencing written into the load_n_flatten() func.

testMany_1 <- transform_fst_save(asso, 10, 1,thousGenPops, TRUE, TRUE) # this call had an issue where the calculation of Fst only seemed to happen for the first chunk loaded. Gotta see why.


testMany_1 <- transform_fst_save(asso, 10, 2,thousGenPops, TRUE, TRUE) # no error on this call.. Gonna see if expected number of fst come out which means fst num should scale close to linearly along with the number of snps processed.

# it runs.. still not getting the num of Fst calcs we should see....
rownames(testMany_1[[4]])


testMany_2 <- transform_fst_save(asso, 3, 20, thousGenPops, TRUE, TRUE)
# 1 more test, then we inspect directly .... we are actually seeing an amount that is pretty reasonable... as many may be discarded for valid reasons... if nothing else the number of columns is correct, which is promising. Still may want to directly inspect the processing to ensure everything is functioning as expected.
#


# what is the range of Hudson Fst?
#

fst_ret <- testMany_2[[4]]

max_hudson <- max(as.matrix(fstPer), na.rm = TRUE) # 0.9950495
max_hudson # NA ... ok that is just not helpful at all. why woudl this be a thing lol?
# ok found the option to ignor NA values sweet.

min_hudson <- min(as.matrix(fstPer), na.rm=T)
min_hudson # -0.01626984 ... nearly 0 anyway... which I suppose means we can just take any value that is not 0 and treat it as 0, as nothing is really that far below, which suggests its not that meaningful. 1 still appears to be the cap. Thus we will just do this, as finding a way to do Wrights Fst seems a pain in the ass.


debug(transform_fst_save) # debugging to see what the heck is going on inside
testMany_3 <- transform_fst_save(asso, 5, 1, thousGenPops, TRUE, TRUE)

# continue debugging next sesh. Step into HusonFst calculation to see where the discarding is happening, ensure it is correct and take note of it.
