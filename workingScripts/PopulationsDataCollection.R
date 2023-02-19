# Data collection

#Setup Environment
source("bootCalls2.R")
load("./WorkingData/GwasAssocitions.rda")
library(GWASpops.pheno2geno)

# getPopsData working version... Second version below is for code edits / func updates.
getPopsData <- function(rsIDChunkList, nChunks, startingChunk, reportNumErrors = TRUE){
  # end chunk = nChunks+startingChunk and will be the next starting chunk and will not be called for in a given run. The final ending chunk must thus be 2499
  retList <- list()

  for (i in 1:nChunks){
    chunk <- (startingChunk + i - 1)
    chnkName <- paste0("Chunk_", chunk, "_CONT")

    retList[[chnkName]] <- tryCatch(
      expr = {
        GWASpops.pheno2geno:::get_ensVariants(rsIDChunkList[[chunk]], population_data = TRUE)
      },

      error = function(e){
        warning(paste0("Error occured for ", chnkName))
        return(chnkName) # this should just return a character vector instead of a list which will be our means of identifying error counts
      }
    )
  }

  fileName <- paste0("chunk", startingChunk, "-", (startingChunk+nChunks-1) , ".rds")

  setwd("D:\\Programming\\R_projects\\Kulathinal_Lab\\GWASpops.pheno2geno\\workingData\\unprocessedChunks")
  save(retList, file = fileName)
  setwd("../")

  if(reportNumErrors){
    numErrors <- sum(sapply(retList, is.character))
    message(paste0("Number of empty chunks returned: ", numErrors))
    message("\nEmpty chunks are returned as character vectors when an error is caught, or when the API fails to return expected data.\n")
  }

  return(retList)
}

# Load GWAS data, save as data object
asso <- data.table::fread("D:\\Kulathinal_files\\GWASc_All_Data_20230114\\gwas_catalog_v1.0.2-associations_e108_r2023-01-14.tsv")

getwd()
setwd('./WorkingData')
save(asso, file ="rsIDlist_GWAS_asso.rda") # rda can take multiple objects... but also single objects.
setwd("../")




# Create chunks of 100 rsIDs each in list

allrsID <- asso$SNPS # 470406 elements
allrsID <- unique(asso$SNPS) # 249792 elements
head(allrsID)

numChunks <- length(allrsID)/100 # 2497.92 == 2498 chunks

class(allrsID) # character vec

# generates list of 2498 items
allrsID <- split(allrsID, ceiling(seq_along(allrsID)/100))

# Naming entries of 'allrsID'

#  first generating vector with names
nameVec <- rep("Chunk_", 2498)
for (i in 1:length(nameVec)){
  nameVec[i] <- paste0(nameVec[i], i)
}

# Then apply those names
names(allrsID) <- nameVec

#saving new data object
save(asso, allrsID, file ="GwasAssocitions.rda")

##############################################################################

# UPDATING allrsID to not include ids specified by chromosome position

allrsID_u <- unique(asso$SNPS)
allchr <- grep("^chr", allrsID_u) #vector with 13192 entries
13192/249792 # = 0.05281194 = 5.3% of the unique SNPs are being cut out of the data set
249792 - 13192 # = 236600 .... should be the number of remaining elements after removal..

allrsID_u <- allrsID_u[!allchr] # this should work if they supported it.. but they don't in R... so we have to take a roundabout solution instead. i.e, cannot subset with inverse indices.. only positive indices.

# Create index list of good indeces, then subset our rsID list with that index vector
sdTest <- setdiff(1:length(allrsID_u), allchr) # generate list of indices with indices from 'allchr' removed
allrsID_u <- allrsID_u[sdTest]



# Generate new Chunk list, naming and saving new object for future use:
#
allrsID <- split(allrsID_u, ceiling(seq_along(allrsID_u)/100)) # 2366 elements now

nameVec <- rep("Chunk_", 2366)
for (i in 1:length(nameVec)){
  nameVec[i] <- paste0(nameVec[i], i)
}
names(allrsID) <- nameVec
save(asso, allrsID, file ="GwasAssocitions.rda") # successfully wrote new data object



# UPDATE 2: removing multientries and setting aside for later.

allrsID_u2 <- unique(asso$SNPS)
allchr <- grep("^chr", allrsID_u2)
sdMask <- setdiff(1:length(allrsID_u2), allchr)
allrsID_u2 <- allrsID_u2[sdMask]

allMultiEntries <- grep(";", allrsID_u2)
sdMask2 <- setdiff(1:length(allrsID_u2), allMultiEntries) # 1258 multientries
allrsID_u2 <- allrsID_u2[sdMask2]
236600 - 1258 # = 235342   .... Should be the remaining number of values within allrsID_u2... success.. numbers are correct


allrsID <- split(allrsID_u2, ceiling(seq_along(allrsID_u2)/100)) #2354 elements now

nameVec <- rep("Chunk_", 2354)
for (i in 1:length(nameVec)){
  nameVec[i] <- paste0(nameVec[i], i)
}
names(allrsID) <- nameVec


# saving updated objects again
save(asso, allrsID, file ="GwasAssocitions.rda")

# saving incices of multiEntry rows
setwd("./workingData")
save(allMultiEntries, allrsID_u2, file ="multiEntrySNPs.rda")
setwd("../")


### Checking for non-expected form rsID entries in allrsID

allrsFlat <- as.character(flatten(allrsID))
goodEntries <- grep("^rs\\d+[^\\s]", allrsFlat) # of 235,342 entries 234,303 were hit...
badEntries <- setdiff(1:length(allrsFlat), goodEntries) # 1039 entries here
rsBad <- allrsFlat[badEntries]
rsBad # we see a wide variety of formats coming out here... good thing I did this kind of positive filtering. Only things that pass a test stay that is. (not sure if this is the appropriate term)

# UPDATE 3: updating the rsID list once again.. should have done this from the start.. Learning

allrsFlat <- allrsFlat[goodEntries]

allrsID <- split(allrsFlat, ceiling(seq_along(allrsFlat)/100)) #2344 elements
nameVec <- rep("Chunk_", 2344)
for (i in 1:length(nameVec)){
  nameVec[i] <- paste0(nameVec[i], i)
}
names(allrsID) <- nameVec

setwd("./workingData")
save(asso, allrsID, file ="GwasAssocitions.rda") # now our rsID list has 2344 elements which have been positively filtered, meaning the format for each entry has been strictly verified using regex.
setwd("../")
getwd()


#### SETTING UP MORE GRANULAR allrsID list as allrsID_ch10

allrsID_ch10 <- split(allrsFlat, ceiling(seq_along(allrsFlat)/10)) # 23431 elements

nameVec <- rep("Chunk_", 23431)
for (i in 1:length(nameVec)){
  nameVec[i] <- paste0(nameVec[i], i)
}
names(allrsID_ch10) <- nameVec

setwd("./workingData")
save(asso, allrsID, allrsID_ch10, file ="rsIDlist_GWAS_asso.rda") # success
setwd("../")
getwd()

###############################################################################





# Function development ------------------------------------------------------
# create empty list,
# for i in 1:nChunks{ call data collecting func for chunk startingChunk+i; store data in empty list; if error arises, report which chunk it failed on by appending value to the list (so add "chunk_startingCHunk+i produced an error)} then continue.
# fileName <- paste0("chunk", startingChunk, "-", startingChunk+nChunks, ".rds")
#
# setwd("./workingData")
# save(retList, file = fileName)
# setwd("../")


# THINGS I NEED TO TEST FOR FUNC TO WORK: ---------------------------------

# 1. error catching generally. Need to ensure my function behaves as desired.
#   - meaning the func can have an error occur but still work just fine
# 2. error catching when rsID is invalid
# 3. countErrors needs to function properly
#

# getPopsData w/ error catching -------------------------------------------



getPopsData <- function(rsIDChunkList, nChunks, startingChunk, reportNumErrors = TRUE){
  # end chunk = nChunks+startingChunk and will be the next starting chunk and will not be called for in a given run. The final ending chunk must thus be 2499
  retList <- list()

  for (i in 1:nChunks){
    chunk <- (startingChunk + i - 1)
    chnkName <- paste0("Chunk_", chunk, "_CONT")

    retList[[chnkName]] <- tryCatch(
      expr = {
        GWASpops.pheno2geno:::get_ensVariants(rsIDChunkList[[chunk]], population_data = TRUE)
      },

      error = function(e){
        warning(paste0("Error occured for ", chnkName))
        return(chnkName) # this should just return a character vector instead of a list which will be our means of identifying error counts
      }
    )
  }

  fileName <- paste0("chunk", startingChunk, "-", (startingChunk+nChunks-1) , ".rds")

  setwd("D:\\Programming\\R_projects\\Kulathinal_Lab\\GWASpops.pheno2geno\\workingData\\unprocessedChunks")
  save(retList, file = fileName)
  setwd("../")

  if(reportNumErrors){
    numErrors <- sum(sapply(retList, is.character))
    message(paste0("Number of empty chunks returned: ", numErrors))
    message("\nEmpty chunks are returned as character vectors when an error is caught, or when the API fails to return expected data.\n")
  }

  return(retList)
}


####### Wrapper Func #############

grabChunks <- function(data, StartChunk = 1,chunksPerCall = 100, numCalls = 240){
  #numCalls at default of 240 should mean that by default this would just call for all of the data

  allrsID_ch10 <- data

  for(i in 1:numCalls){
    startPoint <- (StartChunk + (i - 1)*chunksPerCall) # incrementally updates starting chunk relative to starting point

    if(startPoint > 23300){
      message("All chunks should be grabbed except the last few")
      return()
    }

    getPopsData(allrsID_ch10, nChunks = chunksPerCall, startPoint, reportNumErrors = FALSE)
  }

  message("call finished without automatic termination; i.e. not all chunks grabbed yet, but specified amount should be saved.")
  return()
}


### grab chunks test call

grabChunks(allrsID_ch10, StartChunk = 501, 1)

# DRAFT SPACE -------------------------------------------------------------

for (i in 1:nChunks){
  chunk <- (startingChunk + i - 1)
  chnkName <- paste0("Chunk_", chunk, "_CONT")

  retList[[chnkName]] <- tryCatch(
    expr = {
      GWASpops.pheno2geno:::get_ensVariants(rsIDChunkList[[chunk]], population_data = TRUE)
    },

    error = function(e){
      warning(paste0("Error occured for ", chnkName))
      return(chnkName) # this should just return a character vector instead of a list which will be our means of identifying error counts
    }
  )
}





# TEST SPACE:  ------------------------------------------------------------

sc = 5
fileName <- paste0("chunk", sc, "-", sc+sc, ".rds") # works

GWASpops.pheno2geno:::get_ensVariants()
?createMT



chunk1 <- allrsID[[1]]
class(chunk1) # "character"

GWASpops.pheno2geno:::get_ensVariants()


debug(getPopsData)
undebug(getPopsData)

# TEST CALLS
testDataRetr <- getPopsData(allrsID, nChunks = 2, startingChunk = 1, reportNumErrors = F)
# FIRST RUN RESULTS: first chunk looks ok but not ideal. First list has 99 elements... 1 dropped, not sure why.. second list is entirely empty... which is confusing as well. I will have to run with the debugger to see where we are losing a whole call before integrating the trycatch() block

testDataRetr_2 <- getPopsData(allrsID, nChunks = 3, startingChunk = 1, reportNumErrors = F)


# WITH NON rsID SNP entries removed from 'allrsID' RETEST BASIC FUNC CALL:
#     still need to determine why we aren't getting the same number of IDs back as we call for at times..
#
testDataRetr <- getPopsData(allrsID, nChunks = 2, startingChunk = 1, reportNumErrors = F)
# post updating the 'allrsID' obj call: I still have 0 entries in the second chunk, and that is strange, but I now have 100 in the first chunk which is promising. The third chunk returned 95 results initially, and there was at least 1 element removed from the call set of chunk2 as a single 'chrxx:..." entry was in there.. suggesting it is possible valid data is being discarded due to some issue within the call set of chunk2. When first debugging to manually inspect I was seeing that the setup for calling looked pretty good / as expected so I still have a mystery as to why this data is disappearing/returning nothing. ALSO we know at least 1 element from chunk2 is now in chunk1, and it is returning data successfully, again suggesting some issue within 2 specifically as chunk 3 successfully returned data..
#  GOING TO call for 5 chunks to see if any patterns emerge wrt the missing data in chunk2

testDataRetr <- getPopsData(allrsID, nChunks = 5, startingChunk = 1, reportNumErrors = F)
# RESULT: chunk 2 and 5 are both empty, the rest (1 3 4) are full of 100 entries each.. obviously this is enough lost data that we must attempt to figure out why data is being discarded like this and how to best resolve the issue...
#  IN THE CASE I CANNOT PINPOINT ANY SPECIFIC ISSUE, break problematic chunks into subchuncks of size 10 to possibly identify troublesome entries... maybe that will help me determine why we see this issue and how to resolve it.

# TESTING AFTER POSITIVE FILTERING OF allrsID object:
testDataRetr <- getPopsData(allrsID, nChunks = 5, startingChunk = 1, reportNumErrors = F)
# RESULT: still missing data on chunk 2 and 5 ... Not clear on why that may be...

allrsID[[1]]
allrsID[[2]]
allrsID[[3]]
allrsID[[5]] # nothing looks inherently weird about the data in these chunks. I think it may be best to do sub-chunk calls of groups of 10 then I can see what is going on a little more granularly

chunck_2sub <- allrsID[[2]]
chunk_2sub <- split(chunck_2sub, ceiling(seq_along(chunck_2sub)/10))

chnk2Test <- getPopsData(chunk_2sub, nChunks = 10, startingChunk = 1, reportNumErrors = F) # calls do happen a little faster.. in the case this means we just lose less data, maybe we should split up our data by chunks of 10 like this.. and go more granular where needed? Could be a useful
# RESULT: sub-chunks 4 and 5 failed... others worked fine.. so maybe we can search ensembl for the ids within and see about missing values.. maybe values with no return crash calls?
chunk_2sub[[4]] # [1] "rs11627967" ✔  "rs13059110" ✔ "rs17036328" ✔ "rs7578326" ✔  "rs201450565" ✔"r  s7217780"  ✔ "rs11466328"  ✔
# "rs201386833"✔(as synonym) "rs2963463"✔   "rs201226733" ✔
# OK every ID is in their system in some capacity... every id has population genetics data coming back too.. (at least its not greyed out on their respective pages)... I could call for them individually and see which one give issues and look deeper into it within their system next... any other approaches to identifying the reason for dropped data?

 # TESTING: length 1 calls for failure in chunk2-subchunk4
chunk_2_4sub <- split(chunk_2sub[[4]],ceiling(seq_along(chunk_2sub[[4]])/1 ))
chnk2_sub4Test <- getPopsData(chunk_2_4sub, nChunks = 10, startingChunk = 1, reportNumErrors = F)
# RESULT: all data grabbed successfully.. so this is strange. I can either try and figure out why calls aren't working when grouped... or just take a granular approach when needed... so I think making a more granular allrsID object with groups of 10, thus chunks of much smaller size would make sense.. I doubt I would exceed the limit of calls to the API this way, and can just account for missing data later on by grabbing all missing chunks and further granularizing them to a point where they work. If any data is missing after this point, its been minimized and will be acceptable.
#




# TESTING ERROR PRODUCTION

getErr_ch2 <- chunk_2sub
getErr_ch2[[1]][2] <- "WILL?.th\"i+_-s Get \\Error\\?"
getErr_ch2[[1]][2]
getErr_ch2 # looks good.

errorProductionTest <- getPopsData(getErr_ch2, nChunks = 2, startingChunk = 1, reportNumErrors = F)
# looks like no error is produced from a bad ID... it just filtered the dang thing out without a problem...
#     How can I intentionally produce an error in order to test my error catching code? ... Going to try a worse string id..
#     next idea for error production is to create a bad version of my POST call and to intentional call it on the nth loop to just test for error catching
#     Could also try to call for bigger chunks than it want to take by mangling my function ....

errorProductionTest <- getPopsData(getErr_ch2, nChunks = 2, startingChunk = 1, reportNumErrors = F) # try again after further string mangling.. this worked!!
# RESULT: Bad request... Great! lets see if we can use this as a way to robustify calling for data next session. For now I want to just setup a call to grab as much data as I can while away today.

# Troubleshooting for func ------------------------------------------------

# Finding missing IDs in Chunk 1 and 3

ch1Names <- allrsID[[1]]
ch1Returned <- names(testDataRetr[[1]])
ch1Missing <- ch1Names[ch1Names %in% ch1Returned] # 76 elements
ch1Missing2 <- ch1Names[ch1Returned %in% ch1Names] # 76 elements ??? wtf

mask <- ch1Returned %in% ch1Names # mask of size 99
mask2 <- ch1Names %in% ch1Returned # mask of size 100 ... left side determines the size of the mask for %in%
sum(mask) #76
sum(mask2) #76

?setdiff
tSDiff <- setdiff(ch1Names, ch1Returned)  # size 24.. so yes my other method is working.. meaning the returned set must have 23 elements not in the call set... which in my mind must point to synonyms... sheesh. I don't know how to deal with synonyms for the analysis
#"rs149253773"    "rs116551911"    "rs114002231"    "rs139789464"    "rs144433536"    "rs147151648"    "rs115100928"    "rs115819854"
#[9] "rs115447191"    "rs114738294"    "rs118087489"    "rs146204659"    "rs114385935"    "rs116668069"    "rs116381494"    "rs115287935"
#[17] "rs116442863"    "rs139244745"    "rs138488080"    "rs115729734"    "rs116418332"    "rs115707823"    "rs116260619"    "chr19:10739143"
#
# ^^^ There is also this "chr19:..." value, which shouldn't exist as its not a valid way to call for rsIDs I think. I need to test that against the API too... Holy shit the amount of slippery data within this project is obscene. Nothing is consistent or clean. Fucking biology.
tSDiff
ch1Names
ABsd <- tSDiff
BAsetDiff <- setdiff(ch1Returned, ch1Names) # 23 elements:
BAsetDiff #"rs211468"   "rs34724414" "rs3134977"  "rs1054026"  "rs9267123"  "rs3129112"  "rs3117581"  "rs3132374"  "rs886422"   "rs9257681"  "rs9268123" "rs78533661" "rs3129794"  "rs7766107"  "rs2517617"  "rs3131112"  "rs1150756"  "rs3129817"  "rs2949719"  "rs3095268"  "rs62445633" "rs2395160" "rs1655930"

# OK SYNONYM THOERY IS CONFIRMED:
#  SEARCHED rs149253773... RETURNED rs78533661 from the Ensembl website. This means that preferred synonyms will automatically be returned... I don't know what this actually means for later stages of the analysis... I know it means that matching the GWAS associations table against the returned data may not create neat and complete tables... as I don't have data on rs149253773 for populations.. I only have it for a synonym.. So, hypothetically I need to convert all rsIDs in GWASc to their preferred synonyms before doing a merge.. which means I need to somehow determine preferred synonyms programmatically and do that at some stage of the data wrangling. The synonyms are found within the `synonyms` col of the returned data at least.. so maybe we just need to mind such interactions when making merges later down the line!

# another note: I needd to go through my chuncks of rsIDs.. infact, I need to recreate them, and remove anything that doesn't fit the rsID schema.



# Test: Error catching ----------------------------------------------------


getErr_ch2 <- allrsID_ch10
getErr_ch2[[2]][2] <- "WILL?.th\"i+_-s Get \\Error\\?"

errorProductionTest <- getPopsData(getErr_ch2, nChunks = 2, startingChunk = 1, reportNumErrors = F)# First error catch.. we should get a successful return with the appropriate default value back if error catching works. also should see an error message.

# SUCCESS. Error has been caught. Default element went in. Warning message printed. Now to make our error checker function


errorProductionTest_2 <- getPopsData(getErr_ch2, nChunks = 3, startingChunk = 1) # added in reporting on how many errors / empty returns occur per call for visibility of the process. Want to verify efficacy now.


# Test: Data saving
#  1. test manual calls
#  2. test source of the same calls
#  3. test saving within for loop (using wrapper function)

# CALL THESE MANUALLY, THEN CALL IN 'testDataSaving.R' ... watch for when new data is actually added to dedicated folder
dataSaving_1 <- getPopsData(allrsID_ch10, nChunks = 2, startingChunk = 1)
dataSaving_2 <- getPopsData(allrsID_ch10, nChunks = 2, startingChunk = 3)
dataSaving_3 <- getPopsData(allrsID_ch10, nChunks = 2, startingChunk = 5)
dataSaving_4 <- getPopsData(allrsID_ch10, nChunks = 2, startingChunk = 7)
dataSaving_5 <- getPopsData(allrsID_ch10, nChunks = 2, startingChunk = 9)


gChunksTest_1 <- grabChunks(allrsID_ch10, StartChunk = 1, chunksPerCall = 1, numCalls = 5) # we expect to see chunks 1-5 called for from this, all saved as their own rds files.. do we see them save as the function runs?

gChunksTest_1 <- grabChunks(allrsID_ch10, StartChunk = 1, chunksPerCall = 5, numCalls = 5) # we expect to see chunks 1-25 called for from this saved in 5 .rds files



# Data Collection ---------------------------------------------------------

size10_1000chunks <- getPopsData(allrsID_ch10, nChunks = 1000, startingChunk = 1, reportNumErrors = F) # encountered memory error... Need to determine how to allocate a certain amount of immovable memory to R and how to gaurentee I will have enough memory to retrieve a maximum chunk size.


# Sesh 4 benchmarking and solving memory issues
#
size10Bench_1 <- system.time(getPopsData(allrsID_ch10, nChunks = 100, startingChunk = 1, reportNumErrors = F))
size10Bench_2 <- system.time(getPopsData(allrsID_ch10, nChunks = 100, startingChunk = 101, reportNumErrors = F))
size10Bench_3 <- system.time(getPopsData(allrsID_ch10, nChunks = 100, startingChunk = 201, reportNumErrors = F))
size10Bench_4 <- system.time(getPopsData(allrsID_ch10, nChunks = 100, startingChunk = 301, reportNumErrors = F))
size10Bench_5 <- system.time(getPopsData(allrsID_ch10, nChunks = 100, startingChunk = 401, reportNumErrors = F))

size10Bench_1
size10Bench_2
size10Bench_3
size10Bench_4
size10Bench_5


# Using grab chunks to grab maximal amount of data.

grabChunks(allrsID_ch10, StartChunk = 1201) # this should just grab data until the cows come home. Will duplicate teh 1301-1400 chunk


grabChunks(allrsID_ch10, StartChunk = 1901)

grabChunks(allrsID_ch10, StartChunk = 2301) # windows decided it wanted to update without permission and killed my process.

grabChunks(allrsID_ch10, StartChunk = 3701)

grabChunks(allrsID_ch10, StartChunk = 5901)

grabChunks(allrsID_ch10, StartChunk = 11901)

grabChunks(allrsID_ch10, StartChunk = 12901)

grabChunks(allrsID_ch10, StartChunk = 16101)


36# Session End Notes: ------------------------------------------------------

#Sesh1:
#TODO tomorrow: use get_ensVariants() to do test runs.. workout the error catching code and use bad iDs to test for errors or that scenario. Do a small test run of like 10 chunks and maybe more if time to get a sense for how long each chunk takes on average.. obviously work out the error counting function... and by night get this script calling for data over night / when you're away from home until all data is collected
# AFTER WE GET ALL DATA: working out some modifications to the transformation funcs will be the next step, and that honestly should mostly be a matter of doing some strong asserting of what is supposed to be transformed, i.e. discarding all things which don't suit our schema early on,
# AFTER EDITING TRANSFORM FUNC: setup architecture to unpack-transform-merge all data into managable chunks depending on how big it is... You will just have to find a reasonable way to organize that data such that it can be used for later analysis.
# WHILE WAITING FOR DATA TO BE GRABBED OR MERGED: work on subsequent tasks such as finding the right libraries for fst / cluster analysis.. do practice runs of those packages and setup the code architecture needed to get it poppin. Also work on a report for the background you built.. be brief but really aim this one at the class grade needs.. just do some bs reporting to show the work you have been putting in.. even if the work is probably 50% worrying about being confused.
#
#
#

# Sesh2:
#
# WHERE TO PICKUP NEXT TIME:
#  I have not integrated trycatch() yet.
#  I have not determined why entire Chunks are returning nothing *****
#  ****** I have filtered the bad rows out of our rsID set completely now. (I hope completely.. maybe I should verify by finding anything that doesn't match a regex search pattern of "^rs\d+[^\s]" .. so starts with rs.. then 1 or more digits.. then NOT whitespace... may need to refine pattern) ✔
#       ****** I have not done an exhaustive verification that all rows of rsIDs are in a the standard form of a single rsID per row.... see regex above for means of achieving verification ✔
#
#  I have not benchmarked how long calls take yet and will not until I know I am getting the majority of data needed
#  I have seen that the function works as inteded thus far.. saving objects with appropriate names in the designated folder
#

# Sesh3:
#
# Still need to test and integrate trycatch() .. found error producing call which can be injected into a small run in the middle for testing
# Still need to benchmark larger calls to determine suitable size of calls for periods of time.. ~how many size10 chunks / hour
#
# Successfully filtered out all incorrect format rsIDs, determined they are just skipped over by the API in the case of not working.
# Determined that the missing chunks problem can be at least managed by making the calls smaller, so I intend to call for chunks of size 10 now and then for empty chunks I will go back through with chunks of size 1 to fill in the gaps, anything missed at that degree of granularity can be considered unimportant
# running a test run without error catching for 1000 chuncks of size 10, so 10,000 rsIDs. Not benchmarking it right now though... Would be worth benchmarking my next big run to understand roughly how big a request I should aim to make when I will be away for x amount of time.
#
# So in the current plan I would grab all data... look for all chunk indices which dropped data, further granularize them, then reintegrate that all together for a total data as lists set... then would need to transform that data, then integrate with GWAS catalog while taking the synonyms into account to ensure all GWAS SNPs with matches are found... Will also need to narrow down the GWAS table by good SNPs as well since there is nothing to integrate for any rows/samples which lack quality SNP format.

# looks like calling for 1k groups of 10 may have asked for more memory than was available.. not sure how to increase available memory... but this obviously poses a serious issue to getting the data I want... after talking to Dr. Pond I also learned about the Dr. Hey ... who is a populations geneticist and also about how there is not one single way to calculate Fst.. thus suggesting that I must be somewhat considerate of how I actually weight the factors which differentiate populations in the context of the project.. and thus learning more on this and reporting on it would be a good exercise as a part of the project long term... ...... Guess a next step is to figure out the memory issue. Understand how to allow more write memory for the curl fetching / R generally... not sure what all this means, but opening up maximal memory, and setting up my function to consume minimal memory where it isn't actually necessary is the way forward here. The call was obviously stalled out. Not sure how far it got.. but I would guess about half way or so.

#
# Error in curl::curl_fetch_memory(url, handle = handle) :
#   Operation was aborted by an application callback
# 7.
# curl::curl_fetch_memory(url, handle = handle)
# 6.
# request_fetch.write_memory(req$output, req$url, handle)
# 5.
# request_fetch(req$output, req$url, handle)
# 4.
# request_perform(req, hu$handle$handle)
# 3.
# POST(baseURL, content_type("application/json"), accept("application/json"),
#      body = rsID_Array) at get_ensVariants.R#65
# 2.
# GWASpops.pheno2geno:::get_ensVariants(rsIDChunkList[[chunk]],
#                                       population_data = TRUE)
# 1.
# getPopsData(allrsID_ch10, nChunks = 1000, startingChunk = 1,
#             reportNumErrors = F)
#

#  ^^^ the above error seems to occur consistently when I cancel calls to the API. The reason why I cancel though, is because sometimes the calls just freeze, where the function in still running, but for some reason things are frozen. Not sure its a memory issue at this point. Not sure why the freeze is happening. Not honestly clear on how to gain insight into this sort of bug either... would need to read into API trouble shooting.. if I could somehow view the stage of data requesting my machine is freezed at, whether its saving something in RAM temporarily and lacking space, or waiting for a lost response that is never coming in with no default means of managing such lost responses... I am not sure what all could go wrong, and thus would need to think more, read more and maybe consult Ben if I wanted to get to the bottom of this.

# SESH 4:









































































































