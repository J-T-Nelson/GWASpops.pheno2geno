## code to prepare `DATASET` dataset goes here
##


# populations data --------------------------------------------------------


if(!require(data.table)){
  install.packages('data.table')
  library(data.table, include.only = 'fread')
}

Populations <- fread('./data_and_testData/PopulationFrequencyDataKey.csv', col.names = c('Population_Abbreviation', 'Sample_Count', 'Pop_Ancestry', 'PopAnces_Graph_Labels'))

usethis::use_data(Populations, overwrite = TRUE)


# GWAS testing data -------------------------------------------------------


if(!require(data.table)){
  install.packages('data.table')
  library(data.table)
}
if(!require(tidyverse)){
  install.packages('tidyverse')
  library(tidyverse)
}
if(!require(data.table)){
  install.packages('data.table')
  library(data.table)
}
setwd('./R')
source('pipeline_helperFuncs.R')
source('pipeline_masterFunc.R')
source('data_retrivalFixingLoading.R')
source('get_ensVariants.R')
source('plotting_graphing.R')
setwd('..')

GWAS.asso.study.data <- createMT("./exampleData/air_pollution", varAnnotations = FALSE)
titleVec <- pull(GWAS.asso.study.data, 'Title')
Encoding(titleVec) <- 'UTF-8'
Encoding(titleVec) # seeing some UTF-8 entries, thus assuming success
GWAS.asso.study.data[ , 'Title'] <- titleVec # Nice fast way to fix encoding of a vector in a structure. Just need to verify the conversion preserves the human-read meaning of the contents of the vec.
Encoding(pull(GWAS.asso.study.data, 'Title')) # seeing some UTF-8 entries, thus assuming success

usethis::use_data(GWAS.asso.study.data, overwrite = TRUE)


# masterList creation for dev use --------------------------------------------

masterList <- createMT('./exampleData/air_pollution', population_data = T)
use_data(masterList, internal = T)

# masterList creation for user --------------------------------------------

masterList <- createMT('./exampleData/air_pollution', population_data = T)
testMasterList <- masterList # masterList is coming from sysdata.rda in data folder.
tempMT <- testMasterList[[1]]

test_titles <- pull(testMasterList[[1]], 'Title')
titles <- pull(tempMT, 'Title')
Encoding(titles) <- 'UTF-8'  # this isn't actually registering all of them as UTF-8 afterwords, but does appear to be getting the elements which were specifically causing an issue.
tempMT[,'Title'] <- titles
Encoding(pull(tempMT,'Title')) #has UTF-8 entries
Encoding(test_titles) #does not have UTF-8 entries
testMasterList[[1]] <- tempMT
Encoding(pull(testMasterList[[1]], 'Title')) #has UTF-8 entries

usethis::use_data(testMasterList, overwrite = T)


newData <- matrix(c(1,1,1,2,2,2,3,3,3), nrow=3)
usethis::use_data(newData, overwrite = T)

# TROUBLE SHOOTING CODE ARCHIVE -------------------------------------------


# having issues extracting character column  ------------------------------
#  I believe the issue was related to non-standard evaluation introduced by data.table package. however, the issue seems pretty mysterious as in the creation of 'titles_3' I had different outcomes under what appeared to be the same project environment conditions. Thus I have determined that until further notice the best way to extract vectors from data.frames and similar tables is to use the dplyr::pull() command. It is consistent.

newTitleCol <- as.character(GWAS.asso.study.data[,'Title'])
tempMT[, 'Title'] <- newTitleCol
testMasterList[[1]][, 'Title'] <- newTitleCol


for(title in titles){
  as.character(title)
}

titles <- testMasterList[[1]][, 'Title']
titles
titles_2 <- pull(tempMT, 'Title') # this works for pulling it out properly..

titles_3 <- tempMT[, 'Title'] # this doesn't work.. BECAUSE OF data.table's nonstandard evaluation... after unloading it this works fine to grab what I was looking for the whole time. .. OKK this worked once.. but now its going back to extracting not vectors, but a column.. which cannot easily be coerced into a character vector .. just use pull() ...
unload("data.table")

data.table::fwrite(encodings_1, file = "testFile1.txt")


# -------------------------------------------------------------------------


