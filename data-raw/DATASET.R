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

GWAS.asso.study.data <- createMT("./data_and_testData/GWASc_air_pollution_search", varAnnotations = FALSE)
titleVec <- as.character(GWAS.asso.study.data[ ,'Title'])
Encoding(titleVec) <- 'UTF-8'
GWAS.asso.study.data[ , 'Title'] <- titleVec # Nice fast way to fix encoding of a vector in a structure. Just need to verify the conversion preserves the human-read meaning of the contents of the vec.


usethis::use_data(GWAS.asso.study.data, overwrite = TRUE)


# masterList creation for dev use --------------------------------------------

masterList <- createMT('./exampleData/air_pollution', population_data = T)
use_data(masterList, internal = T)


