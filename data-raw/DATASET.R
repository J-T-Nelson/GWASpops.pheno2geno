## code to prepare `DATASET` dataset goes here
##


# populations data --------------------------------------------------------


if(!require(data.table)){
  install.packages('data.table')
  library(data.table, include.only = 'fread')
}

Populations <- fread('./data_and_testData/PopulationFrequencyDataKey.csv', col.names = c('Population_Abbreviation', 'Sample_Count', 'Pop_Ancestry', 'PopAnces_Graph_Labels'))

    #usethis::use_data(DATASET, overwrite = TRUE)
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
usethis::use_data(GWAS.asso.study.data, overwrite = TRUE)


# Code block of saving objects for internal usage by a package  -----------
#   or by a user in the context of a p

internalDataExample <- sample(100)
use_data(internalDataExample, internal = TRUE)
rm(internalDataExample)
setwd('./R')
load('sysdata.rda')

# -------------------------------------------------------------------------



