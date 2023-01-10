# package development start up script.
#  This script is here to ease the process of beginning development sessions with a fresh Rstudio environment.

cat("Functions loaded: \n1. All.ignore()\n")

# loading packages for development purposes
# NOT LOADING DATA.TABLE as nonstandard evaluation is a bitch
library(devtools)
library(GWASpops.pheno2geno)
library(httr)
library(jsonlite)

library(tidyverse)

#loading custom data for the package
gwasData <- GWAS.asso.study.data
pops <- Populations
masterList <- testMasterList



# Loading all files from package for immediate usage! ---------------------

# source('./R/data_retrivalFixingLoading.R')
# source('./R/get_ensVariants.R')
# source('./R/pipeline_masterFunc.R')
# source('./R/plotting_graphing.R')
# source('./R/pipeline_helperFuncs.R')

# Loading debugging funcs:

# source('workingScripts/Debugging_APIcall_functions.R') # DEPRECATED

All.ignore <- function(fileName){
  use_git_ignore(fileName)
  use_build_ignore(fileName)
}

GWASdataSets <- c('./GWAS.data/alcohol_consumption', './GWAS.data/breast_carcinoma', './GWAS.data/colorectal_cancer', './GWAS.data/inflammatory_bowel_disease', './GWAS.data/Intelligence', './GWAS.data/lung_cancer', './GWAS.data/malabsorption_syndrome', './GWAS.data/neuroticism_measurement', './GWAS.data/prostate_cancer', './GWAS.data/substance_abuse')


# testing if I can load data in a diff way.  ------------------------------

load('./data/transformed_data_for_graphing/airPollutionPops.rds')

load('./data/transformed_data_for_graphing/alcConsumpPops.rds')

load('./data/transformed_data_for_graphing/bCarcinomaPops.rds')

load('./data/transformed_data_for_graphing/colorectalCancerPops.rds')

load('./data/transformed_data_for_graphing/IBFPops.rds')

load('./data/transformed_data_for_graphing/IntPops.rds')

load('./data/transformed_data_for_graphing/malabsorptionSyndPops.rds')

load('./data/transformed_data_for_graphing/lungCancerPops.rds')

load('./data/transformed_data_for_graphing/neuroticismPops.rds')

load('./data/transformed_data_for_graphing/prostateCancerPops.rds')

load('./data/transformed_data_for_graphing/substanceAbusePops.rds')

load('./data/transformed_data_for_graphing/malabsorptionSyndPops.rds')

#success! We can load data this way.. so an R script can interact with the Global Environment a bit differently than a function run within that environment can
#

# -------------------------------------------------------------------------


