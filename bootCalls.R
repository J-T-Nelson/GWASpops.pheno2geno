# package development start up script.
#  This script is here to ease the process of beginning development sessions with a fresh Rstudio environment.

cat("Functions loaded: \n1. All.ignore()")

# loading packages for development purposes
# NOT LOADING DATA.TABLE as nonstandard evaluation is a bitch
library(devtools)
library(tidyverse)
library(GWASpops.pheno2geno)

#loading custom data for the package
gwasData <- GWAS.asso.study.data
pops <- Populations
masterList <- testMasterList


#IDEAS:

All.ignore <- function(fileName){
  use_git_ignore(fileName)
  use_build_ignore(fileName)
}

GWASdataSets <- c('./GWAS.data/alcohol_consumption', './GWAS.data/breast_carcinoma', './GWAS.data/colorectal_cancer', './GWAS.data/inflammatory_bowel_disease', './GWAS.data/Intelligence', './GWAS.data/lung_cancer', './GWAS.data/malabsorption_syndrome', './GWAS.data/neuroticism_measurement', './GWAS.data/prostate_cancer', './GWAS.data/substance_abuse')

