# package development start up script.
#  This script is here to ease the process of beginning development sessions with a fresh Rstudio environment.

# loading packages for development purposes
# NOT LOADING DATA.TABLE as nonstandard evaluation is a bitch
library(devtools)
library(tidyverse)
library(GWASpops.pheno2geno)

#loading custom data for the package
gwasData <- GWAS.asso.study.data
pops <- Populations
masterList <- testMasterList
