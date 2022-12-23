#using this script to store some optimization data that will help me make decisions about what to optimize



# Variant Data ------------------------------------------------------------

# order of data to test:
# 1. malabs               .47 secs
# 2. airPol               .46 secs
# 3. prostateCancer       2.72 secs ... weird merge was needed
# 4. Colorectal           1.64 secs ...
# 5. SubstanceAbus        1.27 secs
# 6. lungCancer           1.78
# 7. breastCarcinoma      2.5
# 8. IBF                  2.2
# 9. alcConsump           2.28
# 10. neuroticism         4.17
# 11. Int                 6.56

# OK nothing is taking long at all. I consider these times acceptable for the current expectations of this package. Moving on.
#
#


# Populations Data Trasformation  -----------------------------------------


# 1. malabs

# -------------------------------------------------------------------------
# popData <- lapply(popData, function(x) bind_rows(x)): 0.19 sec elapsed
# if(popsData) # 1: 0.19 sec elapsed
#   CONT_Table <- GWASpops.pheno2geno:::rsTable(CONT_noNULL): 0.22 sec elapsed
# singlePopTransform Inner: 1.62 sec elapsed
# singlePopTransform Inner: 1.63 sec elapsed
# singlePopTransform Inner: 1.62 sec elapsed
# singlePopTransform Inner: 1.63 sec elapsed
# singlePopTransform Inner: 1.62 sec elapsed
# singlePopTransform Inner: 1.69 sec elapsed
# singlePopTransform Inner: 1.86 sec elapsed
# singlePopTransform Inner: 1.73 sec elapsed
# singlePopTransform Inner: 1.66 sec elapsed
# singlePopTransform Inner: 1.64 sec elapsed
# singlePopTransform Inner: 1.61 sec elapsed
# singlePopTransform Inner: 1.69 sec elapsed
# singlePopTransform Inner: 1.65 sec elapsed
# singlePopTransform Inner: 1.63 sec elapsed
# singlePopTransform Inner: 1.56 sec elapsed
# singlePopTransform Inner: 1.58 sec elapsed
# singlePopTransform Inner: 1.72 sec elapsed
# singlePopTransform Inner: 1.62 sec elapsed
# singlePopTransform Inner: 1.55 sec elapsed
# singlePopTransform Inner: 1.62 sec elapsed
# singlePopTransform Inner: 1.65 sec elapsed
# singlePopTransform Inner: 1.53 sec elapsed
# singlePopTransform Inner: 1.53 sec elapsed
# singlePopTransform Inner: 1.58 sec elapsed
# singlePopTransform Inner: 1.61 sec elapsed
# singlePopTransform Inner: 1.53 sec elapsed
# singlePopTransform Inner: 1.56 sec elapsed
# singlePopTransform Inner: 1.58 sec elapsed
# singlePopTransform Inner: 1.59 sec elapsed
# singlePopTransform Inner: 1.56 sec elapsed
# singlePopTransform Inner: 1.6 sec elapsed
# singlePopTransform Inner: 1.54 sec elapsed
# singlePopTransform Inner: 1.58 sec elapsed
# singlePopTransform Inner: 1.66 sec elapsed
# singlePopTransform Inner: 1.55 sec elapsed
# singlePopTransform Inner: 1.53 sec elapsed
# singlePopTransform Inner: 1.57 sec elapsed
# singlePopTransform Inner: 1.55 sec elapsed
# singlePopTransform Inner: 1.53 sec elapsed
# singlePopTransform Inner: 1.57 sec elapsed
# singlePopTransform Inner: 1.56 sec elapsed
# singlePopTransform Inner: 1.59 sec elapsed
# singlePopTransform Inner: 1.63 sec elapsed
# singlePopTransform Inner: 1.57 sec elapsed
# singlePopTransform Inner: 1.57 sec elapsed
# singlePopTransform Inner: 1.58 sec elapsed
# singlePopTransform Inner: 1.68 sec elapsed
# singlePopTransform Inner: 1.55 sec elapsed
# singlePopTransform Inner: 1.53 sec elapsed
# singlePopTransform Inner: 1.55 sec elapsed
# singlePopTransform Inner: 1.54 sec elapsed
# singlePopTransform Inner: 1.6 sec elapsed
# singlePopTransform Inner: 1.53 sec elapsed
# singlePopTransform Inner: 1.53 sec elapsed
# singlePopTransform Inner: 1.56 sec elapsed
# singlePopTransform Inner: 1.61 sec elapsed
# singlePopTransform Inner: 1.55 sec elapsed
# singlePopTransform Inner: 1.59 sec elapsed
# singlePopTransform Inner: 1.55 sec elapsed
# singlePopTransform Inner: 1.58 sec elapsed
# singlePopTransform Inner: 1.64 sec elapsed
# singlePopTransform Inner: 1.55 sec elapsed
# singlePopTransform Inner: 1.53 sec elapsed
# singlePopTransform Inner: 1.73 sec elapsed
# singlePopTransform Inner: 1.55 sec elapsed
# singlePopTransform Inner: 1.55 sec elapsed
# singlePopTransform Inner: 1.57 sec elapsed
# singlePopTransform Inner: 1.57 sec elapsed
# singlePopTransform Inner: 1.57 sec elapsed
# singlePopTransform Inner: 1.61 sec elapsed
# singlePopTransform Inner: 1.57 sec elapsed
# singlePopTransform Inner: 1.54 sec elapsed
# singlePopTransform Inner: 1.58 sec elapsed
# lapply(Populations$Population_Abbreviation, function(x) GWASpops.pheno2geno:::singlePopTransform(masterList[[2]], targetPopulation = x)): 116.34 sec elapsed
# START: 116.78 sec elapsed


# AFTER MASKING IS USED INSTEAD OF SELECT() -------------------------------------------

# popData <- lapply(popData, function(x) bind_rows(x)): 0.2 sec elapsed
# if(popsData) # 1: 0.2 sec elapsed
#   CONT_Table <- GWASpops.pheno2geno:::rsTable(CONT_noNULL): 0.21 sec elapsed
# singlePopTransform Inner: 1.01 sec elapsed
# singlePopTransform Inner: 1.02 sec elapsed
# singlePopTransform Inner: 1.04 sec elapsed
# singlePopTransform Inner: 0.91 sec elapsed
# singlePopTransform Inner: 0.91 sec elapsed
# singlePopTransform Inner: 0.89 sec elapsed
# singlePopTransform Inner: 0.89 sec elapsed
# singlePopTransform Inner: 0.89 sec elapsed
# singlePopTransform Inner: 0.89 sec elapsed
# singlePopTransform Inner: 0.89 sec elapsed
# singlePopTransform Inner: 0.89 sec elapsed
# singlePopTransform Inner: 0.87 sec elapsed
# singlePopTransform Inner: 0.91 sec elapsed
# singlePopTransform Inner: 0.89 sec elapsed
# singlePopTransform Inner: 0.91 sec elapsed
# singlePopTransform Inner: 0.9 sec elapsed

#loop is already substantially faster!
#

# SECOND MASK ADDED INSTEAD OF SELECT()

# popData <- lapply(popData, function(x) bind_rows(x)): 0.21 sec elapsed
# if(popsData) # 1: 0.21 sec elapsed
#   CONT_Table <- GWASpops.pheno2geno:::rsTable(CONT_noNULL): 0.24 sec elapsed
# singlePopTransform Inner: 0.02 sec elapsed
# singlePopTransform Inner: 0.02 sec elapsed
# singlePopTransform Inner: 0 sec elapsed
# singlePopTransform Inner: 0.01 sec elapsed
# .......
# singlePopTransform Inner: 0.02 sec elapsed
# singlePopTransform Inner: 0 sec elapsed
# singlePopTransform Inner: 0.01 sec elapsed
# singlePopTransform Inner: 0 sec elapsed
# singlePopTransform Inner: 0.02 sec elapsed
# singlePopTransform Inner: 0.01 sec elapsed
# singlePopTransform Inner: 0 sec elapsed
# singlePopTransform Inner: 0.02 sec elapsed
# singlePopTransform Inner: 0.01 sec elapsed
# lapply(Populations$Population_Abbreviation, function(x) GWASpops.pheno2geno:::singlePopTransform(masterList[[2]], targetPopulation = x)): 0.86 sec elapsed
# START: 1.34 sec elapsed


# -------------------------------------------------------------------------
# 2. airPol
# 3. prostateCancer
# 4. Colorectal
# 5. SubstanceAbus       PRE mask in section 2 of singlePopTransform: START: 5.27 sec elapsed  AFTER NEW MASK START: 4.89 sec elapsed


# 6. lungCancer
# 7. breastCarcinoma
# 8. IBF
# 9. alcConsump
# 10. neuroticism
# 11. Int




FL2 <- lapply(popFreqList, \(x) x[x$population == targetPopulation,])

SPT2 <- singlePopTable2[,names(singlePopTable2) != 'population']




# singlePoptransformOG ------------------------------------------------------
#   PURPOSE OF THIS IS SIMPLY TO STORE THE UNBUGGED BUT INCREBILY SLOW ORIGINAL VERSION OF THIS FUNCTION.
#   Function takes population frequency list of data tables and extracts all
#   rows from each table in the list that match a population label of interest.

#' singlePopTransform
#'
#' single Population Transform, takes in a list of tables which posses population data. The data is structured as 1 variant per table. This function itterates through the set of tables within the list to make tables which are instead, 1 population per table. Data must come from
#'
#' Function takes population frequency list of data tables and extracts all rows from each table in the list that match a population label of interest.
#'
#' @param popFreqList list of population frequency data which is in its unedited state from the ensembl variants endpoint
#' @param targetPopulation name of target population to generate a table for; names of populations can be found in included data 'popFreqMeta'
#' @param includeAncestralAlleles if TRUE ancestral alleles will not be excluded from the new table generated
#'
#' @return
#' A data.frame / tibble data frame
#'
#' @example
#' colombianPopData <- singlePopTransform(popFreqList, targetPopulation = '1000GENOMES:phase_3:CLM')
#'
#' @noRd
# singlePopTransformOG <- function(popFreqList, targetPopulation = 'gnomADg:ALL', includeAncestralAlleles = FALSE){
#   tic('singlePopTransform Inner')
#   filteredList <- lapply(popFreqList, \(x) filter(x,x$population == targetPopulation))
#
#   if(!includeAncestralAlleles){ #removing ancestral allele rows by default.
#     filteredList <- lapply(filteredList, \(x) filter(x, x$allele != attr(x, 'Ancestral_Allele')))
#   }
#
#   singlePopTable <- bind_rows(filteredList, .id = 'VariantID') %>% distinct(.keep_all = FALSE)
#   #    ^^ distinct removes duplicate rows
#   attr(singlePopTable, 'population') <- singlePopTable$population[1]
#   attr(singlePopTable, 'Ancestral_Allele') <- NULL
#   attr(singlePopTable, 'VariantID') <- NULL
#   singlePopTable <- select(singlePopTable, !population)
#   toc()
#   return(singlePopTable)
# }




