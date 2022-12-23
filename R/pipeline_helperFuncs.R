# AUTHOR: Jon Tanner Nelson
# LAST UPDATED: 9-26-2022
# CONTENTS: importGWAS_DataTables(), list2table_associations_studies(), singlePopTransform(), GWAScat_sciConvert(), maxVecLength(), bindTables_keyTerm()
# PURPOSE: These functions all aid in the operation of exposed functions in this package, all functions within this file are to be left unexposed to the user.

#' importGWAS_DataTables
#'
#' Streamlines importing data downloaded from the GWAS Catalog by allowing the user to simply enter a folder containing ONLY GWAS data tables.
#'
#' The folder specified must contain only files which can be converted to data.frames by as.data.frame(). It is recommended to simply nest all
#' downloaded GWAS data in a single folder, then to call this function in order to retrieve all files in a single list.
#'
#' @param folderPath The folder path from which to import GWAS tables from
#'
#' @returns
#' A list of data frames.
#'
#' @examples
#' data <- importGWAS_DataTables("D:\\Parent_Directory\\GWAS_Data")
#'
#' @importFrom dplyr filter
#' @importFrom dplyr bind_rows
#' @importFrom dplyr distinct
#' @importFrom dplyr select
#' @importFrom tidyr separate
#' @importFrom data.table fread
#' @importFrom data.table rbindlist
#'
#' @noRd
importGWAS_DataTables <- function(folderPath){
  #import all files from a single folder using fread from data.table

  # Fill a list with all files from specified folder
  files <- list.files(folderPath, full.names = TRUE)
  lenFiles <- length(files)
  dfList <- vector(mode = 'list', length = lenFiles)

  for (FILE in 1:lenFiles) {
    dfList[[FILE]] <- as.data.frame(fread(files[FILE]))
  }

  # save filenames as names of listed tables imported
  fileNames <- list.files(folderPath)
  fileNames <- sapply(fileNames, function(x) sub('\\..*', '', x)) # removes file extensions from the names
  names(dfList) <- fileNames

  return(dfList)
}


# Merge association and study tables into a master GWASc data tabl --------

# merge together all associations and sudies by rbind(), then merge two different df types together with merge

#purpose is to expidite the process of getting data into R in querieable form

#' list2table_associations_studies
#'
#' Takes in a list of data.frames and merges them all into a single data frame
#'
#' This function is build specifically to work with the studies and associations data tables available on the GWAS Catalog (https://www.ebi.ac.uk/gwas/)
#' Used to be called 'assoStudyMerge_fromImportedList()'
#'
#' @param dfList a list of data frames or data frame like objects: (tibble, data.table)
#'
#' @returns
#' A single data.frame
#'
#' @examples
#' masterTable <- list2table_associations_studies(listOfData)
#'
#' @noRd
list2table_associations_studies <- function(dfList){

  #splitting the DFs into two lists to independently combine each into single DFs
  studyList <- dfList[grep('.*studies.*', names(dfList))]
  assoList <- dfList[grep('.*studies.*', names(dfList), invert = TRUE)] #inversion grabs everything that is not fulfilling the 'studies' grep

  studyDF <- rbindlist(studyList)
  assoDF <- rbindlist(assoList)

  Merge <- merge(assoDF, studyDF)

  finalMerge <- separate(Merge, 'Variant and risk allele',
                         into = c('VariantID', 'GWASc_RiskAllele'), sep = '-')

  finalMerge$GWASc_RiskAllele <- sub("</b>", "", sub("<b>","",finalMerge$GWASc_RiskAllele))

  return(finalMerge)
}


# singlePoptransform ------------------------------------------------------
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
singlePopTransform <- function(popFreqList, targetPopulation = 'gnomADg:ALL', includeAncestralAlleles = FALSE){

  filteredList <- lapply(popFreqList, \(x) x[ x$population == targetPopulation, ]) # using mask instead of select.. and it is many time faster.

  if(!includeAncestralAlleles){ #removing ancestral allele rows by default.
    filteredList <- lapply(filteredList, \(x) x[ x$allele != attr(x, 'Ancestral_Allele'), ])
  }

  singlePopTable <- bind_rows(filteredList, .id = 'VariantID') %>% distinct(.keep_all = FALSE)
  #    ^^ distinct removes duplicate rows
  attr(singlePopTable, 'population') <- singlePopTable$population[1]
  attr(singlePopTable, 'Ancestral_Allele') <- NULL
  attr(singlePopTable, 'VariantID') <- NULL
  singlePopTable <- singlePopTable[ ,names(singlePopTable) != 'population' ]

  return(singlePopTable)
}


#' GWAScat_sciConvert
#'
#' Takes in a string of the form "2 x 10-6" and converts it to an equivalent numeric value.
#'
#' Useful for working with the P-values from GWAS catalog, as they come as strings which are naturally un-graphable
#'
#' @param sciString a scientific notation string of the form '3 x 10-8'. deviation from this specific form will not work with this function.
#'
#' @return
#' numeric
#'
#' @example
#' GWAScat_sciConvert('2 x 10-6')
#'
#' @noRd
GWAScat_sciConvert <- function(sciString){

  num <- as.numeric(sub('x.*', '', sciString))
  deciPower <- as.numeric(sub('.*10-', '', sciString))

  p_Val <- num * 10**-deciPower
}


#' maxVecLength
#'
#' maxVecLength splits a generic vector into multiple vectors contained in a list which all have a maximum length which is specified as an argument
#'
#' N/A
#' @param vector vector to be split up according to max length desired
#' @param maxLen the maximum length of subvectors generated
#'
#' @return
#' List of vectors with same type as 'vector'
#'
#' @example
#' rsIDlist <- maxVecLength(rsIDs, maxLen = 100)
#'
#' @noRd
maxVecLength <- function(vector, maxLen = 1){
  #splits vector into multiple vectors nested in a list which all have
  #a maximum length = maxLen
  list <- split(vector, ceiling(1:length(vector)/maxLen))
}



# Not sure if helper func or personal func --------------------------------
#
# personal funcs are those which I maybe found useful at some point, but are not essential to the packages behavior

#' bindTables_keyTerm
#'
#' Binds all DFs within a list which share a key-term in their name together, producing a single DF.
#'
#' This function is particularly useful when importing several tables of data which share variables and can be easily bound together on the basis of a shared name in their title
#'
#' @param dfList a list of data frames with related names
#' @param term the common term among data frames to be bound together
#'
#' @return
#' Single data.frame is returned
#'
#' @example
#' dataFrame <- bindTables_keyTerm(dataFrameList, term = 'studies')
#'
#' @noRd
bindTables_keyTerm <- function(dfList, term = 'association'){

  term <- paste0('.*', term, '.*')
  ListToBind <- dfList[grep(term, names(dfList))]

  outputDF <- rbindlist(ListToBind)
}




