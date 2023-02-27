# data transformation script dev:

# For future data analysis we want data to be organizes into relatively flat structures that are easily queriable and understood
# Building on the scheme set out by my package I will have 3 + 1 primary structures of data saved as objects which can be brought up adhoc when doing analyses.
#
# (if needed any of the objects can be broken up into partitions for the sake of operating them within a my local CP if they become too large... will do so adhoc only)
#
# DS 1: 1 large DF, MasterTable - entire GWAS catalog + ensembl annoations per SNP
# DS 2: large list of DFs, by VariantID - population data contained within
# DS 3: small list of larger DFs, by Population - population data within organized by population instead (this isn't necessary and can be derived from DS 2, so maybe skip it if it seems unmanagable to create)
# DS 4: 1 large DF, Fst List - for all population-pairs all Hudson and/or Write Fst Stats.




# read data from memory,
# transform data using modified package transform func (only modify adhoc)
# calculate fst stats for that chunk
# save as object to memory

# transform_fst_save() takes in the large GWAS assocation table, a number of chunks, and starting chunk and transforms the lists within the chunks into a list of DFs and lists OF DFs.. 4 objects specified above. Fst will be calculated within then the object will be saved.
transform_fst_save <- function(GWAS_associations, numChunks, startChunk = 1){


}


# setup env ---------------------------------------------------------------

getwd()
setwd("D:/Programming/R_projects/Kulathinal_Lab/GWASpops.pheno2geno/")

source("bootCalls2.R")
load("./WorkingData/GwasAssocitions.rda")
library(GWASpops.pheno2geno)


# planning inspecting .... ------------------------------------------------

load("./workingData/unprocessedChunks/chunk1-100.rds")

chuk2 <- load("./workingData/unprocessedChunks/chunk101-200.rds") # doesn't work for renamming

?load
# make small script to read in all chunks and rename them.. unless there is some option to load objects with custom name? (no option... )

# flatten chunk list into single list of rsIDs (get rid of empty entries )
# filter 'asso' (MT) according to names of list elements (variants within a single flat list)
# combine filtered 'asso' (fAsso) with flattened list of rsIDs .... [fAsso element 1 ; variant list element 2]
# pass into transform func.
#
# .... probably create modifed transform func which trucates data forcibly in order to ensure successful transforms.
#
#  After all above done, compose transform_fst_save which perfroms steps above and will reasonably create the 4 specified data structures needed for analyses.
#


# load_n_flatten() --------------------------------------------------------

# realized an error in my processing of the data.. I actually didn't need to rename anything! as I can assign names to the loaded objects adhoc within a list as loading them in... so I will have to manage the first 13 chunks manually somehow... maybe rename manually before processing through here.

load_n_flatten <- function(numChunks, startChunk = 1) {

  setwd("D:\\Programming\\R_projects\\Kulathinal_Lab\\GWASpops.pheno2geno\\workingData\\unprocessedChunks")
  ret_list <- list()

  for(i in 1:numChunks){
    s <- startChunk + (i-1)*100
    end <- s + 99
    dataName <- paste0("chunk", s, "-", end, ".rds")
    load(dataName)
    ret_list[[dataName]] <- retList
  }

  ret_list <- purrr::flatten(ret_list)
  ret_list <- purrr::flatten(ret_list)
  return(ret_list)
}


# -------------------------------------------------------------------------


# test flattening
tl <- list()

load("./workingData/unprocessedChunks/chunk1-100.rds")
tl[["ret1"]] <- retList
load("./workingData/unprocessedChunks/chunk101-200.rds")
tl[["ret2"]] <- retList

tl2 <- purrr::flatten(tl) # first layer of flattening
tl2 <- purrr::flatten(tl2) # second. both succeeded



testLNF <- load_n_flatten(numChunks = 3) # looks good. runs fast.

# now to filter asso and transform

fAsso <- asso[asso$SNPS %in% names(testLNF)] # down to 3523 obs

# tranform
dList <- list(fAsso, testLNF)
tTransform <- ensListTransform_mod(dList, T) # notably no failure of transformation
debug(ensListTransform_mod)
undebug(ensListTransform_mod)

fassoCOPY <- fAsso
data.table::setnames(fassoCOPY, old = 'SNPS', new = 'VariantID') # cannot use normal syntax to reference names on a data.table. 'DF['colname']'
names(fassoCOPY)

# After transform we calculate Fst per SNP, we bind all snps into one DF, transpose and append to the list before saving
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

# -------------------------------------------------------------------------

# DETERMINE HOW TO GET NAMED VEC FROM perAlleleFst_transform()

namedVec <- perAlleleFst_transform(tTransform[[2]]$rs6921580, thousGenPops, F)

nVec <- namedVec[,6, drop = FALSE] # this is it!
##############################################################



# fill rows func test: ----------------------------------------------------

filledTest <- fill_rows_gpt(fstTest) # error
debug(fill_rows_gpt)

# ISSUE: Seeing a row with 515 unique names.. meaning there must be some names where the order is swapped making for 2x possible unique names... need to ensure this cannot happen in code that generates rownames. 496 is the max number that CAN exist for the 32C2 set of name combinations

# START HERE 2-26 ---------------------------------------------------------
#  working towards complete data transformation and fst calculation for all data thus far...
#  finished calling for data, ignoring missing data for now (probably thorughout the rest of the proj.)
#  need to make formula / func for Wrights fst.. reading from wiki is a bit difficult but with a calm mind and some determination I am sure I can work something out
#  need to finish test run of the transformation func 'transform_fst_save()' ... currently hung up on issue detailed just above ... need to ensure that I cannot produce more than the true max number of rows when calculating Fst... seems a bit of a tough one to crack.. as it was only produced once in the test data set I think, and beyond that I am wondering why so many were exactly the right number if this one somehow slipped through the cracks. .. will need a bit of brainstroming for debug strategies.
#  After the first test run of commands works I need to build out transform_fst_save() and run it against the first 24th of the data... see how long it takes.. and see if transformation fails at any point.. do debugging... then hopefully just run the rest of it through and finally see if all the data can be stored in one huge list with the 4 items of data forms desired.



# fillRows func dev ----------------------------------------------------------------

# I need to fill the rows for missing values of population pairs in variants..
#
# making generalized func that can fill in missing rows for DFs with named rows

# chatGPT wrote this version:
fill_rows_gpt <- function(DF_list){

  # find largest row
  max_rows <- max(sapply(DF_list, nrow))

  # get names of largest row as rowNames
  rowNames <- names(which.max(sapply(DF_list, nrow)))

  # for each DF in DF_list add missing rows with NA filled in using rowNames
  for (i in seq_along(DF_list)) {
    if (nrow(DF_list[[i]]) < max_rows) {
      missing_rows <- data.frame(matrix(NA, nrow = max_rows - nrow(DF_list[[i]]), ncol = ncol(DF_list[[i]])))
      row.names(missing_rows) <- setdiff(rowNames, row.names(DF_list[[i]]))
      DF_list[[i]] <- rbind(DF_list[[i]], missing_rows)
    }
  }

  return(DF_list)
}

# refactored for efficiency version:

fill_rows <- function(DF_list){

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



# modifed Transform -------------------------------------------------------


#' ensListTransform_mod
#'
#' @description Transforms list-form data which is produced by `get_ensVariants()` and `createMT(processData = FALSE)` into flat a list of flat tables which can then be used for graphing or viewing data in tabular form.
#'
#' @details
#'
#'
#' @param dataList Data to be transformed. Must be in list format such that a GWAS data table is the first element and the respective data produced from calling Ensembl's REST API Variants endpoint with get_ensVariants() is the second element within the list. Position in the list is critical to successful execution of this function
#' @param popsData populations data transformation option. when TRUE function runs assuming variants from Ensembl REST API have been called with populations option activated, resulting output is different due to this extra population data.
#'
#' @return data.frame or list of data.frames
#'
#' @examples NA
#'
#' @export
ensListTransform_mod <- function(dataList, popsData = F) {
  # dataList is a list with 2 elements, dataList[[1]] = GWAS data table ; dataList[[2]] = Ensembl API data in R list form

  if(!requireNamespace("GWASpops.pheno2geno", quietly = TRUE)){
    library(GWASpops.pheno2geno)
  }
  Populations <- Populations
  #CONT <- purrr::flatten(dataList[[2]]) #removing nested structure such that all sublists are combined into one list within dataList[[2]]

  CONT <- dataList[[2]]
  CONT <- purrr::compact(CONT) # removing empty elements introduced by:
  ## multiAPIcall_variants2 (?)... I think its one of the for loops that are fixing these data elements: EnsVar_synonyms and EnsVar_Clinical_significance.

  GWAS_DF <- dataList[[1]] #storing GWAS data from GWAS files for later.. (similar to createMT())
  data.table::setnames(GWAS_DF, old = 'SNPS', new = 'VariantID')

  if(popsData){
    # grabbing population data and converting into a list of tibbles.
    popData <- sapply(CONT, function(x) x$populations) #OPTIMIZATION: this may be more efficient with masking.. not sure though
    popData <- lapply(popData, function(x) dplyr::bind_rows(x)) # OPTIMIZATION: check if this can run without the anonymous function in lapply() .. I imagine its increasing operations for this call.

    # removes populations from the response content so further operations proceed properly.
    CONT <- lapply(CONT, function(x) x[names(x) != 'populations']) # OPTIMIZATION: Check for function which removes and returns elements from lists... as this call may removed if the original popData <- sapply() call removed and returned
  }

  # removing multimapping by flattening the lists out. (some rsIDs posses multiple mappings against the reference genome(?) or against different data within Ensembl's API databases(?) )
  CONT <- GWASpops.pheno2geno:::fixMultiMapping(CONT)
  CONT <- CONT[!sapply(CONT, is.null)] # this is a quick and dirty solution to the fact that fixMultiMapping() is producing null list entries at the end of its list output. I don't know why this is happening. OPTIMIZATION: DEBUG THE ISSUE MENTIONED IN THIS LINE FOR fixMultiMapping()  .... OPTIMIZATION 2: look to the comment below about $failed mappings being introduced occassionally, check for them within fixMultiMapping if possible and remove the need for additional code out here.

  # infrequently a `$failed` key:value pair is being introduced into lists after flattening out mappings, this indicates that a mapping doesn't map to the reference genome in Ensembl's data base, thus we are removing such entries.
  hasFailed <- sapply(CONT, \(x) rlang::has_name(x, "failed")) # MAKES Boolean mask
  CONT <- CONT[!hasFailed] # USES Boolean mask to filter out entries with `failed` key:value pairs

  CONT <- lapply(CONT, GWASpops.pheno2geno:::null2NA_ENSvariants)

  CONT_Table <- GWASpops.pheno2geno:::rsTable(CONT) #CONT_Table at this point is just EnsVariants. No GWAS data or Pop data.

  #renaming cols so their source is evident in the master table.
  names(CONT_Table) <- paste0('EnsVar_',names(CONT_Table))

  if(popsData){
    # setting ancestral allele attribute on population frequency data.
    popData <- GWASpops.pheno2geno:::AncestralAllele_attr(CONT_Table, popData)
    masterList <- list(CONT_Table, popData)

    # Merging Ensembl variant and GWAS data tables
    masterTable <- tryCatch(
      expr = {
        masterTable <- data.table:::merge.data.table(GWAS_DF, CONT_Table, by.x = 'VariantID', by.y = 'EnsVar_name')
      },

      error = function(e){ # in the case of too many duplicate rows causing the merge to
        # fail initially this option will allow for the merge to proceed.
        masterTable <- data.table:::merge.data.table(GWAS_DF, CONT_Table, by.x = 'VariantID', by.y = 'EnsVar_name', allow.cartesian = T);
        message("merge.data.table performed with `allow.cartesian = TRUE`, therefore many extra rows may be produced. Duplicated rows have been removed");
        masterTable$EnsVar_synonyms <- as.character(masterTable$EnsVar_synonyms);
        masterTable <- masterTable[!duplicated(masterTable)]; #removing many duplicated rows created.
        return(masterTable)
      }
    )
    # Transforming data for single population based data tables
    singlePop_alleleFreqDTs <- lapply(Populations$Population_Abbreviation,
                                      function(x) GWASpops.pheno2geno:::singlePopTransform(masterList[[2]], targetPopulation = x))

    # Populations is a data object that comes with the package. (see ?Populations for more information or inspect the object itself.)
    names(singlePop_alleleFreqDTs) <- Populations$Population_Abbreviation

    masterListFinal <- list(masterTable, masterList[[2]], singlePop_alleleFreqDTs)
    names(masterListFinal) <- c('masterTable', 'PopAlleleFreqData', 'singlePop_alleleFreqDTs')

    return(masterListFinal) ######## END for pops
  }

  #------------- only variant data ------------------------

  # Merging Ensembl variant and GWAS data tables
  masterTable <- tryCatch(
    expr = {
      masterTable <- data.table:::merge.data.table(GWAS_DF, CONT_Table, by.x = 'VariantID', by.y = 'EnsVar_name')
    },

    error = function(e){ # in the case of too many duplicate rows causing the merge to
      # fail initially this option will allow for the merge to proceed.
      masterTable <- data.table:::merge.data.table(GWAS_DF, CONT_Table, by.x = 'VariantID', by.y = 'EnsVar_name', allow.cartesian = T);
      message("merge.data.table performed with `allow.cartesian = TRUE`, therefore many extra rows may be produced. Duplicated rows have been removed");
      masterTable$EnsVar_synonyms <- as.character(masterTable$EnsVar_synonyms);
      masterTable <- masterTable[!duplicated(masterTable)]; #removing many duplicated rows created.
      return(masterTable)
    }
  )
  return(masterTable) ######### END for vars
}



# -------------------------------------------------------------------------





































































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







































