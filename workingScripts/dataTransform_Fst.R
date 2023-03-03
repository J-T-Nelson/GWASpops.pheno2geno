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

# transform_fst_save() takes in the large GWAS association table, a number of chunks, and starting chunk and transforms the lists within the chunks into a list of DFs and lists OF DFs.. 4 objects specified above. Fst will be calculated within then the object will be saved.
transform_fst_save <- function(GWAS_associations,
                               numChunks,
                               startChunk = 1,
                               Fst_populations,
                               return_DS = FALSE,
                               saveData = TRUE){

  # load data from memory and flatten for processing
  variantList <- load_n_flatten(numChunks = numChunks, startSuperChunk = startChunk)

  # compose list then tranform into GWASpops.geno2pheno masterList format
  dataList <- list(GWAS_associations, variantList)
  masterList <- ensListTransform_mod(dataList, TRUE)

  # calculate Fst, delete redudant data vals, and discard multiallelic sites
  fstList <- hudsonFst_alleleList(masterList[[2]], Fst_populations, deleteRedundants = TRUE, discardMultiAllelic =  TRUE)

  # make single table of Fst Value list, then bind to masterList data structure

  fstList <- fill_rows(fstList) # making all sublists compatible for binding together as data.frame
  names <- names(fstList)
  fstDF <- cbind.data.frame(fstList)
  colnames(fstDF) <- names
  fstDF <- as.data.frame(t(fstDF)) # transpose s.t. rows are alleles, cols are population-pairs

  masterList[['Fst_per_allele']] <- fstDF

  # save new data structure in memory
  if(saveData){
    setwd("D:\\Programming\\R_projects\\Kulathinal_Lab\\GWASpops.pheno2geno\\workingData\\fst_GWAS_annotation_lists")
    fileName <- paste0('fullData_', numChunks, '_', startChunk)
    save(masterList, file = fileName)
    setwd("D:\\Programming\\R_projects\\Kulathinal_Lab\\GWASpops.pheno2geno")
  }

  # return nothing if desired, or resulting data structure if desired.
  if(return_DS){
    return(masterList)
  } else{
    return()
  }
}



# load_n_flatten() --------------------------------------------------------

# realized an error in my processing of the data.. I actually didn't need to rename anything! as I can assign names to the loaded objects adhoc within a list as loading them in... so I will have to manage the first 13 chunks manually somehow... maybe rename manually before processing through here.

load_n_flatten <- function(numChunks, startSuperChunk = 1) {

  setwd("D:\\Programming\\R_projects\\Kulathinal_Lab\\GWASpops.pheno2geno\\workingData\\unprocessedChunks")
  ret_list <- list()

  startPoint <- (startSuperChunk - 1)*100 + 1

  for(i in 1:numChunks){
    s <- startPoint + (i-1)*100
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

# fillRows func dev ----------------------------------------------------------------

# I need to fill the rows for missing values of population pairs in variants..
#
# making generalized func that can fill in missing rows for DFs with named rows

# chatGPT wrote this version, its been edited to actually work by me:
fill_rows <- function(DF_list){

  # find largest row
  max_rows <- max(sapply(DF_list, nrow))

  # get names of largest row as rowNames
  rowNames <- row.names(DF_list[[which.max(sapply(DF_list, nrow))]])

  # for each DF in DF_list add missing rows with NA filled in using rowNames
  for (i in seq_along(DF_list)) {
    if (nrow(DF_list[[i]]) < max_rows) {
      missing_rows <- data.frame(matrix(NA, nrow = max_rows - nrow(DF_list[[i]]), ncol = ncol(DF_list[[i]])))
      colnames(missing_rows) <- "Fst_Hudson"
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

  if(is.null(GWAS_DF[['VariantID']])){ # renaming col for compatibility of pipeline functions
    data.table::setnames(GWAS_DF, old = 'SNPS', new = 'VariantID')

  }

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



# Disorganized working notes:  --------------------------------------------

# ISSUE: Seeing a row with 515 unique names.. meaning there must be some names where the order is swapped making for 2x possible unique names... need to ensure this cannot happen in code that generates rownames. 496 is the max number that CAN exist for the 32C2 set of name combinations

# START HERE 2-26 ---------------------------------------------------------
#  working towards complete data transformation and fst calculation for all data thus far...
#  finished calling for data, ignoring missing data for now (probably throughout the rest of the proj.)
#  need to make formula / func for Wrights fst.. reading from wiki is a bit difficult but with a calm mind and some determination I am sure I can work something out
#  need to finish test run of the transformation func 'transform_fst_save()' ... currently hung up on issue detailed just above ... need to ensure that I cannot produce more than the true max number of rows when calculating Fst... seems a bit of a tough one to crack.. as it was only produced once in the test data set I think, and beyond that I am wondering why so many were exactly the right number if this one somehow slipped through the cracks. .. will need a bit of brainstroming for debug strategies.
#  After the first test run of commands works I need to build out transform_fst_save() and run it against the first 24th of the data... see how long it takes.. and see if transformation fails at any point.. do debugging... then hopefully just run the rest of it through and finally see if all the data can be stored in one huge list with the 4 items of data forms desired.







































































































