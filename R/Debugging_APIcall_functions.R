# this script is to ease the reboot process, such that I can access these updated functions without manually running each one


#' get_ensVariants2
#'
#' @description debug func
#'
#' @details debug func
#'
#'
#' @param rsIDs vector of rsIDs (variant IDs of the rs00000000 form)
#' @param population_data when TRUE, activates the option to grab population data for each variant in the request as well.
#'
#' @return list of data.frame(s) / data.table(s)
#'
#' @examples
#' variantData <- get_ensVariants(c("rs11137048","rs6866110", "rs62227671", "rs6122625", "rs57504074"), population_data = TRUE)
#'
#'
#' @export
get_ensVariants2 <- function(rsIDs, population_data = FALSE, processData = TRUE){
  #Returns a table of variant annotations from Ensembl POST variants API endpoint.
  #If population_data is TRUE and thus requested, a list is instead returned
  # due to incompatible data formats. The list contains the variant annotation table as well as a
  # list of the population data tables. There is 1 population table per rsID entered.
  cat("Pinging from get_ensVaraiants2():  ")
  if(GWASpops.pheno2geno:::ensemblPing() == 0){ # check if service is up
    stop("terminating function early, ping unsuccessful\n")
  }

  if(length(rsIDs) > 101){
    # multiAPIcall_variants() splits up large rsID lists and returns the expected objects, as the API only tolerates calls of ~150 rsIDs despite their documentation specifying 1000 per POST call

    if(population_data){
      masterCONT <- GWASpops.pheno2geno:::multiAPIcall_variants2(rsIDs, popData = TRUE, procsData = processData)
    } else{
      masterCONT <- GWASpops.pheno2geno:::multiAPIcall_variants2(rsIDs, procsData = processData)
    }
    return(masterCONT)
  }

  baseURL <- "https://rest.ensembl.org/variation/homo_sapiens"

  if(population_data){
    baseURL <- "https://rest.ensembl.org/variation/homo_sapiens?pops=1"
  }

  rsID_Array <- paste('{ "ids" : [', paste(shQuote(rsIDs, type="cmd"), collapse=", "), "] }", sep = "")

  response <- POST(baseURL, content_type("application/json"), accept("application/json"), body = rsID_Array)

  stop_for_status(response)

  if (http_type(response) != "application/json"){
    stop("API did not return json", call. = FALSE)
  }

  CONT <- content(response)

  # returns list form of data instead of processed data.
  if(!processData){
    return(CONT)
  }

  if(population_data){
    # grabbing population data and converting into a list of tibbles.
    popData <- sapply(CONT, function(x) x$populations)
    popData <- lapply(popData, function(x) bind_rows(x))

    # removes populations from the response content so further operations proceed properly.
    CONT <- lapply(CONT, function(x) x[names(x) != 'populations'])
  }


  CONT_noMultiMapping <- GWASpops.pheno2geno:::fixMultiMapping(CONT)

  CONT_noNULL <- lapply(CONT_noMultiMapping, null2NA_ENSvariants)

  CONT_Table <- GWASpops.pheno2geno:::rsTable(CONT_noNULL)

  #renaming cols so their source is evident in the master table.
  names(CONT_Table) <- paste0('EnsVar_',names(CONT_Table))

  if(population_data){
    # setting ancestral allele attribute on population frequency data.
    popData <- GWASpops.pheno2geno:::AncestralAllele_attr(CONT_Table, popData)
    masterList <- list(CONT_Table, popData)
    return(masterList)
  }

  return(CONT_Table)
}


#' createMT2
#'
#' @description debug func
#'
#' @details debug func
#'
#'
#' @param fileFolderPath NA
#' @param varAnnotations NA
#'
#' @return NA
#'
#' @examples
#' NA
#'
#'
#' @export
createMT2 <- function(fileFolderPath,
                      varAnnotations = TRUE,
                      population_data = FALSE,
                      processData = TRUE){

  #importing GWAS data and smashing into single data.frame
  GWAS_DF_list <- GWASpops.pheno2geno:::importGWAS_DataTables(fileFolderPath)
  GWAS_DF <- GWASpops.pheno2geno:::list2table_associations_studies(GWAS_DF_list)
  uniqueVariantIDs <- unique(GWAS_DF$VariantID) #calling API with repeated IDs is a waste of time as it will be returning the same data multiple times.

  if(!varAnnotations && population_data){ #simplifying user experience by not allowing invalid input and informing them about invalid input.
    warning("varAnnotations must be True to retrieve population data\nSetting varAnnotations to TRUE and proceeding to retreive data.")
    varAnnotations = TRUE
  }

  # allows createMT to grab untransformed data from Ensembl API (which are just lists)
  if(!processData){
    if(varAnnotations && population_data){
      raw_data <- get_ensVariants2(uniqueVariantIDs, population_data = TRUE, processData = processData)

    } else {
      if(varAnnotations){
        raw_data <- get_ensVariants2(uniqueVariantIDs, processData = processData)

      } else { return(GWAS_DF) }
    }

    MTandRawData <- list(GWAS_DF, raw_data)
    return(MTandRawData)
  }

  # calling Ensembl API for both variant and population allele frequency data
  if(varAnnotations && population_data){
    var_pop_list <- GWASpops.pheno2geno:::get_ensVariants(uniqueVariantIDs, population_data = TRUE)
    masterTable <- merge(GWAS_DF, var_pop_list[[1]], by.x = 'VariantID', by.y = 'EnsVar_name')

    # Transforming data for single population based data tables
    singlePop_alleleFreqDTs <- lapply(Populations$Population_Abbreviation,
                                      function(x) GWASpops.pheno2geno:::singlePopTransform(var_pop_list[[2]], targetPopulation = x))

    # Populations is a data object that comes with the package. (see ?Populations for more information or inspect the object itself.)
    names(singlePop_alleleFreqDTs) <- Populations$Population_Abbreviation

    masterList <- list(masterTable, var_pop_list[[2]], singlePop_alleleFreqDTs)
    names(masterList) <- c('masterTable', 'PopAlleleFreqData', 'singlePop_alleleFreqDTs')

    return(masterList)
  }

  # calling Ensembl API for only variant data
  if(varAnnotations){

    #create variant-annotation table and merge into master table
    variant_Anno_Table <- GWASpops.pheno2geno:::get_ensVariants(uniqueVariantIDs)
    masterTable <- merge(GWAS_DF, variant_Anno_Table, by.x = 'VariantID', by.y = 'EnsVar_name')

  } else{
    masterTable <- GWAS_DF #if annotation isn't desired, the master table IS the GWAS data.frame
  }

  return(masterTable)

}


#' multiAPIcall_variants2
#'
#' @description debug func
#'
#' @details debug func
#'
#'
#' @param fileFolderPath NA
#' @param varAnnotations NA
#'
#' @return NA
#'
#' @examples
#' NA
#'
#'
#' @export
multiAPIcall_variants2 <- function(rsIDs, popData = FALSE, procsData = TRUE) {

  # splits the vector of rsIDs into sub-arrays of length 100, all sub-arrays are held in a list
  splitList <- maxVecLength(rsIDs, 100)
  holder <- as.list(vector(length = length(splitList)))

  #the for() loops below are where the API is repeatedly called.
  if(popData){

    if(!procsData){
      for(i in 1:length(splitList)){
        holder[[i]] <- get_ensVariants2(splitList[[i]], population_data = TRUE, processData = procsData)
      }

      return(holder)

    } else {

      for(i in 1:length(splitList)){
        holder[[i]] <- get_ensVariants2(splitList[[i]], population_data = TRUE, processData = procsData)

        #FIXING DATA FOR BINDING LATER ... not all synonyms or clin_sig come out as lists or chars
        holder[[i]][[1]]$EnsVar_synonyms <- as.character(holder[[i]][[1]]$EnsVar_synonyms)
        holder[[i]][[1]]$EnsVar_clinical_significance <- as.character(holder[[i]][[1]]$EnsVar_clinical_significance)
      }

      #code block below is extracting the tables and population allele frequency lists into separate objects and flattening the results of the multiple calls before returning a list of both a masterTable and a popFreqList
      varTableList <- as.list(vector(length = length(holder)))
      populationList <- as.list(vector(length = length(holder)))
      for(j in 1:length(holder)){
        varTableList[[j]] <- holder[[j]][[1]]
        populationList[[j]] <- holder[[j]][[2]]
      }
      varTable <- bind_rows(varTableList)
      populationList <- purrr::flatten(populationList)
      masterList <- list(varTable, populationList)

      return(masterList)
    }
  } else {

    for(i in 1:length(splitList)){
      holder[[i]] <- get_ensVariants2(splitList[[i]], processData = procsData)

      #FIXING DATA FOR BINDING LATER ... not all synonyms or clin_sig come out as lists or chars
      holder[[i]]$EnsVar_synonyms <- as.character(holder[[i]]$EnsVar_synonyms)
      holder[[i]]$EnsVar_clinical_significance <- as.character(holder[[i]]$EnsVar_clinical_significance)
    }

    # if data is unprocessed then: return the holder list without binding rows together, as rows will have inconsistencies in their variables (columns) which causes failure of the script.
    if(!procsData){
      return(holder)
    } else {
      varAnnotationTable <- bind_rows(holder)

      return(varAnnotationTable)
    }
  }
}


#' dbugTransform2
#'
#' @description debug func
#'
#' @details debug func
#'
#'
#' @param fileFolderPath NA
#' @param varAnnotations NA
#'
#' @return NA
#'
#' @examples
#' NA
#'
#'
#' @export
dbugTransform <- function(dataList, popsData = F ) {

  CONT <- purrr::flatten(dataList[[2]]) #removing nested structure to have 1 list containing all rsID objects
  CONT <- purrr::compact(CONT) # removing empty elements introduced by:
  ## multiAPIcall_variants2 ... think its one of the for loops that are fixing these data elements: EnsVar_synonyms and EnsVar_Clinical_significance

  GWAS_DF <- dataList[[1]] #storing GWAS data from GWAS files for later.. (similar to createMT())

  if(popsData){
    # grabbing population data and converting into a list of tibbles.
    popData <- sapply(CONT, function(x) x$populations)
    popData <- lapply(popData, function(x) bind_rows(x))

    # removes populations from the response content so further operations proceed properly.
    CONT <- lapply(CONT, function(x) x[names(x) != 'populations'])
  }


  CONT_noMultiMapping <- GWASpops.pheno2geno:::fixMultiMapping(CONT)

  CONT_noMultiMapping <- CONT_noMultiMapping[!sapply(CONT_noMultiMapping, is.null)] # this is a quick and dirty solution to the fact that fixMultiMapping() is producing null list entries at the end of its list output. I don't know why this is happening.

  CONT_noNULL <- lapply(CONT_noMultiMapping, GWASpops.pheno2geno:::null2NA_ENSvariants)

  CONT_Table <- GWASpops.pheno2geno:::rsTable(CONT_noNULL) #CONT_Table at this point is just EnsVariants. No GWAS data or Pop data.

  #renaming cols so their source is evident in the master table.
  names(CONT_Table) <- paste0('EnsVar_',names(CONT_Table))

  if(popsData){
    # setting ancestral allele attribute on population frequency data.
    popData <- GWASpops.pheno2geno:::AncestralAllele_attr(CONT_Table, popData)
    masterList <- list(CONT_Table, popData)
    #return(masterList)
  }

  #return(CONT_Table)

  #---------------createMT calls after here -------------------------

  if(popsData){
    #var_pop_list <- GWASpops.pheno2geno:::get_ensVariants(GWAS_DF$VariantID, population_data = TRUE)

     # merging Variant data with GWAS table data
    tryCatch(
      expr = {masterTable <- data.table:::merge.data.table(GWAS_DF, masterList[[1]], by.x = 'VariantID', by.y = 'EnsVar_name')},

      error = function(e){ # in the case of too many duplicate rows causing the merge to fail
                           #  initially this option will allow for the merge to proceed.
        masterTable <- data.table:::merge.data.table(GWAS_DF, masterList[[1]], by.x = 'VariantID', by.y = 'EnsVar_name', allow.cartesian = T);
        message("merge.data.table performed with `allow.cartesian = TRUE`, therefore many extra rows may
                be produced with many possibly being duplicated rows due to duplicates in the tables being merged.")
      }
    )
    #  return(masterList) # TEMP RETURN FOR DEBUGGING
    # Transforming data for single population based data tables
    singlePop_alleleFreqDTs <- lapply(Populations$Population_Abbreviation,
                                      function(x) GWASpops.pheno2geno:::singlePopTransform(masterList[[2]], targetPopulation = x))

    # Populations is a data object that comes with the package. (see ?Populations for more information or inspect the object itself.)
    names(singlePop_alleleFreqDTs) <- Populations$Population_Abbreviation

    masterListFinal <- list(masterTable, masterList[[2]], singlePop_alleleFreqDTs)
    names(masterListFinal) <- c('masterTable', 'PopAlleleFreqData', 'singlePop_alleleFreqDTs')

    return(masterListFinal)######## END for pops
  }

  # DEBUGGING RETURN:
  #   return(list(GWAS_DF, CONT_Table))
  # only variant data ------------------------
  # Merging Ensembl variant and GWAS tables
  tryCatch(
    expr = {masterTable <- data.table:::merge.data.table(GWAS_DF, CONT_Table, by.x = 'VariantID', by.y = 'EnsVar_name')},

    error = function(e){ # in the case of too many duplicate rows causing the merge to fail initially this option will allow for the merge to proceed.
              masterTable <- data.table:::merge.data.table(GWAS_DF, CONT_Table, by.x = 'VariantID', by.y = 'EnsVar_name', allow.cartesian = T);
              message("merge.data.table performed with `allow.cartesian = TRUE`, therefore many extra rows may be produced with many possibly being duplicated rows due to duplicates in the tables being merged.")
    }
  )

  return(masterTable) ######### END for vars
}

# masterMerge <- function(list){
#
#   GWAS_DF <- list[[1]]
#   CONT_Table <- list[[2]]
#
#   tryCatch(
#     expr = {masterTable <- merge.data.table(GWAS_DF, CONT_Table, by.x = 'VariantID', by.y = 'EnsVar_name')},
#     error = function(e){
#       masterTable <- merge.data.table(GWAS_DF, CONT_Table, by.x = 'VariantID', by.y = 'EnsVar_name', allow.cartesian = T)}
#
#   )
#   return(masterTable)
# }

