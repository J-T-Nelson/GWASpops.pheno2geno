# This doc is purely for easy storage and access to archived functions and code.
#
# I will store functioning code with comments that I don't




# -------------------------------------------------------------------------

#' createMT
#'
#' @description  createMT is the central all-in-one function for the GWAS-p2g pipeline. It calls upon many helper functions in order to generate organized, labeled data structures in a list or table (depending on which options are active)which can be further processed for visualization or data storage purposes.
#'
#' * "createMT" is short for "create Master Table"
#'
#' @details For population_data to function when TRUE varAnnotations must also be TRUE. This is because the data associated with the population_data argument is optional data which comes from the same endpoint that is called for varAnnotations = TRUE.
#'
#' GWAS tables which are to be processed and used by this function should be association or studies tables. Other table types are likely to break the function as it was specifically designed to take in those two table types and thus depends on their specific variables / values.
#'
#' This function will take substantially longer depending on the # of variant IDs found within the tables imported and when the population_data option is set TRUE. The major bottleneck is due to Ensembl's REST API, and thus cannot be worked around in some cases. varAnnotations = TRUE will also add a significant amount of run time to the function call due to the REST API bottle neck, though substantially less data is transferred when only varAnnotations is true and population_data is FALSE.
#'
#' This function has failed several times in a row in the past during testing due to failures of API calls, thus one may need patience to see the data they wish to see, as reliability / stability of the API being called is beyond control of the user.
#'
#'
#' @param fileFolderPath directory path where one or more GWAS tables are stored.
#' @param varAnnotations TRUE by default. When TRUE Ensembl Variation API endpoint is called to retrieve genetic variant data
#' @param population_data FALSE by default. When TRUE population allele frequency data will be retrieved in addition to basic variant data from variation endpoint
#'
#' @returns A list or data.frame
#'
#' * When varAnnotations & population_data are FALSE a data.frame is returned.
#' * When either or both arguments are TRUE instead, a list containing different data tables will be returned.
#'
#' @examples # Calling createMT with example data available at https://github.com/J-T-Nelson/GWASpops.pheno2geno
#' # Download the 'exampleData' folder from github, then relocate it to the
#' # working directory you're working in currently for this command to run.
#'
#'      #masterList <- createMT("./exampleData/air_pollution", population_data = TRUE)
#'
#' # The call above will take some time and thus serves best as an example
#' # of how to use the function.
#'
#' @importFrom data.table merge.data.table
#'
#' @export
createMT <- function(fileFolderPath,
                     varAnnotations = TRUE,
                     population_data = FALSE){

  #importing GWAS data and smashing into single data.frame
  GWAS_DF_list <- importGWAS_DataTables(fileFolderPath)
  GWAS_DF <- list2table_associations_studies(GWAS_DF_list)

  if(!varAnnotations && population_data){ #simplifying user experience by not allowing invalid input and informing them about invalid input.
    warning("varAnnotations must be True to retrieve population data\nSetting varAnnotations to TRUE and proceeding to retreive data.")
    varAnnotations = TRUE
  }

  # calling Ensembl API for both variant and population allele frequency data
  if(varAnnotations && population_data){
    var_pop_list <- get_ensVariants(GWAS_DF$VariantID, population_data = TRUE)
    masterTable <- merge(GWAS_DF, var_pop_list[[1]], by.x = 'VariantID', by.y = 'EnsVar_name')

    # Transforming data for single population based data tables
    singlePop_alleleFreqDTs <- lapply(Populations$Population_Abbreviation,
                                      function(x) singlePopTransform(var_pop_list[[2]], targetPopulation = x))
    names(singlePop_alleleFreqDTs) <- Populations$Population_Abbreviation

    masterList <- list(masterTable, var_pop_list[[2]], singlePop_alleleFreqDTs)
    names(masterList) <- c('masterTable', 'PopAlleleFreqData', 'singlePop_alleleFreqDTs')

    return(masterList)
  }

  # calling Ensembl API for only variant data
  if(varAnnotations){

    #create variant-annotation table and merge into master table
    variant_Anno_Table <- get_ensVariants(GWAS_DF$VariantID)
    masterTable <- merge(GWAS_DF, variant_Anno_Table, by.x = 'VariantID', by.y = 'EnsVar_name')

  } else{
    masterTable <- GWAS_DF #if annotation isn't desired, the master table IS the GWAS data.frame
  }

  return(masterTable)

}



# -------------------------------------------------------------------------


#' multiAPIcall_variants
#'
#' A helper func to get_ensVariants() which allows larger numbers of rsIDs to be used successfully in calling get_ensVariants()
#'
#' Ensembl API cannot tolerate calls with > ~140 rsID (variant IDs), which is much less than they specify in documentation. Thus this function exists purely to automate the process of making multiple get_ensVariants() calls within a single call of the get_ensVariants().
#'
#' @param rsIDs vector containing variant ids as strings in the form 'rs000000000'
#' @param popData boolean value which inherits from the get_ensVariants() func calling this helper func. Changes behavior to account for extra data retrieved when popData is true.
#'
#' @example masterCONT <- multiAPIcall_variants(rsIDs, popData = TRUE)
#'
#' @return list of data.frame(s) / data.table(s); since this function calls get_ensVariants() it returns the exact same object types.
#'
#' @importFrom dplyr bind_rows
#' @importFrom purrr flatten
#'
#' @noRd
multiAPIcall_variants <- function(rsIDs, popData = FALSE) {

  splitList <- maxVecLength(rsIDs, 100)
  holder <- as.list(vector(length = length(splitList)))

  if(popData){

    for(i in 1:length(splitList)){
      holder[[i]] <- get_ensVariants(splitList[[i]], population_data = TRUE)

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


  } else {

    for(i in 1:length(splitList)){
      holder[[i]] <- get_ensVariants(splitList[[i]])

      #FIXING DATA FOR BINDING LATER ... not all synonyms or clin_sig come out as lists or chars
      holder[[i]]$EnsVar_synonyms <- as.character(holder[[i]]$EnsVar_synonyms)
      holder[[i]]$EnsVar_clinical_significance <- as.character(holder[[i]]$EnsVar_clinical_significance)
    }

    varAnnotationTable <- bind_rows(holder)

    return(varAnnotationTable)
  }
}



# -------------------------------------------------------------------------


#' get_ensVariants
#'
#' @description grabs ensVariants from API and converts to tabular form. Optionally grabs population data as well.
#'
#' @details Primary method to grab variant data from Ensembl REST API in the pipeline.
#'
#' Use this function if you're interested in grabbing data on specific variants from Ensembl's Variants endpoint.
#'
#' Calls the variants POST endpoint, build into the function are several layers of processing which both ensure arbitrarily large calls are executed successfully, as well as that data returned is in a set of tabular semi-flat formatted objects.
#'
#' Requesting population data will substantially increase the time this function takes to complete due to the increase in data being transferred.
#'
#' @param rsIDs vector of rsIDs (variant IDs of the rs00000000 form)
#' @param population_data when TRUE, activates the option to grab population data for each variant in the request as well.
#'
#' @return list of data.frame(s) / data.table(s)
#'
#' @examples
#' variantData <- get_ensVariants(c("rs11137048","rs6866110", "rs62227671", "rs6122625", "rs57504074"), population_data = TRUE)
#'
#' @importFrom httr POST
#' @importFrom httr stop_for_status
#' @importFrom httr accept
#' @importFrom httr http_type
#' @importFrom httr content_type
#' @importFrom httr content
#' @importFrom dplyr bind_rows
#'
#' @export
get_ensVariants <- function(rsIDs, population_data = FALSE){
  #Returns a table of variant annotations from Ensembl POST variants API endpoint.
  #If population_data is TRUE and thus requested, a list is instead returned
  # due to incompatible data formats. The list contains the variant annotation table as well as a
  # list of the population data tables. There is 1 population table per rsID entered.

  if(ensemblPing() == 0){ # check if service is up
    cat("terminating function early, ping unsuccessful\n")
    return(NULL)
  }

  if(length(rsIDs) > 101){
    # multiAPIcall_variants() splits up large rsID lists and returns the expected objects, as the API only tolerates calls of ~150 rsIDs despite their documentation specifying 1000 per POST call

    if(population_data){
      masterCONT <- multiAPIcall_variants(rsIDs, popData = TRUE)
    } else{
      masterCONT <- multiAPIcall_variants(rsIDs)
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

  if(population_data){
    # grabbing population data and converting into a list of tibbles.
    popData <- sapply(CONT, function(x) x$populations)
    popData <- lapply(popData, function(x) bind_rows(x))

    CONT <- lapply(CONT, function(x) x[names(x) != 'populations'])
    # ^^ removes populations from the response content so further operations proceed properly.
  }

  CONT_noMultiMapping <- fixMultiMapping(CONT)

  CONT_noNULL <- lapply(CONT_noMultiMapping, null2NA_ENSvariants)

  CONT_Table <- rsTable(CONT_noNULL)

  #renaming cols so their source is evident in the master table.
  names(CONT_Table) <- paste0('EnsVar_',names(CONT_Table))

  if(population_data){
    # setting ancestral allele atrribute on population frequency data.
    popData <- AncestralAllele_attr(CONT_Table, popData)
    masterList <- list(CONT_Table, popData)
    return(masterList)
  }

  return(CONT_Table)
}
