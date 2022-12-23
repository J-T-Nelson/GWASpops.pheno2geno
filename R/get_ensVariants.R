# AUTHOR: Jon Tanner Nelson
# LAST UPDATED: 10-3-2022
# CONTENTS: get_ensVariants(), EnsVarList2row(), AncestralAllele_attr(), rsList_asTable(), null2NA_ENSvariants(), fixMultiMapping(), rsTable()
# PURPOSE: primary pipeline tool for grabbing and transforming data from Ensembl Variants endpoint. All functions other than get_ensVariants() are helpers for get_ensVariants()


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


# HELPER FUNCS: -------------------------------------------------------------------------------

## multiAPIcall_variants() splits up large rsID lists and returns the expected objects, as the API only tolerates calls of ~150 rsIDs despite their documentation specifying 1000 per POST call

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

#' AncestralAllele_attr
#'
#' AncestralAllele_attr is a setter for a unique attribute necessary for graphing the data.
#'
#' Sets the attribute 'Ancestral_Allele' of population allele frequency tables for single variants by looking up the ancestral allele from the within the data retrieved from Ensembl. By attaching the ancestral alleles for each variant to each variant's table, graphing is made much easier, as ancestral alleles are highly important to understanding the meaning of the data, and thus should be availble when designing visualization of allele frequency data.
#'
#' @param masterTable master table comes from calling get_ensVariants(), must have variable: "EnsVar_ancestral_allele as ancestral allele data comes from Ensembl, not from GWAS direct downloads
#' @param popFreqList A list of population allele frequency data retrieved from Ensembl API. This list is grabbed by get_ensVariants().
#'
#' @example popData <- AncestralAllele_attr(CONT_Table, popData)
#'
#' @return popFreqList with new attribute added to each table within the list
#'
#' @noRd
AncestralAllele_attr <- function(masterTable, popFreqList){

  if('VariantID' %in% names(masterTable)){
    #if statement here allows for data to come straight from get_ensVariants()
    #or alternatively from a merged masterTable where EnsVar_name col has been overwritten with VariantID name
    uniqMasterT <- masterTable[!duplicated(masterTable$VariantID), ]

    for(variant in 1:length(popFreqList)){ # for each variant table in the population allele frequency list assign the attribute
      attr(popFreqList[[variant]], 'Ancestral_Allele') <-
        uniqMasterT[uniqMasterT$VariantID == names(popFreqList[variant]), ]$EnsVar_ancestral_allele

      attr(popFreqList[[variant]], 'VariantID') <- attr(popFreqList[variant], 'name') #adding name attribute to each table as well for graphing titles
    }
    return(popFreqList)

  }else{

    uniqMasterT <- masterTable[!duplicated(masterTable$EnsVar_name), ]

    uniqMasterT <- uniqMasterT[!is.na(uniqMasterT$EnsVar_name), ] # Removing NA rows that may be introduced, they crash later functions by having "Ancestral_Allele" attributes which contain 2 entries. (the second entry always being `NA`) #not sure if this was why doubles were introduced now...

    for(variant in 1:length(popFreqList)){
      # for this ancestral allele assignment, how are we sure we are getting the right allele when there are multiple mappings for a given rsID? ... Shouldn't we NOT use a unique master table and grab all ancestral alleles for a given variant in all its mappings and attach them to the ancestral allele attribute? (possibly with some meta data to point towards source of diff mappings ideally?)
      attr(popFreqList[[variant]], 'Ancestral_Allele') <-
        uniqMasterT[uniqMasterT$EnsVar_name == names(popFreqList[variant]), ]$EnsVar_ancestral_allele

      #DEBUG CODE FOR TEMP USE: IF THIS WORKS ADD IT TO THE ABOVE ALT PATH FOR THIS FUNCTION AS WELL
      if(length(attr(popFreqList[[variant]], 'Ancestral_Allele')) == 0){ #eliminating empty entries which cause issues in transformations.
        attr(popFreqList[[variant]], 'Ancestral_Allele') <- NA
      }

      attr(popFreqList[[variant]], 'VariantID') <- attr(popFreqList[variant], 'name') #adding name attribute to each table as well for graphing titles
    }
    return(popFreqList)
  }
}


#' null2NA_ENSvariants
#'
#' Converts NULL values to NA values within get_ensVariants() function.
#'
#' This function is strictly a helper for get_ensVariants() and is not capable of more generalized use. This function simply goes through the known levels of the list at a late stage of processing data retrieved from Ensembl API calls and fixes the data by converting NULL values to NA.
#'
#' @param rsOBJ look at get_ensVariants() to verify desired input.
#'
#' @example N/A
#'
#' @return N/A
#'
#' @noRd
null2NA_ENSvariants <- function(rsOBJ){
  #converts NULL values to NA ...
  #NULL vals make many object manipulations (data transformations) impossible and thus must be culled

    #using logical mask to grab any keys which have NULL values and then assigning them as NA instead
  rsOBJ[as.logical(lapply(rsOBJ, is.null))] <- NA

    # logical masking to replace NULL values with NA within the "mappings" list which is within a given rsID list object
  rsOBJ[["mappings"]][[1]][ sapply(rsOBJ[["mappings"]][[1]], is.null) ] <- NA

  return(rsOBJ)
}


#' ENsVarList2row
#'
#' Converts response objects (lists) from Ensembl API calls into single rows which can be bound into a data.frame like data structure
#'
#' Significant data manipulation / fixing is necessary to convert data returned from Ensembl "Variation" endpoint to rows from their natural list format. This function handles the various inconsistencies and structural issues which come from the lists returned by Ensembl REST API calls.
#'
#' @param rsO_list rs Object list, a list which is returned by converting JSON data retrieved from Ensembl Variation endpoint via jsonlite::fromJSON()
#'
#' @example row <- EnsVarList2row(rsObjList[[i]]) ## see rsTable() for greater context. (it is within this script!)
#'
#' @return tibble / data.frame
#'
#' @importFrom tibble as_tibble
#' @importFrom purrr pluck
#' @importFrom stringr str_flatten
#' @import magrittr
#' @importFrom dplyr select
#'
#' @noRd
EnsVarList2row <- function(rsO_list){

   if (is.null(rsO_list$clinical_significance)){ # (DATA FIXING) not many IDs have clinical_significance, but those that do mess up table creation due to the extra column.. thus we add in NA values when none are found
    rsO_list$clinical_significance <- NA
  }

  if (length(rsO_list$synonyms) == 0){ # (DATA FIXING) reassigning a value so the tibble is consistent with other tibbles
    rsO_list$synonyms <- 'NONE'
  }

  if (length(rsO_list$synonyms) > 1 || length(rsO_list$clinical_significance) > 1){ # converts many synonyms into a single string for convenient tabulation of data

    # CRITICAL: removing all elements which can be nested lists in the data set such that as_tibble() creates a single row.
    list_noSyn <- rsO_list[-which(names(rsO_list) %in% c("synonyms", "evidence", "mappings", "clinical_significance"))]

    rsTib <- as_tibble(list_noSyn)

    # cleaning disorderly cols before rebinding into a tibble,
    # in this case synonyms list is also condensed into a single string just like evidence
    evidence <- rsO_list %>% pluck('evidence') %>% as.character() %>% str_flatten(collapse = "|")
    synonyms <- rsO_list %>% pluck('synonyms') %>% as.character() %>% str_flatten(collapse = "|")

    mappings <- as.data.frame(rsO_list$mappings)

    if(length(rsO_list$clinical_significance) > 1){
      clinical_significance <- rsO_list %>% pluck('clinical_significance') %>% as.character() %>% str_flatten(collapse = "|")
    return(as_tibble(cbind(evidence, synonyms, rsTib[1,], mappings, clinical_significance)))
    }

    clinical_significance <- as.character(rsO_list$clinical_significance)
    return(as_tibble(cbind(evidence, synonyms, rsTib[1,], mappings, clinical_significance)))

  } else {

    rsTib <- as_tibble(rsO_list)
    cleanedTBL <- select(rsTib, -c(evidence, mappings))

    # cleaning disorderly cols before rebinding into a tibble
    evidence <- rsO_list %>% pluck('evidence') %>% as.character() %>% str_flatten(collapse = "|")
    mappings <- as.data.frame(rsO_list$mappings)

    return(as_tibble(cbind(evidence, cleanedTBL[1,], mappings)))
  }
}


#' rsTable
#'
#' stupid little function binding rows into a table
#'
#' function should be replaced by bind_rows in get_ensVariants() at some point for efficency and clarity of code
#'
#' @param rsObjList N/A
#'
#' @example N/A
#'
#' @return N/A
#'
#' @noRd
rsTable <- function(rsObjList){

  table <- data.frame()

  for(i in 1:length(rsObjList)){
    row <- EnsVarList2row(rsObjList[[i]])
    table <- rbind(table, row)
  }
  return(table)
}


# Multiple Mappings Script ------------------------------------------------
#
#  having multiple mappings for any rsID creates issues in transforming API returned objects into rows
#  Each mapping will correspond to its own row and this script is going to make a new list item for each mapping beyond the first for and rsIDobject

#' fixMultiMapping
#'
#' helper function to get_ensVariants() which finds objects within the list returned by calling Ensembl Variants endpoint that posses multiple mappings and converts each mapping into its own object which allows for flattening and tabulation of data later in processing of data.
#'
#' N/A
#'
#' @param rsO_list rs Object List (variant Object List) look to get_ensVariants() to see the origin of objects meant for this func.
#'
#' @example N/A
#'
#' @return list of lists
#'
#' @noRd
fixMultiMapping <- function(rsO_list){
  #check each rsObj in list for mapping list > 1
  #create new Obj for each additional mapping and add them to the list

  listLen <- length(rsO_list)

  # creating vector containing mappings per object in rsO_list
  mappings <- numeric(listLen)
  for (i in 1:listLen){
    mappings[i] <- length(rsO_list[[i]][['mappings']])
  }

  # if there are multiple mappings anwhere, make them into new objects in a list of upadated size
  # if not, return the original object

  if(sum(mappings) > listLen){

    split_rsO_list <- vector(mode = 'list', length = sum(mappings))
    multiMappingIndices <- which(mappings[]>1)

    listPosition <- 1
    for (index in 1:length(multiMappingIndices)){
      #for each index with multimapping... for each mapping at that in that rsOBJ with multimapping ... make an rsOBJ in the newList with only mapping = forL2-index in that position
      innerLoopCounter <- 1
      for(mapping in 1:length(rsO_list[multiMappingIndices[index]][['mappings']])){

        split_rsO_list[ listPosition ] = rsO_list[multiMappingIndices[index]]
        split_rsO_list[[listPosition]][['mappings']] <- rsO_list[[multiMappingIndices[index]]][['mappings']][innerLoopCounter]

        innerLoopCounter = innerLoopCounter + 1 #innerLoopCounter is in use because the value of 'mapping' was decrementing mysteriously

        listPosition = listPosition + 1
      }
    }

    rsO_listNEW <- rsO_list[-multiMappingIndices]

    for(item in 1:length(rsO_listNEW)){ #adding non-multimapped rsOBJs back into the new list
      split_rsO_list[ listPosition ] <- rsO_listNEW[item]
      listPosition = listPosition + 1
    }
    return(split_rsO_list)

  } else {
      return(rsO_list)
    }
}


#' ensemblPing
#'
#' pings ensembl REST API in order to verify that it is currently operational
#'
#' Useful for constructing functions which call upon the REST API, as this will help narrow down where any errors come from
#'
#' @return 1 or 0 depending on response of ensembl servers.
#' A message detailing the success or failure of ping is also printed to console.
#'
#' @example
#' ensemblePing()
#'
#' @importFrom httr stop_for_status
#' @importFrom httr GET
#' @importFrom httr content_type
#' @importFrom httr content
#'
#' @noRd
ensemblPing <- function() {

  server <- "https://rest.ensembl.org"
  ext <- "/info/ping?"

  request <- GET(paste(server, ext, sep = ""), content_type("application/json"))

  stop_for_status(request)


  result <- content(request)
  if(result == 1){
    cat("ensembl API Ping successful\n")
    return(1)
  } else {
    cat("ensembl API ping failed\n")
    return(0)
  }
}
