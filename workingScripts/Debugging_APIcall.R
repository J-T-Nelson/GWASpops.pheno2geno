# work space for debugging get_EnsVars() ...
# First need to attempt to grab the list data form of another data set and wait for it to fail to learn what error I am catching.




# rewriting func, will integrate into pipeline later.  --------------------



get_ensVariants2 <- function(rsIDs, population_data = FALSE, processData = TRUE){
  #Returns a table of variant annotations from Ensembl POST variants API endpoint.
  #If population_data is TRUE and thus requested, a list is instead returned
  # due to incompatible data formats. The list contains the variant annotation table as well as a
  # list of the population data tables. There is 1 population table per rsID entered.

  if(GWASpops.pheno2geno:::ensemblPing() == 0){ # check if service is up
    cat("terminating function early, ping unsuccessful\n")
    return(NULL)
  }

  if(length(rsIDs) > 101){
    # multiAPIcall_variants() splits up large rsID lists and returns the expected objects, as the API only tolerates calls of ~150 rsIDs despite their documentation specifying 1000 per POST call

    if(population_data){
      masterCONT <- multiAPIcall_variants2(rsIDs, popData = TRUE, )
    } else{
      masterCONT <- multiAPIcall_variants2(rsIDs)
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

  # returns list form of data instead of processed data.
  if(!processData){
    return(CONT)
  }

  CONT_noMultiMapping <- GWASpops.pheno2geno:::fixMultiMapping(CONT)

  CONT_noNULL <- lapply(CONT_noMultiMapping, null2NA_ENSvariants)

  CONT_Table <- GWASpops.pheno2geno:::rsTable(CONT_noNULL)

  #renaming cols so their source is evident in the master table.
  names(CONT_Table) <- paste0('EnsVar_',names(CONT_Table))

  if(population_data){
    # setting ancestral allele atrribute on population frequency data.
    popData <- GWASpops.pheno2geno:::AncestralAllele_attr(CONT_Table, popData)
    masterList <- list(CONT_Table, popData)
    return(masterList)
  }

  return(CONT_Table)
}





# modifying createMT()  ---------------------------------------------------

#WHEN YOU"RE DONE DEBUGGING FIX NAMES AND UPDATE MANUAL PAGES ...
#       change createMT2 -> createMT and get_ensVariants2 -> get_ensVaraints
createMT2 <- function(fileFolderPath,
                      varAnnotations = TRUE,
                      population_data = FALSE,
                      processData = TRUE){

  #importing GWAS data and smashing into single data.frame
  GWAS_DF_list <- GWASpops.pheno2geno:::importGWAS_DataTables(fileFolderPath)
  GWAS_DF <- GWASpops.pheno2geno:::list2table_associations_studies(GWAS_DF_list)

  if(!varAnnotations && population_data){ #simplifying user experience by not allowing invalid input and informing them about invalid input.
    warning("varAnnotations must be True to retrieve population data\nSetting varAnnotations to TRUE and proceeding to retreive data.")
    varAnnotations = TRUE
  }

  # allows createMT to grab untransformed data from Ensembl API (which are just lists)
  if(!processData){
    if(varAnnotations && population_data){
      raw_data <- get_ensVariants2(GWAS_DF$VariantID, population_data = TRUE, processData)

    } else {
        if(varAnnotations){
        raw_data <- get_ensVariants2(GWAS_DF$VariantID, processData)

      } else { return(GWAS_DF) }
    }

    MTandRawData <- list(GWAS_DF, raw_data)
    return(MTandRawData)
  }

  # calling Ensembl API for both variant and population allele frequency data
  if(varAnnotations && population_data){
    var_pop_list <- GWASpops.pheno2geno:::get_ensVariants(GWAS_DF$VariantID, population_data = TRUE)
    masterTable <- merge(GWAS_DF, var_pop_list[[1]], by.x = 'VariantID', by.y = 'EnsVar_name')

    # Transforming data for single population based data tables
    singlePop_alleleFreqDTs <- lapply(Populations$Population_Abbreviation,
                                      function(x) GWASpops.pheno2geno:::singlePopTransform(var_pop_list[[2]], targetPopulation = x))
    names(singlePop_alleleFreqDTs) <- Populations$Population_Abbreviation

    masterList <- list(masterTable, var_pop_list[[2]], singlePop_alleleFreqDTs)
    names(masterList) <- c('masterTable', 'PopAlleleFreqData', 'singlePop_alleleFreqDTs')

    return(masterList)
  }

  # calling Ensembl API for only variant data
  if(varAnnotations){

    #create variant-annotation table and merge into master table
    variant_Anno_Table <- GWASpops.pheno2geno:::get_ensVariants(GWAS_DF$VariantID)
    masterTable <- merge(GWAS_DF, variant_Anno_Table, by.x = 'VariantID', by.y = 'EnsVar_name')

  } else{
    masterTable <- GWAS_DF #if annotation isn't desired, the master table IS the GWAS data.frame
  }

  return(masterTable)

}



# multiAPIcall modificatied -----------------------------------------------


multiAPIcall_variants2 <- function(rsIDs, popData = FALSE, proccData = TRUE) {

  splitList <- maxVecLength(rsIDs, 100)
  holder <- as.list(vector(length = length(splitList)))

  if(popData){

    for(i in 1:length(splitList)){
      holder[[i]] <- get_ensVariants2(splitList[[i]], population_data = TRUE, processData = proccData)

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
      holder[[i]] <- get_ensVariants2(splitList[[i]], processData = proccData)

      #FIXING DATA FOR BINDING LATER ... not all synonyms or clin_sig come out as lists or chars
      holder[[i]]$EnsVar_synonyms <- as.character(holder[[i]]$EnsVar_synonyms)
      holder[[i]]$EnsVar_clinical_significance <- as.character(holder[[i]]$EnsVar_clinical_significance)
    }

    varAnnotationTable <- bind_rows(holder)

    return(varAnnotationTable)
  }
}



# testing raw data calls.  ------------------------------------------------

testData1 <- createMT2('./exampleData/air_pollution', processData = F)

tdecomp1 <- testData1[[1]]
tdecomp2 <- testData1[[2]]

testData2 <- createMT2('./exampleData/air_pollution', processData = F, population_data = T)

# Checking num obs (rsIDs for alcohol consumption) ... sufficiently large dataset. I anticipate failure. So, lets call the API with it and see what errors we get.
alcConsumpRawTable <- createMT2(GWASdataSets[1], processData = F, varAnnotations = F)

alcConsumpRawVars <- createMT2(GWASdataSets[1], processData = F)
alcConsumpRawPops <- createMT2(GWASdataSets[1], processData = F, population_data = T)

#while
