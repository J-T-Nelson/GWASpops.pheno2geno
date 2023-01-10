# AUTHOR: Jon Tanner Nelson
# LAST UPDATED: 10-3-2022
# CONTENTS: createMT()
# PURPOSE: createMT is the central all-in-one function for the GWAS-p2g pipeline. It calls upon many helper functions in order to generate organized labeled data structures in a list which can be further processed for visualization or data storage purposes.


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
#' @param processData TRUE by defualt. When TRUE data returned will be a list of flat tables. When FALSE data returned will be nested lists which are derived from the JSON objects returned by calling the REST API with `get_ensVariants()`.
#'
#' @returns A list of data.frames OR a nested list of lists (depending on param: processData).
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
      raw_data <- get_ensVariants(uniqueVariantIDs, population_data = TRUE)

    } else {
      if(varAnnotations){
        raw_data <- get_ensVariants(uniqueVariantIDs)

      } else {
        return(GWAS_DF)
      }
    }
    MTandRawData <- list(GWAS_DF, raw_data)
    return(MTandRawData)
  }


  # calling Ensembl API for both variant and population allele frequency data
  if(varAnnotations && population_data){
    varPop_dataLists <- get_ensVariants(uniqueVariantIDs, population_data = TRUE)
    cat('API calling complete. Data Transformation beginning.\n')
    allData <- list(GWAS_DF, varPop_dataLists)
    masterList <- ensListTransform(allData, popsData = T)

    return(masterList)
  }

  # calling Ensembl API for only variant data
  if(varAnnotations){

    #create variant-annotation table and merge into master table
    variant_dataLists <- get_ensVariants(uniqueVariantIDs)
    cat('API calling complete. Data Transformation begining.')
    allData <- list(GWAS_DF, variant_dataLists)
    masterTable <- ensListTransform(allData, popsData = F)

  } else {
    masterTable <- GWAS_DF # if variant annotation isn't desired, the master table IS the GWAS data.frame
  }

  return(masterTable)
}























