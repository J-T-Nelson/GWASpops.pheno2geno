# AUTHOR: Jon Tanner Nelson
# LAST UPDATED: 10-3-2022
# CONTENTS: createMT()
# PURPOSE: createMT is the central all-in-one function for the GWAS-p2g pipeline. It calls upon many helper functions in order to generate organized labeled data structures in a list which can be further processed for visualization or data storage purposes.


#' createMT
#'
#' createMT is the central all-in-one function for the GWAS-p2g pipeline. It calls upon many helper functions in order to generate organized, labeled data structures in a list which can be further processed for visualization or data storage purposes.
#'
#' For population_data to function when TRUE varAnnotations must also be TRUE. This is because the data associated with the population_data argument is optional data which comes from the same endpoint that is called for varAnnotations = TRUE.
#' GWAS tables which are to be processed and used by this function should be association or studies tables. Other table types are likely to break the function as it was specifically designed to take in those two table types and thus depends on their specific variables / values.
#' This function will take substantially longer depending on the # of variant IDs found within the tables imported and when the population_data option is set TRUE. The major bottleneck is due to Ensembl's REST API, and thus cannot be worked around in some cases. varAnnotations = TRUE will also add a significant amount of run time to the function call due to the REST API bottle neck, though substantially less data is transferred when only varAnnotations is true and population_data is FALSE.
#' This function has failed several times in a row in the past during testing due to failures of API calls, thus one may need patience to see the data they wish to see, as reliability / stability of the API being called is beyond control of the user.
#'
#' @param fileFolderPath directory path where one or more GWAS tables are stored.
#' @param varAnnotations TRUE by default. When TRUE Ensembl Variation API endpoint is called to retrieve genetic variant data
#' @param population_data FALSE by default. When TRUE population allele frequency data will be retrieved in addition to basic variant data from variation endpoint
#'
#' @returns When varAnnotations & population_data are FALSE a data.frame is returned. When either or both arguments are TRUE instead a list containing different data tables will be returned.
#'
#' @example N/A
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























