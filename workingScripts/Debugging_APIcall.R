# work space for debugging get_EnsVars() ...
# First need to attempt to grab the list data form of another data set and wait for it to fail to learn what error I am catching.




# rewriting func, will integrate into pipeline later.  --------------------



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
      masterCONT <- multiAPIcall_variants2(rsIDs, popData = TRUE, procsData = processData)
    } else{
      masterCONT <- multiAPIcall_variants2(rsIDs, procsData = processData)
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
      raw_data <- get_ensVariants2(GWAS_DF$VariantID, population_data = TRUE, processData = processData)

    } else {
        if(varAnnotations){
        raw_data <- get_ensVariants2(GWAS_DF$VariantID, processData = processData)

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

     # Populations is a data object that comes with the package. (see ?Populations for more information or inspect the object itself.)
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



# testing raw data calls.  ------------------------------------------------



testData1 <- createMT2('./exampleData/air_pollution', processData = F)

tdecomp1 <- testData1[[1]]
tdecomp2 <- testData1[[2]]

testData2 <- createMT2('./exampleData/air_pollution', processData = F, population_data = T)

# Checking num obs (rsIDs for alcohol consumption) ... sufficiently large dataset. I anticipate failure. So, lets call the API with it and see what errors we get.
alcConsumpRawTable <- createMT2(GWASdataSets[1], processData = F, varAnnotations = F) #GWASdataSets is a part of our bootscripts now, just makes it a little easier to access the data sets programmatically

alcConsumpRawVars <- createMT2(GWASdataSets[1], processData = F)
alcConsumpRawPops <- createMT2(GWASdataSets[1], processData = F, population_data = T) # recycling error again, attempting to fix the functions for grabbing API data without fixing/transforming the data at all for population_data = TRUE

#while
?data_frame



# 11-13-2022 sesh ---------------------------------------------------------

# STRATEGY: ------------------------------
#
# I want to provoke API failure and see what error message is thrown, I suppose a better was to do this would be to look into the POST() function or whatever function is actually calling the API and to see what error messages can be produced... then I can begin implementing error catching into the functions more robustly.
# Once we know what the errors are we need to write the error catching code then find the biggest data to call the API for to provoke failure such that we can test for the error catching to work.. which may be a bit time consuming... will require some good console output to verify.
#
#  Once we have fixed the API calling failures we move onto debugging the data transformation steps by taking tables like 'alcConsumpRawTable' (unprocessed API data) and processing them with debugger to catch errors and fix code.
#  After debugging transformation for all data sets we integrate changes into the original functions, update the documentation for the those functions, and retest everything. (may be useful to learn unit testing at this point to speed up the testing process and build on that experience.)
#
#  FINALLY update tutorial or man pages such that the ability to grab un-transformed data is spoken about such that a user could transform data for themselves if they're interested or whatever they want with that data form.

alcConsumpRawTable <- createMT2('GWAS.data/alcohol_consumption', processData = F) # <- testing the APIs stability (attempting to)

?POST

createMT2('GWAS.data/alcohol_consumption', varAnnotations = F, processData = F)

# importing all tables into global environment for visual inspection. Want to determine which are most likely to cause error messages... (finding the bigger data sets)
allGWAStables <- list()
for(i in 1:10){ allGWAStables[[i]] <- createMT2(GWASdataSets[i], varAnnotations = F)}

alcConsumpAlldata <- createMT2('GWAS.data/alcohol_consumption', processData = F, population_data = T) # <- seeking API failure for error catching ....
# ^^ no API failure yet, even after calling in many times over several hours... a bit annoying. I think I will just generically write the API code in... Population data was discarded though? why?


airPollutionAlldata <- createMT2('exampleData/air_pollution', processData = F, population_data = T)


GWASdataSets

alcConsumpVars <- createMT2(GWASdataSets[1], processData = F)
alcConsumpAlldata <- createMT2(GWASdataSets[1], processData = F, population_data = T)

bCarcinomaVars <- createMT2(GWASdataSets[2], processData = F)
bCarcinomaAlldata <- createMT2(GWASdataSets[2], processData = F, population_data = T) #testing for failure against a slightly bigger data set. The population data is no longer being discarded.

colorectalCancerVars <- createMT2(GWASdataSets[3], processData = F)
colorectalCancerAlldata <- createMT2(GWASdataSets[3], processData = F, population_data = T)

IBFVars <- createMT2(GWASdataSets[4], processData = F)
IBFAlldata <- createMT2(GWASdataSets[4], processData = F, population_data = T)

IntVars <- createMT2(GWASdataSets[5], processData = F)
IntAlldata <- createMT2(GWASdataSets[5], processData = F, population_data = T)

lungCancerVars <- createMT2(GWASdataSets[6], processData = F)
lungCancerAlldata <- createMT2(GWASdataSets[6], processData = F, population_data = T)

malabsorptionSyndVars <- createMT2(GWASdataSets[7], processData = F)
malabsorptionSyndAlldata <- createMT2(GWASdataSets[7], processData = F, population_data = T)

neuroticismVars <- createMT2(GWASdataSets[8], processData = F)
neuroticismAlldata <- createMT2(GWASdataSets[8], processData = F, population_data = T)

prostateCancerVars <- createMT2(GWASdataSets[9], processData = F)
prostateCancerAlldata <- createMT2(GWASdataSets[9], processData = F, population_data = T)

substanceAbuseVars <- createMT2(GWASdataSets[10], processData = F)
substanceAbuseAlldata <- createMT2(GWASdataSets[10], processData = F, population_data = T)


AirPollutionVars <- createMT2('exampleData/air_pollution', processData = F)
AirPollutionAlldata <- createMT2('exampleData/air_pollution', processData = F, population_data = T)

# Saving data -------------------------------------------------------------

save(alcConsumpRawVars, file = 'data/debugging_raw_data/alcConsumpVars.rds')
save(alcConsumpVars, file = 'data/debugging_raw_data/alcConsumpVars.rds')
save(alcConsumpAlldata, file = 'data/debugging_raw_data/alcConsumpAlldata.rds')

save(bCarcinomaVars, file = 'data/debugging_raw_data/bCarcinomaVars.rds')
save(bCarcinomaAlldata, file = 'data/debugging_raw_data/bCarcinomaAlldata.rds')

save(colorectalCancerVars, file = 'data/debugging_raw_data/colorectalCancerVars.rds')
save(colorectalCancerAlldata, file = 'data/debugging_raw_data/colorectalCancerAlldata.rds')

save(IBFVars, file = 'data/debugging_raw_data/IBFVars.rds')
save(IBFAlldata, file = 'data/debugging_raw_data/IBFAlldata.rds')

save(IntVars, file = 'data/debugging_raw_data/IntVars.rds')
save(IntAlldata, file = 'data/debugging_raw_data/IntAlldata.rds')

save(lungCancerVars, file = 'data/debugging_raw_data/lungCancerVars.rds')
save(lungCancerAlldata, file = 'data/debugging_raw_data/lungCancerAlldata.rds')

save(malabsorptionSyndVars, file = 'data/debugging_raw_data/malabsorptionSyndVars.rds')
save(malabsorptionSyndAlldata, file = 'data/debugging_raw_data/malabsorptionSyndAlldata.rds')

save(neuroticismVars, file = 'data/debugging_raw_data/neuroticismVars.rds')
save(neuroticismAlldata, file = 'data/debugging_raw_data/neuroticismAlldata.rds')

save(prostateCancerVars, file = 'data/debugging_raw_data/prostateCancerVars.rds')
save(prostateCancerAlldata, file = 'data/debugging_raw_data/prostateCancerAlldata.rds')

save(substanceAbuseVars, file = 'data/debugging_raw_data/substanceAbuseVars.rds')
save(substanceAbuseAlldata, file = 'data/debugging_raw_data/substanceAbuseAlldata.rds')

save(AirPollutionVars, file = 'data/debugging_raw_data/AirPollutionVars.rds')
save(AirPollutionAlldata, file = 'data/debugging_raw_data/AirPollutionAlldata.rds')

# 11-17-2022 --------------------------------------------------------------
#  gonna save any unprocessed data I grab from the API as data objects.. Takes too long to make them fresh. Each will be saved first, then used for debugging.
#
#  While waiting for data to come in from the API I will make the code for processing them into useable tabulated data, and debug that transformation function until all 10 data sets are functional, from there I will integrate the changes into the actual scripts and test those scripts against the same data but from createMT()... at this point I will assume the pipeline will be robust to anything from GWAS catalog and only will do further debugging upon request from users.
#


# getEnsVars transformation code ------------------------------------------


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
  # setting ancestral allele atrribute on population frequency data.
  popData <- GWASpops.pheno2geno:::AncestralAllele_attr(CONT_Table, popData)
  masterList <- list(CONT_Table, popData)
  return(masterList)
}

return(CONT_Table)


# createMT transformation code --------------------------------------------


# calling Ensembl API for both variant and population allele frequency data
if(varAnnotations && population_data){
  var_pop_list <- GWASpops.pheno2geno:::get_ensVariants(GWAS_DF$VariantID, population_data = TRUE)
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
  variant_Anno_Table <- GWASpops.pheno2geno:::get_ensVariants(GWAS_DF$VariantID)
  masterTable <- merge(GWAS_DF, variant_Anno_Table, by.x = 'VariantID', by.y = 'EnsVar_name')

} else{
  masterTable <- GWAS_DF #if annotation isn't desired, the master table IS the GWAS data.frame, and it only contains data from the GWAS catalog.
}

return(masterTable)


# Constructing independent transformation stack function.  ----------------
#   Thinking all the transformations inside getEnsVars needs to happen first, then that output is put into further transformation within createMT() .... I wonder if at the higher architectural level it makes more sense to have all transformation occur within a single function, where get ensVars simply returns unprocessed data, and depending on what a user asks for when calling create MT the transforming function will be called appropriately. ...

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

  CONT_noNULL <- lapply(CONT_noMultiMapping, null2NA_ENSvariants)

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

    masterTable <- merge(GWAS_DF, masterList[[1]], by.x = 'VariantID', by.y = 'EnsVar_name')

    # Transforming data for single population based data tables
    singlePop_alleleFreqDTs <- lapply(Populations$Population_Abbreviation,
                                      function(x) GWASpops.pheno2geno:::singlePopTransform(masterList[[2]], targetPopulation = x))

    # Populations is a data object that comes with the package. (see ?Populations for more information or inspect the object itself.)
    names(singlePop_alleleFreqDTs) <- Populations$Population_Abbreviation

    masterListFinal <- list(masterTable, masterList[[2]], singlePop_alleleFreqDTs)
    names(masterListFinal) <- c('masterTable', 'PopAlleleFreqData', 'singlePop_alleleFreqDTs')

    return(masterListFinal)######## END for pops
  }

  # only variant data ------------------------
  # Merging Ensembl variant and GWAS tables
  masterTable <- merge(GWAS_DF, CONT_Table, by.x = 'VariantID', by.y = 'EnsVar_name')

    return(masterTable) ######### END for vars
}




# Testing dbugTransform() -------------------------------------------------
#
# NEED TO RELOAD DATA IN... I learned how to install my package after updating it.. all locally.
workHere <- function(HERE) {}


devtools::install("YOUR FUCKING PACKAGE LOCATION/Path") # <- this is how you locally update your package.... its annoying and I think you should unattach it first to not break anything like I did the first time.
unload('GWASpops.pheno2geno')
devtools::install("D:\\Programming\\R_projects\\Kulathinal_Lab\\GWASpops.pheno2geno") # even after unloading my package... it was not able to reload without restarting RStudio.. kind of annoying, but relatively trivial if following best practices.

load('data/debugging_raw_data/alcConsumpVars.rds') #loading my data

alcConsumpVarTransformed <- dbugTransform(alcConsumpRawVars, popsData = F)

debug(dbugTransform)
undebug(dbugTransform)

 # OK, there are too many nested levels in the list of rsID lists right now, which is coming from multiAPIcall... and thus we need to fix this nesting then retest.
 #


bCarc_flatten <- bCarcinomaVars[[2]]

bCarc_flatten2 <- purrr::flatten(bCarc_flatten) #perfect

alcConsumpVarFlatten <- purrr::flatten(alcConsumpRawVars[[2]])


l <- alcConsumpVarFlatten[[89]] # <- this thing is empty?
l2 <- alcConsumpVarFlatten[[1]]
l3 <- alcConsumpVarFlatten[[88]] #also empty? ... its an empty character vector?
l4 <- alcConsumpVarFlatten[[85]]

class(l)
names(alcConsumpVarFlatten) # GOOD. I am seeing that EnsVar_synonyms and EnsVar_Clinical_significance are both being added to the list here and there as empty character vectors for some reason.. I think the process of flattening is doing this.

# actually those two rows are being added systemically somehow... presumably from my createMT2() script... at the end of every sub list of rsID objects.
#
# multiAPIcall_variants2 is the source... think its one of the for loops that are fixing these data elements... ****going to go for a fast and dirty solution of removing all empty elements in a list.

alcConsumpVarFlatten2 <- purrr::compact(alcConsumpVarFlatten)
names(alcConsumpVarFlatten2) # ok good.



# NOTES FOR NOVEL FUNCTION USE: -------------------------------------------

# in order to check whether a list contains an element use purr::has_element. Further, refer to the KEY names as below, and refer to values by entering the values themselves in quotations. SO if MAF contains the number 56. entering "56" will return TRUE for purr::has_element(c, "56")..
# BY DEFAULT purrr::has_element(c, "MAF") will return false, as its searching in the value set, not the key set.
purrr::has_element(c, c$mappings)
#[1] FALSE
 View(c)
purrr::has_element(c, c$MAF)
#[1] TRUE




# 11-21-2022 .. continuing Debugging --------------------------------------
#    Currently still having some troubles transforming. There appears to be some name inconsistencies arising in the current data set, so I want to workout which variable is dropped or added somewhere by inspecting my data.

alcLists <- purrr::flatten(alcConsumpRawVars[[2]]) #removing nested structure to have 1 list containing all rsID objects
alcLists <- purrr::compact(alcLists)

alcFilt1 <- sapply(alcLists, length)
a <- alcFilt1 > 11
sum(a) # = 0 ... nothing has more than 11 variables, meaning that something must have less, or there could be the odd case where a common variable is dropped and a rare one is introduced... I can rule this out by getting all unique names from the whole thing

uniqueNames <- sapply(alcLists, names)
uniqueNames <- purrr::flatten_chr(uniqueNames)
uniqueNames <- unique(uniqueNames)
uniqueNames # only 11 total, so we now know there cannot be any odd cases as aforementioned

b <- alcFilt1 < 10
sum(b) # 0 here too... therefore... only 10 and 11 length for the whole thing.
unique(alcFilt1) # <- this is actually a better way to see whats going on. Only 10 and 11.
# ok. so the data seems pretty regualr tbh. I guess I will have to dig in with the debugger to see where mismatching names is coming up.


#"\"risk factor\""         "allele_string"           "ambiguity"               "ancestral_allele"        "assembly_name"
# ^^ this is our issue. Question is.. where the heck is it coming from? It is totally new. but didn't come up when I went through the list initially. at position 176

alcLists[176]
alcLists[175]
alcLists[177]
alcLists[178]
rsIDsInspect <- c('rs4480324', 'rs12673949', 'rs1937522', 'rs62325464')
filtList <- alcLists[rsIDsInspect]
# no idea where "rsik factor" is coming in from ... searched the keys and values of all surrounding entries in the list.. and I am fairly certain its processed in the natural order of the list meaning it is indeed one of these objects I grabbed which is creating the inconsistency.. must trace by debug now.


# getting an issue merging data now because of 'duplicated key values'? thinking there may be repeating rsIDs in my list... and thus reason to filter rsID lists before transformation

rsIDsUnique <- unique(names(alcLists))
?duplicated
dupes <- duplicated(names(alcLists))
sum(dupes) # 480 dupes. ... the question is.. are we getting repetitious data? also more importantly, are we calling Ensemble with duplicated rsIDs?


alcConsumpGWAStable <- alcConsumpRawVars[[1]]
IDs <- alcConsumpGWAStable$VariantID
uniIDs <- unique(IDs)
# without unique applied, we are seeing 1772 vs 1064 ids... more questions: am I filtering for uniqueness before doing an API call? # question 2: Are entire rows duplicated in the GWAS table? if not, why do we see duplicated IDs? (I know it could be coming in from different studies.. and honestly since the data is read straight in from GWAS downloads I don't see any reason to alter its contents at a base level, as its not my doing that there are duplications.. )
# The fact is that calling the Ensemble API for the same rsIDs more than once is pointless, as the same data will be returned each time, so filtering before calling is an absolute necessity. beyond that, I suppose I just have to put in some safe guards against duplicated data in the transform stage... christ.


# testing unique ID calling after moding createMT2 ------------------------

alcConsumpVars <- createMT2(GWASdataSets[1], processData = F)
# ^^ no issues

alcVars <- purrr::flatten(alcConsumpVars[[2]])
uniNamesAlc <- unique(names(alcVars))
d <- duplicated(alcVars)
sum(d) # 21 ... so even though we ensure only unique IDs are called for, some non unique IDs are returned. Meaning they may just be synonyms for each other? (this makes sense assuming the systems prioritizes one name for a variant over another and thus always returns the same name for an entire set of synonyms.) ... though I am assuming this behavior.
# lets try to transform the new alcVars data

alcConsumpVarTransformed <- dbugTransform(alcConsumpVars, popsData = F)
# it worked... presumably removing the majority of duplicates from the variant data by not asking for so much duplicate data the merged table doesn't exceed the set threshold now.
# As expected, the table has grown to some extent... which makes sense in the case of multiple mappings IIRC.. I do




































