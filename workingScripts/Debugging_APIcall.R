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

save(substanceAbuseVars, file = 'data/debugging_raw_data/substanceAbuseVars.rds')# MARKER... HERE IS WHERE NEW DATA CALLS STOPPED ... THINK I SAVED EVERYTHING, though theres some chance I didn't... SO assuming its all good but possibly these last two may not be.
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

# ALTERNATIVE MEANS OF INSTALLING PACKAGE AFTER UPDATING ITS CONTENTS: Build > More (dropdown) > Clean and Install


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

# Testing new tryCatch() block for initial merge failure...
alcConsumpVarTransformedRAW <- dbugTransform(alcConsumpRawVars, popsData = F)
# SUCCESS!! . Wonderful day. Now merging has a failsafe.
# HOWEVER, the failsafe merge is > 2x the size of the merge where only unique vars was called for.


# how to check for duplicated rows in R again? Checked my notes, and nothing matched... So I am not sure.. I don't immediately need to know which rows are and arent duplicates... I could write a function to do this.. I would be wise to at least search around a bit before doing so though... would be wiser to stay on task and not worry about number of duplicate rows until it actually matters.



# testing transform for population data now: ------------------------------
workHERE11_23 <- function(g){}

library(data.table)

load("data/debugging_raw_data/alcConsumpAlldata.rds")

alcConsumpAllTransformed <- dbugTransform(alcConsumpAlldata, popsData = T)
debug(dbugTransform)
undebug(dbugTransform)

load("data/debugging_raw_data/AirPollutionAlldata.rds")
airPollutionAllTransformed <- dbugTransform(AirPollutionAlldata, popsData = T)
airPollutionTempReturn <- dbugTransform(AirPollutionAlldata, popsData = T)
# OK still getting error with finding merge.data.table... maybe I will just import `merge` as well?
  # AFTER PUTTING DEBUG FUNCTIONS INSIDE OF THE PROJECT: still getting the infinite loop issue I think... Going to once again check for regular createMT efficiacy.. wondering if adding new scripts to the package may have imparted this infinite looping function on creatMT() as well? .. if this executes successfully. I will.. maybe... go in and look at the data.frame that is supposed to be returned by one loops of the lapply() causing issues.. what other leads exist?
  #    # execution was successful for createMT()... So I am going to just let the transform run for 20 mins and see if its still going after that... presumably it is actually stuck in a loop as the total execution of transformation and API calling was less than 20 mins for createMT(). I suppose we can look at what is being returned by the infinite looping return in the case ... ok it just worked... took like 4.5 mins. I think this means we can go back to working on the ancestral allele issue from before... I also think we may want to do some benchmarking once all debugging is said and done.


# BIG ISSUE ^^^: some kind of infinite loop.. or freezing occurring Upon the first(?) return of :
# 'singlePop_alleleFreqDTs <- lapply(Populations$Population_Abbreviation, function(x) GWASpops.pheno2geno:::singlePopTransform(masterList[[2]], targetPopulation = x))'
#   We are just having the function freeze to a halt... I don't know what is going on.... I think a restart and reproducing this bug comes first. maybe I can verify the normal pipelines behavior next by grabbing an install from github and running it for the airPoll data set .. then like debug it to spot any differences? ....

# TESTING BASE FUNCTIONS...
library(data.table)
airPollTest <- createMT('exampleData/air_pollution', population_data = T) # THIS SUCCEEDED WITHOUT DEBUGGER, but GOT STUCK IN THE SAME INIFINITE LOOP AS OUR `dbugTransform()` WHEN IN THE DEBUGGER.

debug(createMT)
# ^^ this works fine. I am pretty confused about the freezing in my dbugTransfrom func. I just don't understand what conditions can lead to some indefinite processing. I don't see any loop repeatedly executing.. its upon hitting a return function that the freeze happens.. which is pretty strange. I feel pretty lost on this. My best lead right now is to compare the DS going into the function in my regular pipline vs the current dev pipeline.. any inconsistencies MAY possible point me to my issue... as the base call is the same at the line causing the freezing.
#

#OK . so without debugging it (createMT) worked fine. WITH debugging on, I got stuck in the same infinite freeze point as before.. at the same return for singlePopTransform. So I am wondering if something about the debugging environment is causing this issue. Something else I noticed, is that we aren't actually seeing any obvious differences in the data going into that function from our debug command compared to the creatMT version... so I don't think its the data in. Gonna try and run the dbugTransform WITHOUT debugger on.. and see if it can transform airpol data again.

# it is looking like we are getting frozen at the same spot again. I don't know if there is any way to gather information about what the heck is happening under the hood to cause this freeze, which is bad, this is the worst sort of error, one where our system doesn't even sort of tell us what is causing it, where its happening, or what the system is doing amidst the error. No transparency error type. My only lead on how to manage this right now is to clean up my working envrionment... basically I want to edit bootCalls.R s.t. the chance for calling a func outside of the packages functions isn't possible... I am sure this will lead to other issues. but maybe cleaning the environment will avoid this weird infinite loop thing.

# OK I PUT DEBUGGING FUNCS IN THE PACKAGE TEMPORARLY S.T. THEY HAVE ACCESS TO NON-EXPORTED FUNCS FOR DEBUGGING AND ENVIRONMENT CONSISTENCY PUROSES


# error message below occurred for this function. Due to time I will now be generating the new, non duplicate filled data while doing other stuff today.

# Error in `filter()`:
#   ! Problem while computing `..1 = x$allele != attr(x, "Ancestral_Allele")`.
# ✖ Input `..1` must be of size 1, not size 2.
# Run `rlang::last_error()` to see where the error occurred.
# Warning messages:
#   1: Unknown or uninitialised column: `population`.
# 2: Unknown or uninitialised column: `population`.
# 3: In x$allele != attr(x, "Ancestral_Allele") :
#   longer object length is not a multiple of shorter object length
# 4: In x$allele != attr(x, "Ancestral_Allele") :
#   longer object length is not a multiple of shorter object length
# 5: In x$allele != attr(x, "Ancestral_Allele") :
#   longer object length is not a multiple of shorter object length



# 11-23-2022 --------------------------------------------------------------
# DEBUGGING POPULATION DATA ISSUES..

# seeking multiple mappings in alc data to understand how to deal with multiple ancestral alleles

alcConsumpEnsData <- purrr::flatten(alcConsumpAlldata[[2]])
mappingNum <- sapply(alcConsumpEnsData, \(x) length(x$mappings))
mappingG2 <- mappingNum[mappingNum > 1]
mappingG2indeces <- which(mappingNum > 1)
mappingG2 # ok so we have a bunch which are actually 3 or more.. and two with 6 and 7 mappings.. woah.
getWeirdMappings <- alcConsumpEnsData[mappingG2indeces[10:15]]
getWeirdMappingsAll <- alcConsumpEnsData[mappingG2indeces]

# inspecting ancestral alleles to get a feel for what I should account for with the assignment script
# SEEING: some alleles only have NULL values.. (is this a problem? Should we be filling with some alternate default value?)

onlyMappings <- lapply(getWeirdMappingsAll, \(x) x$mappings)
onlyAA <- lapply(onlyMappings, \(x) for(ele in 1:length(x)){ })

AAlist <- as.list(vector(length = 53))
num = 0
for(i in onlyMappings){ # would be nice to not need a whole for loop to access these nested elements and to get them packed into their own arrays for each element of the base list 'onlyMappings'
  num = num + 1
  subList <- list()
  for(ele in 1:length(i)){
    subList[ele] <- purrr::pluck(onlyMappings[[num]][[ele]], "ancestral_allele")
  }
  AAlist[[num]] <- subList
 }

warnings()

a <- pluck(onlyMappings[[1]][[1]], 'ancestral_allele')
?pluck

fAAlist <- flatten(AAlist)

#Not seeing any rsID with more than 1 ancestral allele, though there still exists the issue of NULL values. I think it would be wise to allow for assignment of multiple ancestral alleles... maybe by smashing them into a string with some divider character, and in the cases where there is only NULL AAs, we would instead assign some default value. Finally if there is already an AA and NULL values are encountered, we should discard the NULL value. Also I guess we need to deal with if we find a non-null value and the default value has already been assigned somehow.


# why am I discarding ancestral alleles from population data once again? ... I think it is because its not interesting to look at ancestral alleles in graphs? Yeah... when graphing a set of variants for a single population... if we graph both the ancestral alleles as well as the alleles which are being reported by a given variantID (as variantIDs describe abnormalities in the genome, and thus a variants ancestral allele is the nucleiotide(s) being replaced.) then the graph doesn't really make sense right away.. as we should hypothetically have 100% of frequency for both ancestral allele and variant allele. Thus we just graph the variant frequency for each allele in a given population to see how frequent the occurence of that variant is in that population.
#
#

# OK thinking I will just ensure there is only 1 ancestral allele per rsID... thus the line using it for comparison will execute succesfully.
# If we encounter an issue where there is more than one reported ancestral allele in some data set.. maybe we think of a solution. Really not sure what would be best, there should only be one reported for any given variant.. we could just discard the row entirely? .... or try to find some way to suggest the actual AA. Or just choose at random and discard at random.. (not a valid solution imo)



workHERE11_23_NUM2 <- function(g){}

library(data.table)

load("data/debugging_raw_data/alcConsumpAlldata.rds")

alcConsumpAllTransformed <- dbugTransform(alcConsumpAlldata, popsData = T)
debug(dbugTransform)
undebug(dbugTransform)

load("data/debugging_raw_data/AirPollutionAlldata.rds")
airPollutionAllTransformed <- dbugTransform(AirPollutionAlldata, popsData = T)
airPollutionTempReturn <- dbugTransform(AirPollutionAlldata, popsData = T)



# we have identified where NAs are coming from to an extent in attr_ancestralAllele... will be attempting to fix
# # THIS IS THE CODE TO REMOVE THE NA ROW FROM THE DF... fixMT <- uniqMasterT[!is.na(uniqMasterT$EnsVar_name), ]
# Fixing the AncAllele problem allowed the code to execute successfully.. it was just FUCKING slow. So we will likely need to work on some optimization down the line.

# post fix stated above ^^
a <- alcConsumpAlldata[[2]][[1]]
b <- alcConsumpAlldata[[2]][[2]]
alcConsumpPopDataReduced <- list(a,b) # grabbing a single sublist.. reducing the amount of Ensembl data by 10x
alcConsumpPopDataReduced <- list(alcConsumpAlldata[[1]],alcConsumpPopDataReduced)
alcConsumpAllTransformed2 <- dbugTransform(alcConsumpPopDataReduced, popsData = T) # success.




# RUNNING TRANSFORM FUNCS FOR DEBUGGING ON ALL DATA -----------------------
# STRATEGY. Need to first ensure that all Vars calls are working in the transform function. Then I need to do the same for populations. This is essentially linear debugging where I will handle issue case by case... Not sure what to do about data inconsistencies which may arise outside of this tested data set.

debug(dbugTransform)
undebug(dbugTransform)

options(error = browser)
options(error = recover)
options(error = NULL)

source('bootCalls.R')

load("data/debugging_raw_data/alcConsumpVars.rds")
alcConsumpVarTrans <- dbugTransform(alcConsumpVars)
load("data/debugging_raw_data/alcConsumpAlldata.rds")
alcConsumpAllTrans <- dbugTransform(alcConsumpAlldata, popsData = T) # error here... once again, Ancestral Allele issue... weird. Wonder if I didn't save my progress last time...  (12-15-2022) Looks like I didn't manage to recomplie the package which is why the executed version of the dbugTransform() func doesn't contain the code I would expect it to. I really need to define startup and windown processes that are a necessary part of my workflow whenever I am working on code.

load("data/debugging_raw_data/bCarcinomaVars.rds")
bCarcinomaVarTrans <- dbugTransform(bCarcinomaVars) # ERROR HERE: `Error in `*tmp*`[[1]] : subscript out of bounds` ... I believe this error is coming from a `failed` name:key pair being introduced in this call. Only 2 instances of it across ~1.3k lists. Thus its certainly uncommon. Going to discard this data, as it seems to represent variants which lack mappings onto reference genomes being called against within Ensembl's database, and thus the data is too incomplete to be considered.
# FIXED ^^ The solution is to remove the entries in CONT_noMultiMapping after MultiMappings are removed (flattened out)... using a mask which is built by scanning each sublist in the superlist for names that match failed.
# mask1 <- sapply(CONT_noMultiMapping, \(x) "failed" names(x))
# badentries <- CONT_noMultiMapping[mask1]
# goodEntries <- CONT_noMultiMapping[!mask1]
# Orginal solution was sufficient. Deleting code of second solution (which was almost identical)


load("data/debugging_raw_data/bCarcinomaAlldata.rds")
alcConsumpAllTrans <- dbugTransform(bCarcinomaAlldata, popsData = T)

load("data/debugging_raw_data/colorectalCancerVars.rds")
colorectalCancerVarTrans <- dbugTransform(colorectalCancerVars)
load("data/debugging_raw_data/colorectalCancerAlldata.rds")
colorectalCancerAllTrans <- dbugTransform(colorectalCancerAlldata, popsData = T)

load("data/debugging_raw_data/IBFVars.rds")
IBFVarTrans <- dbugTransform(IBFVars)
load("data/debugging_raw_data/IBFAlldata.rds")
IBFAllTrans <- dbugTransform(IBFAlldata, popsData = T)

load("data/debugging_raw_data/IntVars.rds")
IntVarTrans <- dbugTransform(IntVars)
load("data/debugging_raw_data/IntAlldata.rds")
IntAllTrans <- dbugTransform(IntAlldata, popsData = T)

load("data/debugging_raw_data/lungCancerVars.rds")
lungCancerVarTrans <- dbugTransform(lungCancerVars)
load("data/debugging_raw_data/lungCancerAlldata.rds")
lungCancerAllTrans <- dbugTransform(lungCancerAlldata, popsData = T)

load("data/debugging_raw_data/malabsorptionSyndVars.rds")
malabsorptionSyndVarTrans <- dbugTransform(malabsorptionSyndVars)
load("data/debugging_raw_data/malabsorptionSyndAlldata.rds")
malabsorptionSyndAllTrans <- dbugTransform(malabsorptionSyndAlldata, popsData = T)

load("data/debugging_raw_data/neuroticismVars.rds")
neuroticismVarTrans <- dbugTransform(neuroticismVars)
load("data/debugging_raw_data/neuroticismAlldata.rds")
neuroticismAllTrans <- dbugTransform(neuroticismAlldata, popsData = T)

load("data/debugging_raw_data/prostateCancerVars.rds")
prostateCancerVarTrans <- dbugTransform(prostateCancerVars) # FINAL ERROR FOR Var data -- ERROR HERE Error in dbugTransform(prostateCancerVars) :object 'masterTable' not found ...
# FIXED ^^ . The issue was in my understanding of the behavior of tryCatch() blocks, specifically the somewhat strange activity of the `error = ...` portion.
# The problem of these duplicate rich tables still exists however and is not trivial from my base judgement. I would need to determine if any of the merged tables posses any true duplicate rows and if they do I would need to learn where they were introduced and prevent that from occurring. Ideally my transformation should never fall back into the Cartesian expansion for the merging. So I guess I can write a function which is used to test tables for duplicate rows first... and then I can start debugging further. I really do think solving this issue is worth the time it will take, as otherwise downstream issues would be introduced.
load("data/debugging_raw_data/prostateCancerAlldata.rds")
prostateCancerAllTrans <- dbugTransform(prostateCancerAlldata, popsData = T)

load("data/debugging_raw_data/substanceAbuseVars.rds")
substanceAbuseVarTrans <- dbugTransform(substanceAbuseVars)
load("data/debugging_raw_data/substanceAbuseAlldata.rds")
substanceAbuseAllTrans <- dbugTransform(substanceAbuseAlldata, popsData = T)

load("data/debugging_raw_data/AirPollutionVars.rds")
AirPollutionVarTrans <- dbugTransform(AirPollutionVars)
load("data/debugging_raw_data/AirPollutionAlldata.rds")
AirPollutionAllTrans <- dbugTransform(AirPollutionAlldata, popsData = T)


# MAJOR PROGRESS. ALL VARIANT TRANSFORMATIONS NOW WORK. TIME TO MOVE ONTO WORKING ON POPULATION TRANSFORMATION ISSUES.

# -------------------------------------------------------------------------


# Testing scoping issue hypothesis for trycatch issue ---------------------

dog <- list(LETTERS)

tryCatch(
  expr = {cat <-LETTERS},
)
# ^^ assignment works fine here
#  cat <- <- fog  ..... this is a generic reproducible error


cat <- LETTERS*2 # this produces: `Error in LETTERS * 2 : non-numeric argument to binary operator`

catss <- c(1,2,3,4,6) # functions on its own no problem

rm(cat, catss)
tryCatch(
  expr = {cat <- LETTERS*2},
  error = function(e) {catss <- c(1,2,3,4,6); message('error has executed')}
)
# ^^ despite the error clearly executing I am not seeing any cats assignment occuring in my environment. I honestly thought we encountered this before and it was somehow fixed.


tryCatch(
  expr = {cat <- LETTERS*2},
  error = function(e) {catss <- c(1,2,3,4,6); message('error has executed')},
  finally = {return(catss)}
)
# ^^ : `Error in tryCatch(expr = { : object 'catss' not found`


objectt <- tryCatch(
  expr = {cat <- LETTERS*2},
  error = function(e) {catss <- c(1,2,3,4,6); message('error has executed');return(catss)}
)
# ^^ this works. But honestly it is super weird. what about when cat does execute? what does object become?

rm(objectt)
objectt <- tryCatch(
  expr = {cat <- LETTERS},
  error = function(e) {catss <- c(1,2,3,4,6); message('error has executed');return(catss)}
)

# yep. objectt became the value of whatever is naturally returned by the tryCatch block. So it seems like its built into these try catch blocks that ... something is returned by them... and so they're a sort of 'handles one object only' sort of code block. Which should work for my purposes.

# so do I need the return statement to return the object created in the code of the error block?
rm(objectt, cat, dog)
objectt <- tryCatch(
  expr = {cat <- LETTERS*2},
  error = function(e) {catss <- c(1,2,3,4,6); message('error has executed')}
)
# ^^ objectt was returned as NULL ... so YES we do need a return statement inside of the error block. I really need to learn more about error catching. sheesh.



# 12-19-2022 WORKSPACE ----------------------------------------------------
 # previous workspace is now notes dense, therefore I will be copying necessary bits over here and working off of it

debug(dbugTransform)
undebug(dbugTransform)

options(error = browser)
options(error = recover)
options(error = NULL)

source('bootCalls.R')
library(tictoc)


load("data/debugging_raw_data/alcConsumpVars.rds")
alcConsumpVarTrans <- dbugTransform(alcConsumpVars)
load("data/debugging_raw_data/alcConsumpAlldata.rds")
alcConsumpAllTrans <- dbugTransform(alcConsumpAlldata, popsData = T) # error here... once again, Ancestral Allele issue... weird.
#^^ FIXED .. runtime is over 15 mins though.

load("data/debugging_raw_data/bCarcinomaVars.rds")
bCarcinomaVarTrans <- dbugTransform(bCarcinomaVars)
load("data/debugging_raw_data/bCarcinomaAlldata.rds")
bCarcinomaAllTrans <- dbugTransform(bCarcinomaAlldata, popsData = T) # bug here. ... for some reason output was still produced here. NAMING ISSUE... resolved
# FIXED ^^^ Same fix as lungCancer worked here.
# Error in `filter()` at GWASpops.pheno2geno/R/pipeline_helperFuncs.R:119:4:
#   ! Problem while computing `..1 = x$allele != attr(x, "Ancestral_Allele")`.
# ✖ Input `..1` must be of size 1, not size 0.


load("data/debugging_raw_data/colorectalCancerVars.rds")
colorectalCancerVarTrans <- dbugTransform(colorectalCancerVars)
load("data/debugging_raw_data/colorectalCancerAlldata.rds")
colorectalCancerAllTrans <- dbugTransform(colorectalCancerAlldata, popsData = T) # works

load("data/debugging_raw_data/IBFVars.rds")
IBFVarTrans <- dbugTransform(IBFVars)
load("data/debugging_raw_data/IBFAlldata.rds")
IBFAllTrans <- dbugTransform(IBFAlldata, popsData = T) # works 50+ warnings again

load("data/debugging_raw_data/IntVars.rds")
IntVarTrans <- dbugTransform(IntVars)
load("data/debugging_raw_data/IntAlldata.rds")
IntAllTrans <- dbugTransform(IntAlldata, popsData = T) # works .. large data.. extremely slow. No warnings

load("data/debugging_raw_data/lungCancerVars.rds")
lungCancerVarTrans <- dbugTransform(lungCancerVars)
load("data/debugging_raw_data/lungCancerAlldata.rds")
lungCancerAllTrans <- dbugTransform(lungCancerAlldata, popsData = T) # ERROR
# Error in `filter()` at GWASpops.pheno2geno/R/pipeline_helperFuncs.R:119:4:
#   ! Problem while computing `..1 = x$allele != attr(x, "Ancestral_Allele")`.
# ✖ Input `..1` must be of size 1, not size 0.
# FIXED ^^^ added a line to AncestralAllele_Attr() which converts empty entries for the attribute to NA entries.

load("data/debugging_raw_data/malabsorptionSyndVars.rds")
malabsorptionSyndVarTrans <- dbugTransform(malabsorptionSyndVars)
load("data/debugging_raw_data/malabsorptionSyndAlldata.rds")
malabsorptionSyndAllTrans <- dbugTransform(malabsorptionSyndAlldata, popsData = T) # works

load("data/debugging_raw_data/neuroticismVars.rds")
neuroticismVarTrans <- dbugTransform(neuroticismVars)
load("data/debugging_raw_data/neuroticismAlldata.rds")
neuroticismAllTrans <- dbugTransform(neuroticismAlldata, popsData = T) # works

load("data/debugging_raw_data/prostateCancerVars.rds")
prostateCancerVarTrans <- dbugTransform(prostateCancerVars)
load("data/debugging_raw_data/prostateCancerAlldata.rds")
prostateCancerAllTrans <- dbugTransform(prostateCancerAlldata, popsData = T) # works

load("data/debugging_raw_data/substanceAbuseVars.rds")
substanceAbuseVarTrans <- dbugTransform(substanceAbuseVars)
load("data/debugging_raw_data/substanceAbuseAlldata.rds")
substanceAbuseAllTrans <- dbugTransform(substanceAbuseAlldata, popsData = T) # works

load("data/debugging_raw_data/AirPollutionVars.rds")
AirPollutionVarTrans <- dbugTransform(AirPollutionVars)
load("data/debugging_raw_data/AirPollutionAlldata.rds")
AirPollutionAllTrans <- dbugTransform(AirPollutionAlldata, popsData = T) # works





#ERRRORS:
# Error in `filter()` at GWASpops.pheno2geno/R/pipeline_helperFuncs.R:119:4:
# ! Problem while computing `..1 = x$allele != attr(x, "Ancestral_Allele")`.
# ✖ Input `..1` must be of size 1, not size 2.
# Run `rlang::last_error()` to see where the error occurred.
# Warning messages:
#   1: Unknown or uninitialised column: `population`.
# 2: Unknown or uninitialised column: `population`.
# 3: In x$allele != attr(x, "Ancestral_Allele") :
#   longer object length is not a multiple of shorter object length
# 4: In x$allele != attr(x, "Ancestral_Allele") :
#   longer object length is not a multiple of shorter object length
# 5: In x$allele != attr(x, "Ancestral_Allele") :
#   longer object length is not a multiple of shorter object length





# Optimization work -------------------------------------------------------


debug(dbugTransformOPT)
undebug(dbugTransformOPT)

source('bootCalls.R')
library(tictoc)


# VARIANT DATA TRANSFORMS

load("data/debugging_raw_data/malabsorptionSyndVars.rds")
malabsorptionSyndVarTrans <- dbugTransformOPT(malabsorptionSyndVars)

load("data/debugging_raw_data/AirPollutionVars.rds")
AirPollutionVarTrans <- dbugTransformOPT(AirPollutionVars)

load("data/debugging_raw_data/prostateCancerVars.rds")
prostateCancerVarTrans <- dbugTransformOPT(prostateCancerVars)

load("data/debugging_raw_data/colorectalCancerVars.rds")
colorectalCancerVarTrans <- dbugTransformOPT(colorectalCancerVars)

load("data/debugging_raw_data/substanceAbuseVars.rds")
substanceAbuseVarTrans <- dbugTransformOPT(substanceAbuseVars)

load("data/debugging_raw_data/lungCancerVars.rds")
lungCancerVarTrans <- dbugTransformOPT(lungCancerVars)

load("data/debugging_raw_data/bCarcinomaVars.rds")
bCarcinomaVarTrans <- dbugTransformOPT(bCarcinomaVars)

load("data/debugging_raw_data/IBFVars.rds")
IBFVarTrans <- dbugTransformOPT(IBFVars)

load("data/debugging_raw_data/alcConsumpVars.rds")
alcConsumpVarTrans <- dbugTransformOPT(alcConsumpVars)

load("data/debugging_raw_data/neuroticismVars.rds")
neuroticismVarTrans <- dbugTransformOPT(neuroticismVars)

load("data/debugging_raw_data/IntVars.rds")
IntVarTrans <- dbugTransformOPT(IntVars)




# POPULATION DATA TRANSFORMS

load("data/debugging_raw_data/malabsorptionSyndAlldata.rds")
malabsorptionSyndAllTrans <- dbugTransformOPT(malabsorptionSyndAlldata, popsData = T)

load("data/debugging_raw_data/AirPollutionAlldata.rds")
AirPollutionAllTrans <- dbugTransformOPT(AirPollutionAlldata, popsData = T)

load("data/debugging_raw_data/prostateCancerAlldata.rds")
prostateCancerAllTrans <- dbugTransformOPT(prostateCancerAlldata, popsData = T)

load("data/debugging_raw_data/colorectalCancerAlldata.rds")
colorectalCancerAllTrans <- dbugTransformOPT(colorectalCancerAlldata, popsData = T)

load("data/debugging_raw_data/substanceAbuseAlldata.rds")
substanceAbuseAllTrans <- dbugTransformOPT(substanceAbuseAlldata, popsData = T)

load("data/debugging_raw_data/lungCancerAlldata.rds")
lungCancerAllTrans <- dbugTransformOPT(lungCancerAlldata, popsData = T)

load("data/debugging_raw_data/bCarcinomaAlldata.rds")
bCarcinomaAllTrans <- dbugTransformOPT(bCarcinomaAlldata, popsData = T) # 50+ warnings ... wonder if its connected to the name mismatches in GWASc table and popFreqLists

load("data/debugging_raw_data/IBFAlldata.rds")
IBFAllTrans <- dbugTransformOPT(IBFAlldata, popsData = T) # 50+ warnings

load("data/debugging_raw_data/alcConsumpAlldata.rds")
alcConsumpAllTrans <- dbugTransformOPT(alcConsumpAlldata, popsData = T)# 50+ warnings

load("data/debugging_raw_data/neuroticismAlldata.rds")
neuroticismAllTrans <- dbugTransformOPT(neuroticismAlldata, popsData = T)

load("data/debugging_raw_data/IntAlldata.rds")
IntAllTrans <- dbugTransformOPT(IntAlldata, popsData = T)




# order of data to test:
# 1. malabs
# 2. airPol
# 3. prostateCancer
# 4. Colorectal
# 5. SubstanceAbus
# 6. lungCancer
# 7. breastCarcinoma
# 8. IBF
# 9. alcConsump
# 10. neuroticism
# 11. Int




# 12-22-2022
#  Question to start today is to try and get rid of warnings first... or to just start reintegrating code? ... reintegrating code is more fundamentally important... If I thought the warnings were important to maintaining data integrity I would prioritize that. But I don't have a way of testing for data integrity to be honest, and it feels a bit beyond the current scope of my ability / this project in my mind. So I may just go straight into integration, and then I will think about the warnings at a later point maybe.

# First test all functions against data for success.. Then begin reintegration... ensure any functions called (like createMT2() ) have their changes migrated over first... then change the names in debug functions to that... then reintegrate properly.? ... may take more planning.





# Integration work:  ------------------------------------------------------


debug(dbugTransformOPT)
undebug(dbugTransformOPT)

source('bootCalls.R')


# VARIANT DATA TRANSFORMS

load("data/debugging_raw_data/malabsorptionSyndVars.rds")
malabsorptionSyndVarTrans <- ensListTransform(malabsorptionSyndVars)

load("data/debugging_raw_data/AirPollutionVars.rds")
AirPollutionVarTrans <- ensListTransform(AirPollutionVars)

load("data/debugging_raw_data/prostateCancerVars.rds")
prostateCancerVarTrans <- ensListTransform(prostateCancerVars)

load("data/debugging_raw_data/colorectalCancerVars.rds")
colorectalCancerVarTrans <- ensListTransform(colorectalCancerVars)

load("data/debugging_raw_data/substanceAbuseVars.rds")
substanceAbuseVarTrans <- ensListTransform(substanceAbuseVars)

load("data/debugging_raw_data/lungCancerVars.rds")
lungCancerVarTrans <- ensListTransform(lungCancerVars)

load("data/debugging_raw_data/bCarcinomaVars.rds")
bCarcinomaVarTrans <- ensListTransform(bCarcinomaVars)

load("data/debugging_raw_data/IBFVars.rds")
IBFVarTrans <- ensListTransform(IBFVars)

load("data/debugging_raw_data/alcConsumpVars.rds")
alcConsumpVarTrans <- ensListTransform(alcConsumpVars)

load("data/debugging_raw_data/neuroticismVars.rds")
neuroticismVarTrans <- ensListTransform(neuroticismVars)

load("data/debugging_raw_data/IntVars.rds")
IntVarTrans <- ensListTransform(IntVars)


# POPULATION DATA TRANSFORMS

load("data/debugging_raw_data/malabsorptionSyndAlldata.rds")
malabsorptionSyndAllTrans <- ensListTransform(malabsorptionSyndAlldata, popsData = T)

load("data/debugging_raw_data/AirPollutionAlldata.rds")
AirPollutionAllTrans <- ensListTransform(AirPollutionAlldata, popsData = T)

load("data/debugging_raw_data/prostateCancerAlldata.rds")
prostateCancerAllTrans <- ensListTransform(prostateCancerAlldata, popsData = T)

load("data/debugging_raw_data/colorectalCancerAlldata.rds")
colorectalCancerAllTrans <- ensListTransform(colorectalCancerAlldata, popsData = T)

load("data/debugging_raw_data/substanceAbuseAlldata.rds")
substanceAbuseAllTrans <- ensListTransform(substanceAbuseAlldata, popsData = T)

load("data/debugging_raw_data/lungCancerAlldata.rds")
lungCancerAllTrans <- ensListTransform(lungCancerAlldata, popsData = T)

load("data/debugging_raw_data/bCarcinomaAlldata.rds")
bCarcinomaAllTrans <- ensListTransform(bCarcinomaAlldata, popsData = T)

load("data/debugging_raw_data/IBFAlldata.rds")
IBFAllTrans <- ensListTransform(IBFAlldata, popsData = T)

load("data/debugging_raw_data/alcConsumpAlldata.rds")
alcConsumpAllTrans <- ensListTransform(alcConsumpAlldata, popsData = T)

load("data/debugging_raw_data/neuroticismAlldata.rds")
neuroticismAllTrans <- ensListTransform(neuroticismAlldata, popsData = T)

load("data/debugging_raw_data/IntAlldata.rds")
IntAllTrans <- ensListTransform(IntAlldata, popsData = T)









# post integration testing ------------------------------------------------

# I have attempted to integrate the changes by modifying the 3 key functions createMT(), get_ensVariants(), and ensListTransform() ... also multiAPIcall_variants() was modified as well. ... before all transformation was taking place inside of the master func as well as the API get calling func... now its happening in its own space entirely which is much better organizationally. The question is however, did I put all the pieces together correctly, does the path execute properly. Going to take some time to test, as I will be hitting the API as a part of the whole process now.

# ERROR: Can't recycle `rs2161719` (size 10) to match `rs6900057` (size 11). FIXED###########

debug(createMTfinal)

AirPollutionVars <- createMTfinal('exampleData/air_pollution')
# Error in rsO_list[[i]][["mappings"]] : subscript out of bounds .... looks like the data is being over flattened before processing in fixMultiMapping... FIXED
#
AirPollutionAlldata <- createMTfinal('exampleData/air_pollution', population_data = T)

alcConsumpVars <- createMTfinal(GWASdataSets[1])
alcConsumpAlldata <- createMTfinal(GWASdataSets[1], population_data = T)

bCarcinomaVars <- createMTfinal(GWASdataSets[2])
bCarcinomaAlldata <- createMTfinal(GWASdataSets[2], population_data = T)

colorectalCancerVars <- createMTfinal(GWASdataSets[3])
colorectalCancerAlldata <- createMTfinal(GWASdataSets[3], population_data = T)

IBFVars <- createMTfinal(GWASdataSets[4])
IBFAlldata <- createMTfinal(GWASdataSets[4], population_data = T)

IntVars <- createMTfinal(GWASdataSets[5])
IntAlldata <- createMTfinal(GWASdataSets[5], population_data = T)

lungCancerVars <- createMTfinal(GWASdataSets[6])
lungCancerAlldata <- createMTfinal(GWASdataSets[6], population_data = T)

malabsorptionSyndVars <- createMTfinal(GWASdataSets[7])
malabsorptionSyndAlldata <- createMTfinal(GWASdataSets[7], population_data = T) #  ERROR
# Error in rsO_list[[i]] : subscript out of bounds
# 3.
# GWASpops.pheno2geno:::fixMultiMapping(CONT) at Debugging_APIcall_functions.R#553
# 2.
# ensListTransform(allData, popsData = T) at Debugging_APIcall_functions.R#676
# 1.
# createMTfinal(GWASdataSets[7], population_data = T)

neuroticismVars <- createMTfinal(GWASdataSets[8])
neuroticismAlldata <- createMTfinal(GWASdataSets[8], population_data = T)

prostateCancerVars <- createMTfinal(GWASdataSets[9])
prostateCancerAlldata <- createMTfinal(GWASdataSets[9], population_data = T)

substanceAbuseVars <- createMTfinal(GWASdataSets[10])
substanceAbuseAlldata <- createMTfinal(GWASdataSets[10], population_data = T)





# 1-7-23 Saving data for graphing -----------------------------------------

save(AirPollutionVars, file = 'data/transformed_data_for_graphing/airPollution.rds')
save(AirPollutionAlldata, file = 'data/transformed_data_for_graphing/airPollutionPops.rds')

save(alcConsumpVars , file = 'data/transformed_data_for_graphing/alcConsump.rds')
save(alcConsumpAlldata  , file = 'data/transformed_data_for_graphing/alcConsumpPops.rds')

save(bCarcinomaVars , file = 'data/transformed_data_for_graphing/bCarcinoma.rds')
save(bCarcinomaAlldata  , file = 'data/transformed_data_for_graphing/bCarcinomaPops.rds') # an error arose while this was executing.. saving anyway...

save(colorectalCancerVars , file = 'data/transformed_data_for_graphing/colorectalCancer.rds')
save(colorectalCancerAlldata  , file = 'data/transformed_data_for_graphing/colorectalCancerPops.rds')

save(IBFVars , file = 'data/transformed_data_for_graphing/IBF.rds')
save(IBFAlldata  , file = 'data/transformed_data_for_graphing/IBFPops.rds')

save(IntVars , file = 'data/transformed_data_for_graphing/Int.rds') # START BACK HERE
save(IntAlldata  , file = 'data/transformed_data_for_graphing/IntPops.rds')

save(lungCancerVars , file = 'data/transformed_data_for_graphing/lungCancer.rds')
save(lungCancerAlldata  , file = 'data/transformed_data_for_graphing/lungCancerPops.rds')

save(malabsorptionSyndVars , file = 'data/transformed_data_for_graphing/malabsorptionSynd.rds')
save(malabsorptionSyndAlldata  , file = 'data/transformed_data_for_graphing/malabsorptionSyndPops.rds')

save(neuroticismVars , file = 'data/transformed_data_for_graphing/neuroticism.rds')
save(neuroticismAlldata  , file = 'data/transformed_data_for_graphing/neuroticismPops.rds')

save(prostateCancerVars , file = 'data/transformed_data_for_graphing/prostateCancer.rds')
save(prostateCancerAlldata  , file = 'data/transformed_data_for_graphing/prostateCancerPops.rds')

save(substanceAbuseVars , file = 'data/transformed_data_for_graphing/substanceAbuse.rds')
save(substanceAbuseAlldata, file = 'data/transformed_data_for_graphing/substanceAbusePops.rds')





# NOVEL ERROR from API call? on "bCarcinomaVars <- createMTfinal(GWASdataSets[2])"
# Error in curl::curl_fetch_memory(url, handle = handle) :
# Operation was aborted by an application callback
# 9.
# curl::curl_fetch_memory(url, handle = handle)
# 8.
# request_fetch.write_memory(req$output, req$url, handle)
# 7.
# request_fetch(req$output, req$url, handle)
# 6.
# request_perform(req, hu$handle$handle)
# 5.
# POST(baseURL, content_type("application/json"), accept("application/json"),
#      body = rsID_Array) at Debugging_APIcall_functions.R#743
# 4.
# get_ensVariantsFinal(splitList[[i]]) at Debugging_APIcall_functions.R#841
# 3.
# GWASpops.pheno2geno:::multiAPIcall_variantsFinal(rsIDs) at Debugging_APIcall_functions.R#730
# 2.
# get_ensVariantsFinal(uniqueVariantIDs) at Debugging_APIcall_functions.R#685
# 1.
# createMTfinal(GWASdataSets[2])

# IMPORTANT NOTE ^^^^^ this error didn't replicate upon a second run of the call... therefore there really is potential for the underlying functions I have built upon to be points of failure on an inconsistent basis, and that should be documented for users

# 2nd novel error from API calling:
# Error in get(helpTopicsName, envir = .rs.toolsEnv()) :
# object '.completions.helpTopics' not found
# In addition: Warning message:
#   Expected 2 pieces. Additional pieces discarded in 14 rows [65, 120, 121, 122, 123, 125, 126, 127, 128, 129, 130, 131, 132, 1945].





# NOVEL UNEXPECTED WARNINGS:  ---------------------------------------------

# CALL: colorectalCancerVars <- createMTfinal(GWASdataSets[3])
# Warning message:
# Expected 2 pieces. Additional pieces discarded in 18 rows [203, 675, 676, 690, 691, 692, 693, 694, 695, 696, 925, 926, 927, 928, 929, 930, 954, 955].


# CALL: colorectalCancerAlldata <- createMTfinal(GWASdataSets[3], population_data = T)
# Warning messages:
# 1: Expected 2 pieces. Additional pieces discarded in 18 rows [203, 675, 676, 690, 691, 692, 693, 694, 695, 696, 925, 926, 927, 928, 929, 930, 954, 955].
# 2: In doTryCatch(return(expr), name, parentenv, handler) :
#   restarting interrupted promise evaluation









# FINAL INTEGRATION -------------------------------------------------------

# need to check that integration is successful. Ideally we would once again run through every data set, however I will settle for just a hand full.. as the very exact same code SHOULD exist between the current '*final()' funcs and thus we just need to put things where they go, rebuild the package and verify a few quick test cases.





AirPollutionVars <- createMT('exampleData/air_pollution')
AirPollutionAlldata <- createMT('exampleData/air_pollution', population_data = T)

alcConsumpVars <- createMT(GWASdataSets[1])
alcConsumpAlldata <- createMT(GWASdataSets[1], population_data = T)

bCarcinomaVars <- createMT(GWASdataSets[2])
bCarcinomaAlldata <- createMT(GWASdataSets[2], population_data = T)


