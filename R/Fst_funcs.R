# Fst calculating functions:



# -------------------------------------------------------------------------


HudsonFst <- function(n1, n2, p1, p2){
  return( ((p1-p2)**2 - (p1*(1-p1))/(n1-1) - (p2*(1-p2))/(n2-1)) / (p1*(1-p2) + p2*(1-p1)) )
}


# -------------------------------------------------------------------------


# Wright's Fst func -------------------------------------------------------

WrightFst <- function(){

}


# -------------------------------------------------------------------------



# hudsonFst_alleleList() - accepts an allele list of the form returned by GWASpops.pheno2geno::createMT() as "PopAlleleFreqData",
#                          returns a list of DFs which contain calculated hudson Fst as well as the values used to do the calculation in each row. Each row is per population pair.
#                          `deleteRedundants` deletes the data used to calculate Fst, leaving only the population pair indicator and the value itself. This is useful for saving memory, as the values are highly redundant. (I have not come up with a way to store them more efficiently yet, as I lack sufficient insight into the future use of this data)

hudsonFst_alleleList <- function(alleleList, populationsDF, deleteRedundants = FALSE, discardMultiAllelic = TRUE){

  captureList <- list()

  for(i in 1:length(alleleList)){
    tableName <- names(alleleList[i])
    captureList[[tableName]] <- perAlleleFst_transform(alleleList[[i]], populationsDF, deleteRedundants)
  }
  if(discardMultiAllelic){
    # discarding entries which were multiallelic in the capture list
    captureList <- purrr::discard(captureList, is.logical)
  }

  return(captureList)
}


# -------------------------------------------------------------------------



perAlleleFst_transform <- function(alleleDF, populations, deleteRedundants = FALSE){

  if(length(unique(alleleDF$allele)) > 2){ # no calculations for multiallelic sites. This method of Fst calculation isn't suitable to non-biallelic sites.
    return(NA)
  }

  # Extract data of interest from alleleDF
  ancestralAllele <- attr(alleleDF, "Ancestral_Allele")

  if(is.na(ancestralAllele)){ # sometimes ancestral Allele is NA, which makes the rest of this function impossible to execute. Calculating ancestral allele by finding allele with highest frequency.
    ancestralAllele <- calc_ancestralAllele(alleleDF)
  }

  if( !(ancestralAllele %in% unique(alleleDF$allele)) ) { # reassign AA if assignment of AA is somehow wrong, (ancestral allele not found in data.frame)
    ancestralAllele <- calc_ancestralAllele(alleleDF)     #  Wrong assignment can come directly from data sources, not necessarily my own code
  }

  alleleDF <- alleleDF[alleleDF$population %in% populations$Population_Abbreviation & alleleDF$allele != ancestralAllele , ] # filtering down to minor allele.
  # ^^ here we discard non 'ancestral alleles' because we are really looking for the minor alleles to compare fst wrt to.

  # Digest DF in to create DF out:

  DF_rows <- nrow(alleleDF) # number for efficient pairwise iteration

  if(DF_rows == 0){ # when no pops of interest exist for a given variant, we just cancel the function and return nothing
    return(NULL)
  }

  rowHolder <- list()
  for(i in 1:(DF_rows-1)){ # i correlates to a population
    for(j in (1+i):DF_rows){ # j correlates to the second population used to pair with i's population
      rName <- paste0(alleleDF$population[i], "-X-",alleleDF$population[j])
      row <- c(rName,
               populations[populations$Population_Abbreviation == alleleDF$population[i]]$Sample_Count,
               alleleDF$frequency[i],
               populations[populations$Population_Abbreviation == alleleDF$population[j]]$Sample_Count,
               alleleDF$frequency[j])

      rowHolder[[rName]] <- row
    }
  }

  #setup return DF, name cols, fix col types, assign attributes
  retDF <- as.data.frame(t(dplyr::bind_rows(rowHolder)))

  names(retDF) <- c("pop_pair", 'n1', "p1", 'n2', 'p2')
  retDF['n1'] <- as.numeric(retDF[['n1']])
  retDF['n2'] <- as.numeric(retDF[['n2']])
  retDF['p1'] <- as.numeric(retDF[['p1']])
  retDF['p2'] <- as.numeric(retDF[['p2']])

  attr(retDF, "Ancestral_Allele") <- ancestralAllele
  attr(retDF, "VariantID") <- attr(alleleDF, "VariantID")

  fstVec <- numeric(nrow(retDF))
  for(i in 1:nrow(retDF)){
    fstVec[i] <- HudsonFst(retDF[i,2],retDF[i,4],retDF[i,3],retDF[i,5])
  }
  retDF['Fst_Hudson'] <- fstVec

  if(deleteRedundants){ #removing all but population pairs and Fst value to save memory, populations pairs are stored as row names
    retDF <- retDF[,6, drop = FALSE] # drop = FALSE ensures we don't lose the row names in coercion
  }

  return(retDF)
}



# -------------------------------------------------------------------------

calc_ancestralAllele <- function(population_alleleDF){

  alleleSet <- unique(population_alleleDF$allele)

  highestSum <- 0
  AA <- ''
  for(allele_character in alleleSet){
    tempSum <- sum(population_alleleDF[population_alleleDF$allele == allele_character, ]$frequency)
    if(tempSum > highestSum){
      highestSum <- tempSum
      AA <- allele_character
    }
  }
  return(AA)
}




# Old calc_ancestralAllele ------------------------------------------------
# this version is just poorly written. Wasn't thinking clearly about the needs of the func when I wrote it.

#
# calc_ancestralAllele <- function(population_alleleDF){
#
#   A_mag <- sum(population_alleleDF[population_alleleDF$allele == "A", ]$frequency)
#   C_mag <- sum(population_alleleDF[population_alleleDF$allele == "C", ]$frequency)
#   G_mag <- sum(population_alleleDF[population_alleleDF$allele == "G", ]$frequency)
#   T_mag <- sum(population_alleleDF[population_alleleDF$allele == "T", ]$frequency)
#   tempVec <- c(A_mag, C_mag, G_mag, T_mag)
#
#   ancestralAllele <- switch(which.max(tempVec),
#                             '1' = "A",
#                             '2' = "C",
#                             '3' = "G",
#                             '4' = "T")
#   return(ancestralAllele)
# }


# -------------------------------------------------------------------------
