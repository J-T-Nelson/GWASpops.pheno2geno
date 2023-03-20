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


# updated 3-17-23: deals with ancestral alleles of frequency = 1.0 properly now by adding rows for the minor allele with 0.0 frequency.
#                   - also rounds all negative Fst to 0,
#                   - as well as register fst as 0 where NaN is generated.
#                   (NaN generated when two pops have 0.0 for minor allele freq)
#
perAlleleFst_transform <- function(alleleDF, populations, deleteRedundants = FALSE){

  if(length(unique(alleleDF$allele)) > 2){ # no calculations for multiallelic sites. This method of Fst calculation isn't suitable to non-biallelic sites.
    return(NA)
  }

  # Extract data of interest from alleleDF
  ancestralAllele <- attr(alleleDF, "Ancestral_Allele")

  if(is.na(ancestralAllele)){ # sometimes ancestral Allele is NA, which makes the rest of this function impossible to execute. Calculating ancestral allele by finding allele with highest frequency.
    ancestralAllele <- calc_ancestralAllele(alleleDF)
  } else {

    if( !(ancestralAllele %in% unique(alleleDF$allele)) ) { # reassign AA if assignment of AA is somehow wrong, (ancestral allele not found in data.frame)
      ancestralAllele <- calc_ancestralAllele(alleleDF)     #  Wrong assignment can come directly from data sources, not necessarily my own code
    }
  }

  # adding rows where they are missing when an allele's ancestral allele is fixed at 1.0
  if(any( (alleleDF$frequency == 1) & (alleleDF$allele == ancestralAllele) )){
    nonAncestralAllele <- unique(alleleDF$allele[ alleleDF$allele != ancestralAllele ] )
    if(length(nonAncestralAllele) == 0){ # some data doesn't actually report on the non-ancestral allele... rs1555226898 specifically. Due to complexity of adding in new rows, (and difficulties predicting downstream affects) I am removing such examples
      return(NA)
    }
    mutateRows <- alleleDF[ alleleDF$frequency == 1 , ]
    mutateRows$frequency <- 0
    mutateRows$allele_count <- 0
    mutateRows$allele <- nonAncestralAllele
    alleleDF <- rbind.data.frame(alleleDF, mutateRows)
  }

  alleleDF <- alleleDF[alleleDF$population %in% populations$Population_Abbreviation & alleleDF$allele != ancestralAllele , ] # filtering down to minor allele.
  # ^^  discard non 'ancestral alleles' because we are looking for the minor alleles to compare fst wrt to.

  # Digest DF in to create DF out:

  DF_rows <- nrow(alleleDF) # number for efficient pairwise iteration

  if(DF_rows < 2){ # when no pops of interest exist for a given variant or only 1 row remains, cancel the function and return nothing
    return(NA)
  }

  rowHolder <- list()
  for(i in 1:(DF_rows-1)){ # i correlates to a population
    for(j in (1+i):DF_rows){ # j correlates to the second population used to pair with i's population
      rName <- paste0(alleleDF$population[i], "-X-",alleleDF$population[j])

      row <- c(populations[populations$Population_Abbreviation == alleleDF$population[i]]$Sample_Count,
               alleleDF$frequency[i],
               populations[populations$Population_Abbreviation == alleleDF$population[j]]$Sample_Count,
               alleleDF$frequency[j])

      rowHolder[[rName]] <- row
    }
  }

  #setup return DF, name cols, fix col types, assign attributes
  retDF <- as.data.frame(t(dplyr::bind_rows(rowHolder)))

  names(retDF) <- c('n1', "p1", 'n2', 'p2')
  retDF['n1'] <- as.numeric(retDF[['n1']])
  retDF['n2'] <- as.numeric(retDF[['n2']])
  retDF['p1'] <- as.numeric(retDF[['p1']])
  retDF['p2'] <- as.numeric(retDF[['p2']])

  attr(retDF, "Ancestral_Allele") <- ancestralAllele
  attr(retDF, "VariantID") <- attr(alleleDF, "VariantID")

  fstVec <- numeric(nrow(retDF))
  for(i in 1:nrow(retDF)){
    fstVec[i] <- HudsonFst(retDF[i,1],retDF[i,3],retDF[i,2],retDF[i,4])
  }

  fstVec[is.nan(fstVec)] <- 0 # for cases where p1 and p2 are 0 we get NaN out of 'HudsonFst()'
  fstVec[fstVec < 0] <- 0 # bringing all negatives up to 0 for accurate averaging and realistic Fst values
  retDF['Fst_Hudson'] <- fstVec

  if(deleteRedundants){ #removing all but population pairs and Fst value to save memory, populations pairs are stored as row names
    retDF <- retDF[,5, drop = FALSE] # drop = FALSE ensures we don't lose the row names in coercion
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
