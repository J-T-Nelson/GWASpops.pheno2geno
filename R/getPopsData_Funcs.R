# get pops Data Funcs
# -------------------------------------------------------------------------

getPopsData <- function(rsIDChunkList, nChunks, startingChunk, reportNumErrors = TRUE){

  retList <- list()

  for (i in 1:nChunks){
    chunk <- (startingChunk + i - 1)
    chnkName <- paste0("Chunk_", chunk, "_CONT")

    retList[[chnkName]] <- tryCatch(
      expr = {
        GWASpops.pheno2geno:::get_ensVariants(rsIDChunkList[[chunk]], population_data = TRUE)
      },

      error = function(e){
        warning(paste0("Error occured for ", chnkName))
        return(chnkName) # this should just return a character vector instead of a list which will be our means of identifying error counts
      }
    )
  }

  fileName <- paste0("chunk", startingChunk, "-", (startingChunk+nChunks-1) , ".rds")

  setwd("D:\\Programming\\R_projects\\Kulathinal_Lab\\GWASpops.pheno2geno\\workingData\\unprocessedChunks")
  save(retList, file = fileName)
  setwd("../")

  if(reportNumErrors){
    numErrors <- sum(sapply(retList, is.character))
    message(paste0("Number of empty chunks returned: ", numErrors))
    message("\nEmpty chunks are returned as character vectors when an error is caught, or when the API fails to return expected data.\n")
  }

  return(retList)
}


# -------------------------------------------------------------------------

# grabChunks() is a wrapper for getPopsData which will continuously call for chunks until hitting the num calls limit, or until it has reached the spcified stopping point which comes when all chunks are retrieved.

grabChunks <- function(data, StartChunk = 1,chunksPerCall = 100, numCalls = 240){
  #numCalls at default of 240 should mean that by default this would just call for all of the data

  allrsID_ch10 <- data

  for(i in 1:numCalls){
    startPoint <- (StartChunk + (i - 1)*chunksPerCall) # incrementally updates starting chunk relative to starting point

    if(startPoint > 23300){
      message("All chunks should be grabbed except the last few")
      return()
    }

    getPopsData(allrsID_ch10, nChunks = chunksPerCall, startPoint, reportNumErrors = FALSE)
  }

  message("call finished without automatic termination; i.e. not all chunks grabbed yet, but specified amount should be saved.")
  return()
}

# -------------------------------------------------------------------------




