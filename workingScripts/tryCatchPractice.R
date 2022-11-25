# Practicing with trycatch()


vec1 <- c("5", "5", "5", "fice", "6")
vec2 <- c(0, 0, 0, 0, 0)

for(i in 1:length(vec1)){
  vec2[i] <- bad.numeric(vec1[i])
}

bad.numeric <- function(value){
  v <- as.numeric(value)
  if(is.na(v)){
    stop("Ya funked up brotha, that shit can't be made into a number!")
  }
  return(v) #.. having no return results in a loop stopping error, which is what I need. For some reason having bad.numeric throw an error doesn't stop the execution though?
}



# working on tryCatch now -------------------------------------------------

# for any value that throws an error with bad.numeric() we should see 20 in that position and the loop should complete.
tryCatchTest <- function(vector1, vector2){

  for(i in 1:length(vector1)){
    tryCatch(
      expr = {vector2[i] <- bad.numeric(vector1[i])},

      error = function(e){ vector2[i] <- 20; message('the ERROR condition executed!') },
      warning = function(W){ vector2[i] <- 20; message('the warning condition executed!')}
            )
  }
}

v <- tryCatchTest(vec1, vec2)

vec3 <- c("5", "5", "5", "7", "6")
v2 <- tryCatchTest(vec1, vec3) #verified that the tryCatch block isn't breaking basic funcitonality.

debug(tryCatchTest)
v <- tryCatchTest(vec1, vec2)

vec4 <- as.numeric(vec1) # its vectorized... but why cant I do a single value in my function? after one call of as.numeric on a single value in my vector, I am getting a fully processed output.. not value by value

val1 <- as.numeric(vec1[1L]) #this is working for outputting a single value ... not sure why the funny business is happening in the for loop...

tryCatchTest2 <- function(vector1){

  vector2 <- list()

  for(i in 1:length(vector1)){
    tryCatch(
      expr = {vector2[i] <- bad.numeric(vector1[i])},

      error = function(e){ vector2[i] <- 20; message('the ERROR condition executed!') },
      warning = function(W){ vector2[i] <- 20; message('the warning condition executed!')},
      finally = { vector2[i] <- 20; cat("the finally block executed!!\n") }
    )
  }
  return(vector2)
}


vec5 <- tryCatchTest2(vec1) # assignment doesn't appear to be possible within the error / warning blocks. It is possible within the finally block though... Even in the debugger I never see the list updated for the value that should be added when an error or warning occurs. I am a bit confused about how to properly use the tryCatch() if I can dynamically respond to the error with more code.

debug(tryCatchTest2)
undebug(tryCatchTest2)
?as.vector

tryCatchTest3 <- function(vector1){

  vector2 <- list()

  for(i in 1:length(vector1)){
    tryCatch(
      expr = {vector2[i] <- bad.numeric(vector1[i])},

      error = function(e){  message('the ERROR condition executed!'); return(20) },# hoping the return will assign 20 to the list
      warning = function(W){ message('the warning condition executed!'); return(20) }
    )
  }
  return(vector2)
} # getting closer I think to understanding how to use this tool... I obviously need more practice, as I cannot get it to do what I want... which is simple assignment of value to a for loop specified position in a list of some default value upon catching errors or warnings.

vec6 <- tryCatchTest3(vec1)




# 11-21-2022 TESTING APPLIED USE FOR MERGE --------------------------------

?data.table::merge
# no faith in this working.. I just don't think you can execute code normally from the error portion of a tryCatch block... I don't really know how to use them generally.
masterMerge <- function(list){

  GWAS_DF <- list[[1]]
  CONT_Table <- list[[2]]

  tryCatch(
    expr = {masterTable <- merge.data.table(GWAS_DF, CONT_Table, by.x = 'VariantID', by.y = 'EnsVar_name')},
    error = function(e){
              masterTable <- merge.data.table(GWAS_DF, CONT_Table, by.x = 'VariantID', by.y = 'EnsVar_name', allow.cartesian = T)}

  )
  return(masterTable)
}


testMM <- masterMerge(alcConsumpVarTransformed) # so this worked. Wasn't expecting success...


# First lets make sure the call we want to execute works at all.
library(data.table)
masterTable <- data.table::merge(alcConsumpVarTransformed[[1]], alcConsumpVarTransformed[[2]], by.x = 'VariantID', by.y = 'EnsVar_name', allow.cartesian = T) # getting very strange error "Error: 'merge' is not an exported object from 'namespace:data.table'"

masterTable <- merge(alcConsumpVarTransformed[[1]], alcConsumpVarTransformed[[2]], by.x = 'VariantID', by.y = 'EnsVar_name', allow.cartesian = T) # works SUCCESS?

masterTableTEST <- merge(alcConsumpVarTransformed[[1]], alcConsumpVarTransformed[[2]], by.x = 'VariantID', by.y = 'EnsVar_name')
# ^^ getting the same error as before... ... meaning we are def using the data.table version... the fact we cannot explicitely call it by its full name is super weird though.

methods(merge)

class(alcConsumpVarTransformed[[1]]) # data.table data.frame
class(alcConsumpVarTransformed[[2]]) # tbl_df tbl data.frame
# OK so I am assuming that we are seeing automatic calling of the data.table::merge() method instead of base::merge() because the first argument being passed to the call is in fact a data.table and there are some class methods occurring under the hood here which mediate preferences in function calls??? (no idea how classes are funcitoning for R under the hood.. because like... I don't ever really intentionally use them.)


masterTable2 <- merge.data.table(alcConsumpVarTransformed[[1]], alcConsumpVarTransformed[[2]], by.x = 'VariantID', by.y = 'EnsVar_name', allow.cartesian = T) # SUCCESS
masterTable3 <-  merge.data.frame(alcConsumpVarTransformed[[1]], alcConsumpVarTransformed[[2]], by.x = 'VariantID', by.y = 'EnsVar_name', allow.cartesian = T) # ? ALSO WORKS? wtf.. shouldnt the allow.cartesian be a data.table only thing? Maybe they both are routing to the same function based entirely off of class? .... seems the result is identical.



data.table:::merge.data.table # this will print out the code behind this function.
base:::merge

getAnywhere("merge.data.table") # similarly we are printing out the code... but with some meta-data too?

?getAnywhere

# FROM DETAILS OF MAN PAGE FOR MERGE: merge is a generic function in base R. It dispatches to either the merge.data.frame method or merge.data.table method depending on the class of its first argument. Note that, unlike SQL, NA is matched against NA (and NaN against NaN) while merging.
# ^^ this suggests my assumption about class based calling is indeed correct.
















