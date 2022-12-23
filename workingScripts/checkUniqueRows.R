# Going to be testing for a good way to obtain only unique or only duplicated rows of a data frame

class(gwasData)
?unique

testDF <- gwasData[1:20,]
dupetestDF <- rbind(testDF, testDF[1:5,])
# should be 5 duplicated rows

dupRows <- duplicated(dupetestDF)
sum(dupRows) # 5
dupRows

# lets now test a larger DF

dupeRowsGWAS <- duplicated(gwasData)
sum(dupeRowsGWAS) # 0

dupeRowsPC <- duplicated(prostateCancerVarTrans) # ERROR: Error in forderv(x, by = query$by, sort = FALSE, retGrp = TRUE) :
                                                    #Column 26 passed to [f]order is type 'list', not yet supported.

class(prostateCancerVarTrans)
class(prostateCancerVarTrans[1,]) # not testing each column... still treating it as a whole DF

class(prostateCancerVarTrans$VariantID)

for(i in 1:40) {

}

?col

DFversion <- as.data.frame(prostateCancerVarTrans)
class(DFversion) #"data.frame"
class(DFversion[1,]) #"data.frame"
class(DFversion[,1:10]) #"data.frame"



classesOfDFlist <- sapply(DFversion, class) # this is good
classesOfDF <- as.character(sapply(DFversion, class)) # this is good
classesOfDF

# time to quickly replace the list col and try testing uniqueness again

EnsVar_synonyms <- as.character(DFversion$EnsVar_synonyms)

#can I actually do it in place?

DFversion$EnsVar_synonyms <- as.character(DFversion$EnsVar_synonyms)

class2 <- sapply(DFversion, class)
class2
DFversion$ # <- using this notation is useful for identifying classes too, as at least lists have a different symbol associated with them

duplicated(DFversion)
dupeRowsPC <- duplicated(DFversion)
sum(dupeRowsPC) # 1560
# so this is very ugly ... we should just use this trick to remove them in our scripts.

PCnoDupes <- prostateCancerVarTrans[!dupeRowsPC ] # this works.




# checking if this works for data.table

prostateCancerVarTrans$EnsVar_synonyms <- as.character(prostateCancerVarTrans$EnsVar_synonyms)
prostateCancerVarTrans$ # success no lists.



























































































































