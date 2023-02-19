# drafting the filtering by trait
#

colnames(Associations)
assoTraits <- Associations[,`DISEASE/TRAIT`]


filt <- assoTraits %in% 'Breast cancer'
filt <- assoTraits %in% 'Breast cancer'
BC <- assoTraits[filt]

assoTraits <- tolower(assoTraits)

head(assoTraits)

A <- Associations
testGrep <- A[grep("Breast", A$`DISEASE/TRAIT`), ] # this is a great tool for partial string matching

?grep
returnsss <- grep("Breast", A$`DISEASE/TRAIT`)
returnsss # a bunch of indices

cfTest <- grep("[Cc]ystic", A$`DISEASE/TRAIT`) # 203 things returned
cfTest <- grep("[Cc]ystic [Ff]ib", A$`DISEASE/TRAIT`) # 87
CFrows <- A[cfTest,]

cysticFib <- assoTraits[assoTraits %in% 'cystic fibrosis'] # empty

EL <- vector('list')
class(EL)
EL['dog'] <- 'homjom'
EL['dogs'] <- 5
EL
EL['dog'] <- 'flipflap' # reassignment
EL

regexList[5]
length(regexList[[1]]) #1
length(regexList[[5]]) #2
regexList[[5]][1]
regexList[[5]][2]

df <- data.frame()
A$MAPPED_TRAIT


# verify bind_rows function -----------------------------------------------

CFrows2 <- CFrows
cfrowdouble <- bind_rows(CFrows, CFrows2) # correct number of rows are here 174

#Check for recursive behavior

cfrowdouble <- bind_rows(cfrowdouble, CFrows2) # correct num, 261
87*2 #174
87*3 # 261

# remove dupe rows:  ------------------------------------------------------

ur <- unique(cfrowdouble) # without anything else this just removed the unique rows.. back to 87 rows
?unique # ok I am def using the base R version... there is a data.table version.


# Trait Look Up Func ------------------------------------------------------


filterByTrait <- function(assoTable, regexList = "none", nameList) {

   # ensuring data.table class for bind_rows later
  assoTable <- as.data.frame(assoTable)


  retList <- vector('list')

  for(i in 1:length(regexList)){ #go through all lookups of interest, specified by contents of regexList
    retList[[nameList[i]]] <- data.frame() # make list-entry to fill with data looked up

    if(length(regexList[[i]])>1){ #if list contains more than one pattern for a trait to look up, hit all patterns
      for(j in 1:length(regexList[[i]])){
        if(j==2){ # 2nd entries in lists should be abbreviations of targeted conditions
          DF <- assoTable[grep(regexList[[i]][j], assoTable$`DISEASE/TRAIT`), ]
          DF2 <- assoTable[grep(regexList[[i]][j], assoTable$MAPPED_TRAIT), ]
        } else{
            DF <- assoTable[grep(regexList[[i]][j], assoTable$`DISEASE/TRAIT`, ignore.case = T), ]
            DF2 <- assoTable[grep(regexList[[i]][j], assoTable$MAPPED_TRAIT, ignore.case = T), ]
        }
        DF$DATE <- as.character(DF$DATE)
        DF2$DATE <- as.character(DF2$DATE)
        DF$`DATE ADDED TO CATALOG` <- as.character(DF$`DATE ADDED TO CATALOG`)
        DF2$`DATE ADDED TO CATALOG` <- as.character(DF2$`DATE ADDED TO CATALOG`) # character conversions are necessary for data integrity
        retList[nameList[i]] <- bind_rows(retList[nameList[i]], DF, DF2) # this adds rows

      }

    } else{ #when only 1 pattern, hit table for that pattern against 2 cols of interest
      DF <- assoTable[grep(regexList[[i]][1], assoTable$`DISEASE/TRAIT`, ignore.case = T), ]
      DF2 <- assoTable[grep(regexList[[i]][1], assoTable$MAPPED_TRAIT, ignore.case = T), ]
      DF$DATE <- as.character(DF$DATE)
      DF2$DATE <- as.character(DF2$DATE)
      DF$`DATE ADDED TO CATALOG` <- as.character(DF$`DATE ADDED TO CATALOG`)
      DF2$`DATE ADDED TO CATALOG` <- as.character(DF2$`DATE ADDED TO CATALOG`) # character conversions are necessary for data integrity
      retList[[nameList[i]]] <- bind_rows(retList[nameList[i]], DF, DF2)
    }
  }

  retList <- lapply(retList, base::unique) # remove duplicate rows from all DFs in list
  return(retList)

}


# filterByTrait() Test ----------------------------------------------------
library(tidyverse)

test1 <- filterByTrait(A, regexList = regexList,DnameList)
t1UNI <- lapply(test1, base::unique)# errors... werid coercion of dates was causing issues. fixed
t11 <- test1[[1]]
t11UNI <- unique(t11)
debug(filterByTrait)



# ERROR: ran for 1 sec before throwing this.. need to investigate with debugger FIXED
#
# Error in retList[nameList[i]] <- data.frame() :
# replacement has length zero
# In addition: Warning message:
#   In retList[nameList[i]] <- bind_rows(retList[nameList[i]], DF, DF2) :
#   number of items to replace is not a multiple of replacement length


# Reporting processing.  --------------------------------------------------
nrow
numHits <- sapply(t1UNI, nrow)
numHits
names(numHits)
numHitsN <- unname(numHits)
numHitsN
?sort
sort(numHits, decreasing = T)

results <- t1UNI
results <- results[numHitsN > 0] # filter out the DFs w/ no rows
?split
tSplit <- split(results[[1]], results[[1]]$MAPPED_TRAIT)
TStraits <- sapply(tSplit, function(x) unique(x$MAPPED_TRAIT))
unname(TStraits)

names(results[1])

captureList <- list()
for(i in 1:length(results)){
  tempS <- split(results[[i]], results[[i]]$MAPPED_TRAIT)
  subTraits <- sapply(tempS, function(x) unique(x$MAPPED_TRAIT))
  subTraits <- unname(subTraits)
  captureList[[ names(results[i]) ]] <- subTraits
}
captureList

tibCap <- tibble(captureList)
DFcap <- as.data.frame(captureList)
?as.data.frame

testSave <- captureList

save(testSave, file = "subTraitsList.rds")
rm(testSave)
load("subTraitsList.rds") # looks like I can just load in saved objects to Rmd files..

# grep all relevant DISEASE/TRAITS for all monogenics
# get a list DFs containing results for each set of monogenic trait searches
# split each DF in list of DFs by diseases/traits ... manually discard irrelevant ones
# merge all relevant data back together
# report on num associations and num traits for each monogenic disease.
# Assemble a little report in some way that makes sense... Rmd? What would a data analyst do? Maybe ppt? Maybe watch a vid or two on the subject to make sense of it

# list looks good right now.
regexList <- list("cystic.?fibrosis",
               "hemophilia.?A",
               "retinitis.?pigmentosa",
               "Hemophilia.?B",
               c("H?e?r?e?d?i?t?a?r?y?.?breast.??a?n?d?.??ovarian.?cancer.?syndrome",
               "HBOC"),
               "Huntington.?disease",
               "Neurofibromatosis.?type.?1?",
               "Duchenne.?muscular.?dystrophy",
               "sickle.?cell.?a?n?e?m?i?a?",
               c("Tuberous.?sclerosis.?c?o?m?p?l?e?x?",
               "TSC"),
               c("Charcot.?Marie.?Tooth.?disease",
               "CMT"),
               "Marfan.?syndrome",
               "Fanconi.?anemia",
               "Titin.?related.?l?i?m?b?.?g?i?r?d?l?e?.?muscular.?dystrophy.?R?1?0?",
               "ataxia.?telangiectasia",
               c("Severe.?combined.?immunodeficiency",
               "SCID"),
               c("Lynch.?syndrome",
               "HNPCC",
               "h?e?r?e?d?i?t?a?r?y?.?nonpolyposis.?colorectal.?cancer"),
               "Congenital.?hypothyroidism",
               c("Familial.?adenomatous.?polyposis",
               "FAP"),
               c("Phenylketonuria",
               "PKU"),
               c("Rett.?syndrome",
               "RTS",
               "cerebroatrophic.?hyperammonemia"),
               c("Microphthalmia.?anophthalmia.?coloboma",
               "MAC"),
               "Joubert.?syndrome",
               c("Primary.?ciliary.?dyskinesia",
               "PCD"),
               c("Autosomal.?dominant.?polycystic.?kidney.?disease",
               "ADPKD"),
               c("Wilson.?disease",
               "[Hh]epatolenticular.?[Dd]egeneration"),
               c("Limb.?girdle.?muscular.?dystrophy",
               "LGMD"),
               c("Arryhthmogenic.?right.?entricular.?cardiomyopathy",
               "ARVC"),
               c("Chronic.?granulomatous.?d?i?s?e?a?s?e?",
               "CGD"),
               c("Hereditary.?spastic.?paraplegia",
               "HSP"),
               c("Familial.?thoracic.?aortic.?aneurysm.?and.?dissection",
               "f?a?m?i?l?i?a?l?.?TAAD"),
               c("E?a?r?l?y?.?infantile.?epileptic.?encephalopathy",
               "EIEE",
               "ohtahara.?syndrome"),
               "Alport.?Syndrome")

DnameList <- c("Cystic Fibrosis ",
              "Hemophilia A",
              "Retinitis pigmentosa",
              "Hemophilia B",
              "Hereditary breast and ovarian cancer syndrome",
              "Huntington disease",
              "Neurofibromatosis type 1",
              "Duchenne muscular dystrophy",
              "Sickle cell anemia",
              "Tuberous sclerosis complex",
              "Charcot-Marie-Tooth disease / Hereditary motor",
              "Marfan syndrome",
              "Fanconi anemia",
              "Titin-related limb-girdle muscular dystrophy R10",
              "ataxia-telangiectasia",
              "Severe combined immunodeficiency",
              "Lynch syndrome",
              "Congenital hypothyroidism",
              "Familial adenomatous polyposis",
              "Phenylketonuria",
              "Rett syndrome",
              "Microphthalmia-anophthalmia-coloboma",
              "Joubert syndrome",
              "Primary ciliary dyskinesia",
              "Autosomal dominant polycystic kidney disease",
              "Wilson disease",
              "Limb-girdle muscular dystrophy",
              "Arryhthmogenic right entricular cardiomyopathy",
              "Chronic granulomatous disease",
              "Hereditary spastic paraplegia",
              "Familial thoracic aortic aneurysm and",
              "Early infantile epileptic encephalopathy",
              "Alport Syndrome")
