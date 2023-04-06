# Monogenic diseases associated SNPs table generation, top SNPs by popsVSmetapop table generation, top SNPs by 1% of summed FST for each pop table generation, intersection of tables 2 and 3 generation.


# Monogenic Table generation ----------------------------------------------


load("./workingData/subTraitsList.rds") # not useful. object is called 'testSave'


# original function used to do regex grab of monogenics
#
#
filterByTrait <- function(assoTable, regexList = "none", nameList) {

  # ensuring data.frame class for bind_rows later
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
        retList[nameList[i]] <- dplyr::bind_rows(retList[nameList[i]], DF, DF2) # this adds rows

      }

    } else{ # when only 1 pattern, hit table for that pattern against 2 cols of interest
      DF <- assoTable[grep(regexList[[i]][1], assoTable$`DISEASE/TRAIT`, ignore.case = T), ]
      DF2 <- assoTable[grep(regexList[[i]][1], assoTable$MAPPED_TRAIT, ignore.case = T), ]
      DF$DATE <- as.character(DF$DATE)
      DF2$DATE <- as.character(DF2$DATE)
      DF$`DATE ADDED TO CATALOG` <- as.character(DF$`DATE ADDED TO CATALOG`)
      DF2$`DATE ADDED TO CATALOG` <- as.character(DF2$`DATE ADDED TO CATALOG`) # character conversions are necessary for data integrity
      retList[[nameList[i]]] <- dplyr::bind_rows(retList[nameList[i]], DF, DF2)
    }
  }

  retList <- lapply(retList, base::unique) # remove duplicate rows from all DFs in list
  return(retList)

}


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


# testing filterByTrait() :


original_monogenics <- filterByTrait(asso, regexList, DnameList) # looks like the original DS we made. Saving
save(original_monogenics, file="./workingData/monogenic_dataStructures/original_monogenics_list.rds")
load("./workingData/monogenic_dataStructures/original_monogenics_list.rds")
novel_monogenics_1 <- data.table::fread("./workingData/monogenic_disease_genes_table/Stylianos_2001__PGD.csv", header = TRUE)[,1:3]
compendium_monogenics <- data.table::fread("./workingData/monogenic_disease_genes_table/Compendium_of_causative_genes_common_monogenic_disorders.csv")
rm(novel_monogenics_2)

# ok now we need to differentially process each of these three data sources for the same desired output.
# The objective is to grab rows from the SNPanno_GWAS_ensembl table for each of our 3 sources.
#   1. Each one has its own process, we will thus make 3 sub-table,
#   2. then will bind them together,
#   3. then filter down to the unique rows
#   4. then merge in Fst data
#   5. then save the novel data structure and proceed with other tasks.


# originial_monogenics:

original_monogenics <- original_monogenics[sapply(original_monogenics, nrow) > 0]
length(original_monogenics) # 8
str(original_monogenics[[1]]) # tibble

# going to use a 2 feature key to grab rows of interest. First flatten the list, then filter out unneeded features, then use the key to grab data from our table of interest SNPanno_GWAS_ensembl

# adding in necessary meta-data for later combination into final monogenic table
ogMono_addDiseaseNames <- lapply(1:length(original_monogenics), function(i){
    original_monogenics[[i]][,"monogenic_disease_name"] <- rep(names(original_monogenics[i]), nrow(original_monogenics[[i]])) })

original_monogenics <- lapply(1:length(original_monogenics), function(i){cbind(original_monogenics[[i]], ogMono_addDiseaseNames[i])})
rename <- c(colnames(original_monogenics[[1]])[1:39], "monogenic_disease_name")
rename
for(i in 1:length(original_monogenics)){ colnames(original_monogenics[[i]]) <- rename}
colnames(original_monogenics[[1]]) # success .. toughest rename of my life... like 20 minutes

origMonoFlat <- dplyr::bind_cols(purrr::flatten(original_monogenics))
origMonoFlat <- dplyr::bind_rows(original_monogenics)

head(origMonoFlat) # not allowed to view origMonoFLat for some reason...
class(origMonoFlat)
colnames(origMonoFlat)
origMonoFlatKey <- origMonoFlat[,c("SNPS", "DISEASE/TRAIT", "monogenic_disease_name")] # looks good.
nrow(origMonoFlatKey) # 502
origMonoFlatKey[, "monogenic_searched_resource"] <-
  rep("A. Nesterova, Monogenic rare diseases in biomedical databases and text mining", nrow(origMonoFlatKey))

mask_1 <- SNPanno_GWAS_ensembl$VariantID %in% origMonoFlatKey$SNPS & SNPanno_GWAS_ensembl$`DISEASE/TRAIT` %in% origMonoFlatKey$`DISEASE/TRAIT`
sum(mask_1) # 456 .. less than the number of rows, which means some of the candidates may have been lost in data retrieval, most are still here though
456/502 # 0.9083665
monoSubTable_1 <- SNPanno_GWAS_ensembl[mask_1,] # good
dim(monoSubTable_1) # 456  57 ... need to add on disease names and source

monoSubTable_1[, "monogenic_searched_resource"] <-
  rep("A. Nesterova, Monogenic rare diseases in biomedical databases and text mining", nrow(monoSubTable_1))
monoSubTable_1[, "monogenic_disease_name"] <-
  rep("SEE DISEASE/TRAIT", nrow(monoSubTable_1))

# monoSubTable_1 <- dplyr::left_join(monoSubTable_1, origMonoFlatKey[2:4], by = c("DISEASE/TRAIT"="DISEASE/TRAIT"))
# dim(monoSubTable_1) # 24249    59 ... unexcusable amount of rows added. Why.
# nrow(unique(monoSubTable_1))

# Error in forderv(x, by = by, sort = FALSE, retGrp = TRUE) :
#   Column 44 passed to [f]order is type 'list', not yet supported.
#
# Due to error above, not sure how to check for duplicated values here.

# going to leave the disease associated empty.. no way to easily do this and the info isn't so important

## novel_monogenics_1:

# since we don't have terms for our diseases of interest in term of the 'TRAITS/DISEASES' column terminology we will have to search purely on the genes within the novel_monogenics_1 table.
# This may grab some unrelated values, thus, perhaps it would be wise to merge in the associated disease terms as well as the paper source for them in order to improve data completeness and interpret-ability

length(novel_monogenics_1$V2) #96
length(unique(novel_monogenics_1$V2)) #94
mask_2 <- SNPanno_GWAS_ensembl$MAPPED_GENE %in% novel_monogenics_1$V2
sum(mask_2) # 1474
monoSubTable_2 <- SNPanno_GWAS_ensembl[mask_2,]

# merge in disease names from source as well as searched_resource

colnames(novel_monogenics_1) <- c("monogenic_disease_name", "Gene", "monogenic_searched_resource")
dim(monoSubTable_2)# 1474 57

monoSubTable_2 <- merge(monoSubTable_2, novel_monogenics_1, by.x = "MAPPED_GENE", by.y = "Gene") # more rows than wanted .. typical error generated.
monoSubTable_2 <- dplyr::left_join(monoSubTable_2, novel_monogenics_1, by = c("MAPPED_GENE"="Gene"))
dim(monoSubTable_2) # 1588   59 ... still wrong.

?left_join
nrow(novel_monogenics_1) #96
nrow(unique(novel_monogenics_1)) #96
length(unique(novel_monogenics_1$monogenic_disease_name)) # 94
length(unique(novel_monogenics_1$Gene)) # 94
# So we have two repeated values in two cols, which dont create repeating rows, which means we necessarily have some novel rows introduced, meaning some Genes are associated with more than one disease term and vise-versa

# I think the added rows must presumably be generted because of the unique rows which posses duplicate terms by column. Thus we will go forward beleiving our data is good.

# We need to remember to add the same data cols to our other monoSubTables as well to complete the data. Then from there we can filter it all down and save the product



## Compendium_monogenics:
##
colnames(compendium_monogenics)
compendiumFiltered <- compendium_monogenics[,c("% cases caused by majority gene", "Majority Gene (alone)")]
compendiumFiltered[,"monogenic_searched_resource"] <- rep("Tucker L. Apgar, Charles Sanders, Compendium of causative genes and their encoded proteinsfor common monogenic disorders", nrow(compendium_monogenics))

compendiumFiltered[,"monogenic_disease_name"] <- compendium_monogenics$`Disease or group of diseases`

# removing any genes from this list where '% caused by majority gene" is below 75% ... Don't want genes poorly associted with disesase in the analysis
compendiumFiltered <- compendiumFiltered[compendiumFiltered$`% cases caused by majority gene` > .74,]

mask_3 <- SNPanno_GWAS_ensembl$MAPPED_GENE %in% compendiumFiltered$`Majority Gene (alone)`
sum(mask_3) # 1754

monoSubTable_3 <- SNPanno_GWAS_ensembl[mask_3,]
nrow(monoSubTable_3) # 1754

#joining, but excluding '% causeed by majority gene' col for later binding of rows
monoSubTable_3 <- dplyr::left_join(monoSubTable_3, compendiumFiltered[,2:4], by = c("MAPPED_GENE"="Majority Gene (alone)"))
nrow(monoSubTable_3) # 1906  ... once again, presumably we are seeing duplications due to duplicate genes or disease names in compendiumFiltered, thus requiring more than 1 row for complete mapping.
dim(monoSubTable_3) # 1906   59


# Combine subtables:
#

monogenic_variants_complete <- rbind(monoSubTable_1, monoSubTable_2, monoSubTable_3)
nrow(unique(monogenic_variants_complete))

## COMPLETED Creation of monogenic table


# -------------------------------------------------------------------------











































