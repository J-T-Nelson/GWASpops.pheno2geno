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


#### still need to add Fst Values to this table, possibly filter out cols which aren't meaningful for its purpose.

colnames(monogenic_variants_complete)
length(unique(monogenic_variants_complete$VariantID)) # 1111
allUniqueMonogenicVars <- unique(monogenic_variants_complete$VariantID)


monogenic_variants_complete <- monogenic_variants_complete[,c("VariantID", "DISEASE/TRAIT", "MAPPED_GENE", 'CONTEXT', 'INTERGENIC', 'RISK ALLELE FREQUENCY',
                                                              'MAPPED_TRAIT', 'EnsVar_most_severe_consequence', 'Parent term', 'monogenic_searched_resource',
                                                              'monogenic_disease_name')]

save(monogenic_variants_complete, file = "./workingData/monogenic_disease_genes_table/monogenic_variants_complete.rds")

## COMPLETED Creation of monogenic table


# -------4-9-23 session start------------------------------------------------------------------
library(ggplot2)
library(ggrepel)

# Metapop vs subpops table:

diseaseFst_popsVSall # starting data
save(diseaseFst_popsVSall, file = "./workingData/diseaseFst_popsVSall.rds")

# count number of FST above .4
sum(diseaseFst_popsVSall >= .4) # 1207
sum(diseaseFst_popsVSall[,1] >= .4) # 6 ... Checking to ensure we are looking at all above .4 across cols.. we are.

# summary statistics and boxplot:

summary(diseaseFst_popsVSall)
highFst_count_perPop <- data.frame(apply(diseaseFst_popsVSall, MARGIN = 2, \(x){ sum(x >= .4)}))
highFst_count_perPop <- data.frame(highFst_count_perPop, superPops) # now the super populations are displayed as well
colnames(highFst_count_perPop) <- c("Count", "Sub_Population", "Super_Population")



bar_plot <- ggplot(highFst_count_perPop, aes(x = rownames(highFst_count_perPop), y = Count)) +
  geom_bar(stat = 'identity',color = "black", fill = "dodgerblue") +
  labs(title = "Fst Disease Variants >= 0.4",
       x = "Sub-population",
       y = "Count") +
  theme_minimal(base_size = 14)

print(bar_plot)

bar_plot_2 <- ggplot(highFst_count_perPop[order(highFst_count_perPop$Super_Population),],
                     aes(x = rownames(highFst_count_perPop),y = Count, fill = Super_Population)) +
  geom_bar(stat = 'identity',color = "black") +
  labs(title = "Fst Disease Variants >= 0.4",
       x = "Sub-population",
       y = "Count") +
  theme_minimal(base_size = 14)

print(bar_plot_2)

getwd()
setwd("./GeneralPlots/presentation_plots/report11_12/")
ggsave(filename = "barPlot_diseaseAbove4.png", scale = 3, device = 'png', bg = 'white') # white argument worked to fix the background, scale 3 works for this plot

# ?ggsave -- WANTED TO KNOW WHAT I COULD USE FOR THE BACKGROUND COLOR... DOCUMENTATION IS SPARSE
# Here are some examples of values you can use for the bg argument:
#
#   Color name: You can use any standard color name recognized by R, such as "white", "black", "red", "blue", "green", etc.
# Hexadecimal color code: You can use a hexadecimal color code, such as "#FF5733", "#2E86C1", or "#28B463".
# RGB color: You can use the rgb() function to create a color by specifying the red, green, and blue components, e.g., rgb(0.5, 0.2, 0.8).



plotFstAboveThreshold_bar <- function(fstData, threshold, output_fileName){

  highFst_count_perPop <- data.frame(apply(fstData, MARGIN = 2, \(x){ sum(x >= threshold)}))
  highFst_count_perPop <- data.frame(highFst_count_perPop, superPops)
  colnames(highFst_count_perPop) <- c("Count", "Sub_Population", "Super_Population")

  titleN <- paste0("Fst Disease Variants >= ", threshold)

  bar_plot_2 <- ggplot(highFst_count_perPop[order(highFst_count_perPop$Super_Population),],
                       aes(x = rownames(highFst_count_perPop),y = Count, fill = Super_Population)) +
    geom_bar(stat = 'identity',color = "black") +
    labs(title = titleN,
         x = "Sub-population",
         y = "Count") +
    theme_minimal(base_size = 14)

  fileN <- paste0(output_fileName, ".png")
  ggsave(filename = fileN, scale = 3, device = 'png', bg = 'white')
  #print(bar_plot_2)
}

thresholds <- c(.60, .55, .50, .45, .35, .3, .25, .2)
for(val in thresholds){
  fName <- gsub('\\.', '_', paste0('barPlot_diseaseAbove', val))
  plotFstAboveThreshold_bar(diseaseFst_popsVSall, val, fName)

}



# create monogenic Fst all and Fst disease table.
#
#  Fst all table:

allUniqueMonogenicVars[1:10]

monogenicFst_all <- full_fst_fixed[rownames(full_fst_fixed) %in% allUniqueMonogenicVars,]
monogenicFst_all <- monogenicFst_all[, grep('ALL', colnames(monogenicFst_all))]
dim(monogenicFst_all) # 835  31 ... 835 of 1111 names retrieved.

monogenicFst_disease <- diseaseFst_popsVSall[rownames(diseaseFst_popsVSall) %in% allUniqueMonogenicVars, ]
dim(monogenicFst_disease) # 283  31 .... substantially less within the disease subset. Meaning we may have a pretty ineffective subset of traits.. or a lot of traits associated with variants associated with monogenic terms / genes, don't fall into association with disease alone, but instead fall into association with non-disease traits as well.

# Lets get some numeber



plotFstAboveThreshold_bar_monogenic_all <- function(fstData, threshold, output_fileName){

  highFst_count_perPop <- data.frame(apply(fstData, MARGIN = 2, \(x){ sum(x >= threshold)}))
  highFst_count_perPop <- data.frame(highFst_count_perPop, superPops)
  colnames(highFst_count_perPop) <- c("Count", "Sub_Population", "Super_Population")

  titleN <- paste0("All Monogenic Fst Variants >= ", threshold)

  bar_plot_2 <- ggplot(highFst_count_perPop[order(highFst_count_perPop$Super_Population),],
                       aes(x = rownames(highFst_count_perPop),y = Count, fill = Super_Population)) +
    geom_bar(stat = 'identity',color = "black") +
    labs(title = titleN,
         x = "Sub-population",
         y = "Count") +
    theme_minimal(base_size = 14)

  fileN <- paste0(output_fileName, ".png")
  ggsave(filename = fileN, scale = 3, device = 'png', bg = 'white')
}

thresholds <- c(.50, .4, .3, .2, .1)
for(val in thresholds){
  fName <- gsub('\\.', '_', paste0('barPlot_monogenic_AllAbove', val))
  plotFstAboveThreshold_bar_monogenic_all(monogenicFst_all, val, fName)

}


plotFstAboveThreshold_bar_monogenic_disease <- function(fstData, threshold, output_fileName){

  highFst_count_perPop <- data.frame(apply(fstData, MARGIN = 2, \(x){ sum(x >= threshold)}))
  highFst_count_perPop <- data.frame(highFst_count_perPop, superPops)
  colnames(highFst_count_perPop) <- c("Count", "Sub_Population", "Super_Population")

  titleN <- paste0("Disease Monogenic Fst Variants >= ", threshold)

  bar_plot_2 <- ggplot(highFst_count_perPop[order(highFst_count_perPop$Super_Population),],
                       aes(x = rownames(highFst_count_perPop),y = Count, fill = Super_Population)) +
    geom_bar(stat = 'identity',color = "black") +
    labs(title = titleN,
         x = "Sub-population",
         y = "Count") +
    theme_minimal(base_size = 14)

  fileN <- paste0(output_fileName, ".png")
  ggsave(filename = fileN, scale = 3, device = 'png', bg = 'white')
}

thresholds <- c(.50, .4, .3, .2, .1)
for(val in thresholds){
  fName <- gsub('\\.', '_', paste0('barPlot_monogenic_diseaseAbove', val))
  plotFstAboveThreshold_bar_monogenic_disease(monogenicFst_disease, val, fName)

}

# Boxplot-set graph:


# Trying gpt4's suggestion on how to achieve my desired outcome

head(diseaseFst_popsVSall)

library(tidyverse)

# Reshape the data frame to long format
long_data <- gather(diseaseFst_popsVSall, key = "column_name", value = "value")

# Create a boxplot for each column with column names on the x-axis
box_plot <- ggplot(long_data, aes(x = column_name, y = value)) +
  geom_boxplot() +
  theme_minimal(base_size = 14) +
  labs(title = "Disease Fst Summary - Boxplots Per Population" ,x = "Sub-Population", y = "Fst Value")


# Display the boxplot
print(box_plot)


# Adding super populations for ordering and colorization of boxplots
superPops
colnames(superPops) <- c("sub_population", 'super_population')
colnames(long_data) <- c("sub_population", 'Fst')
nrow(long_data) # 1119007
long_data2 <- merge(long_data, superPops, by='sub_population')

nrow(long_data2) # 1119007

# Order the factor levels of the sub_population variable based on the super_population variable
long_data3$sub_population <- factor(long_data2$sub_population, levels = unique(long_data2$sub_population[order(long_data2$super_population)]))
### ^^ this line was used to order the rows by the super population, which was not working by using DF[order(DF$colToOrderBy),], like we normally would expect


# Create a boxplot for each column with column names on the x-axis
box_plot <- ggplot(long_data3, aes(x = sub_population, y = Fst, fill = super_population)) +
  geom_boxplot() +
  theme_minimal(base_size = 14) +
  labs(title = "Disease Fst Summary - Boxplots Per Population" ,x = "Sub-Population", y = "Fst Value")


# Display the boxplot
print(box_plot)
ggsave(filename = 'DiseaseFst_subPops_boxplots.png', scale = 3, bg = 'white')

long_data3 <- long_data2

# generate top SNPs per population
#
# split DF by col, grab top 20 values per, grab all unique rownames, per rowname create string containing all populations which have occurance of the top 20.
#
# lets try this process with 20, see the final result, and redo with a larder number if needed.


diseaseFst_popsVSall

# split by column
list_diseaseFst_cols <- lapply(names(diseaseFst_popsVSall), function(x){return(data.frame(diseaseFst_popsVSall[, x, drop = F]))})
names(list_diseaseFst_cols) <- colnames(diseaseFst_popsVSall)

# grap top 20 per col ... preserve row names ... gpt4 helped here,

list_diseaseFst_cols_20 <- lapply(list_diseaseFst_cols, function(x){
  sorted_indices <- order(x[, 1], decreasing = TRUE)
  sorted_data <- x[sorted_indices, , drop = FALSE]
  top_20_data <- sorted_data[1:20, , drop = FALSE]
  rownames(top_20_data) <- rownames(sorted_data)[1:20]
  return(top_20_data)
  })
#   ^^^ lots of extra lines to preserve row names.. wonder if there is an option to preserve them with a simpler call

unique_rowNames_20 <- sapply(list_diseaseFst_cols_20, function(x){ return(rownames(x))})
unique_rowNames_20 <- as.character(unique_rowNames_20)
length(unique_rowNames_20) # 620
unique_rowNames_20 <- unique(unique_rowNames_20)
length(unique_rowNames_20) # 300 # that is about what I expected.

#
rowNames_20 <- as.character(sapply(list_diseaseFst_cols_20, function(x){ return(rownames(x))}))
rowNames_20_mat <- sapply(list_diseaseFst_cols_20, function(x){ return(rownames(x))})


# helpful code from gpt4 again... my approach I was conceiving was worse. Didn't think to use which() across the cols of the matrix.. I didn't know that I could but it makes sense.
population_labels <- sapply(unique_rowNames_20, function(row_name) {
  column_indices <- which(rowNames_20_mat == row_name, arr.ind = TRUE)[, 2]
  column_names <- unique(colnames(rowNames_20_mat)[column_indices])
  return(paste(column_names, collapse = "|"))
})
population_labels


# subset disease Fst table with unique names, bind pop_labels

features_of_interest <- c("VariantID", "DISEASE/TRAIT", "MAPPED_GENE", 'CONTEXT', 'INTERGENIC', 'RISK ALLELE FREQUENCY',
                          'MAPPED_TRAIT', 'EnsVar_most_severe_consequence', 'Parent term')

top20variants_perPopulation <- SNPanno_GWAS_ensembl[SNPanno_GWAS_ensembl$VariantID %in% unique_rowNames_20, ]
top20variants_perPopulation <- as.data.frame(top20variants_perPopulation)[,features_of_interest]
length(unique(top20variants_perPopulation$VariantID)) #300 ... we see many non disease traits in here.. Lets filter that out and see if all traits remain.. if not, we know for a fact at some point in this data processing we have an error going on.


diseaseParentTerms <- c('Cancer', 'Other disease', 'Digestive system disorder', 'Cardiovascular disease', 'Neurological disorder', 'Immune system disorder', 'Metabolic disorder')
od <- top20variants_perPopulation[top20variants_perPopulation$`Parent term` %in% diseaseParentTerms,] # down to 586 rows .. from 2758
length(unique(od$VariantID)) # 300, no loss. early filtering good. We can move forward with this table.
top20variants_perPopulation <- od

getwd()
setwd("./workingData/")
save(top20variants_perPopulation, file = "top20variants_perPopulation_table.rds")

# Session end. Time to write our report with all the new figures and the table we have... Gonna have to figure out the best way to display the table.

# gonna make vector of genes for panther, and see if I cant get it to work
setwd("D:/Programming/R_projects/Kulathinal_Lab/GWASpops.pheno2geno/GeneralPlots/presentation_plots/report11_12/")
unique_genes <- unique(top20variants_perPopulation$MAPPED_GENE)
write_lines(unique_genes, file = "top20Variants_associatedGenes.txt") # the lines aren't consistent and clean, need one gene per line

# fixing formatting:

cleaned_uni_genes <- unlist(lapply(unique_genes, function(line){
  words <- strsplit(gsub(",", "", line), "\\s+")[[1]]
  words <- words[words != ""]
  words <- words[words != "-"]
  words <- words[words != "NA"]
  return(words)
}))

length(cleaned_uni_genes) # 354
cleaned_uni_genes <- unique(cleaned_uni_genes)
length(cleaned_uni_genes) # 304

write_lines(cleaned_uni_genes, file = "top20Variants_associatedGenes.txt")

# Now we see about the panther thing


# Getting pictures, figures, tables for report 11-12

colnames(top20variants_perPopulation)
top20variants_perPopulation <- top20variants_perPopulation[,c("VariantID", 'MAPPED_GENE', 'EnsVar_most_severe_consequence','DISEASE/TRAIT', 'MAPPED_TRAIT', 'Parent term')]

# need to bind the pop labels I made to these variants:
length(population_labels)
population_labels[1:20]
length(top20variants_perPopulation[,1]) # 586
class(population_labels)
popLabels <- as.data.frame(population_labels)
popLabels['VariantID'] <- rownames(popLabels)
View(popLabels)
t20 <- merge(top20variants_perPopulation, popLabels, by = 'VariantID')

nrow(t20) # looks great, now we need to see if there is a reasonable way to get the Fst values listed somewhere too..
# worked out Fst by grabbing grid prints of each super population

View(t20)

write.csv(t20, file = 'top_20_variants_per_population.csv')
getwd()
setwd("./workingData/")
save(t20, file = 'top20VarsPerPop_withPopLabels.rds')


# -------------------------------------------------------------------------


# trying to display a bunch of data frames in a grid:


# rounding values for better look in grid display:
list_diseaseFst_cols_20_rounded <- lapply(list_diseaseFst_cols, function(x){
  sorted_indices <- order(x[, 1], decreasing = TRUE)
  sorted_data <- x[sorted_indices, , drop = FALSE]
  top_20_data <- round(sorted_data[1:20, , drop = FALSE], 4)
  rownames(top_20_data) <- rownames(sorted_data)[1:20]
  return(top_20_data)
})

listOfTables <- lapply(list_diseaseFst_cols_20_rounded, gridExtra::tableGrob)
gridExtra::grid.arrange(grobs = listOfTables, ncol = 3) # not quite. Lets try less
gridExtra::grid.arrange(grobs = listOfTables[1:9], ncol = 3)
gridExtra::grid.arrange(grobs = listOfTables[1:3], ncol = 3)

?grid.arrange

# AFRICAN:
gridExtra::grid.arrange(grobs = listOfTables[c(1,2,4,12, 17,22,23,31)], ncol = 8)

# gridExtra::grid.arrange(grobs = listOfTables[c(17,22,23,31)], ncol = 4)
# getwd()
# ggsave(filename = 'table print 1 test.png', scale = 2, device = 'png', bg = 'white')
# ?gridExtra

# ADMIXED AMERICAN:
gridExtra::grid.arrange(grobs = listOfTables[c(3, 10,24,25,27)], ncol = 5)

# SOUTH ASIAN:
gridExtra::grid.arrange(grobs = listOfTables[c(5,16,19,26,28,29)], ncol = 6)


# EAST ASIAN:
gridExtra::grid.arrange(grobs = listOfTables[c(6,8,9,11,20,21)], ncol = 6)


# EUROPEAN:
gridExtra::grid.arrange(grobs = listOfTables[c(13,14,15,18,30)], ncol = 5)


######################################################################


unique(t20$`DISEASE/TRAIT`) # 196 results
nrow(t20)
length(unique(t20$VariantID))

sum(nchar(t20$population_labels) == 3) # not working ... I have duplicated rows somehow... fuck.

sum(duplicated(t20)) #158 duplicated rows ... they came in when we removed a bunch of features that would have otherwise differentiated the rows .. lets remove and proceed

t20 <- t20[!duplicated(t20),] # down to 428 rows
sum(nchar(t20$population_labels) == 3)# 250

write.csv(t20, file = 'top_20_variants_per_population.csv')
dim(t20)
colnames(t20)
unique(t20[,c(1,7)])
nrow(t20)

t20_uni_vars_to_pops <- t20[nchar(t20$population_labels) == 3,]
nrow(t20_uni_vars_to_pops) # 250
t20_uni_vars_to_pops <- t20_uni_vars_to_pops[,c(1,7)]
View(t20_uni_vars_to_pops)
sum(duplicated(t20_uni_vars_to_pops)) # 74
t20_uni_vars_to_pops <- t20_uni_vars_to_pops[!duplicated(t20_uni_vars_to_pops), ]
nrow(t20_uni_vars_to_pops) # 176
length(unique(t20_uni_vars_to_pops$VariantID)) # 176 .. alirght

write.csv(t20_uni_vars_to_pops, file = 'one2one_top_20_variants_per_population.csv')



# -------------------------------------------------------------------------

mono_with_fst <- monogenic_variants_complete[monogenic_variants_complete$VariantID %in% rownames(monogenicFst_disease) ,] # from 3950 to 1742

length(unique(mono_with_fst$VariantID)) #283

length(unique(mono_with_fst$`DISEASE/TRAIT`)) # 705 ... thats a lot




# 4-19-23 ... Refining the monogenic table. Finding high FST mongenic SNPs --------

# first extract only traits which are associated with disease parent terms. X
# then search disease FST against metapop for top monogenic SNPs per population (found top 10 per pop) X
# then create table for monogenics similar to top 20 variants per population table
# then work on ontological phylogeny tree graph

monogenic_variants_complete
features_of_interest
diseaseParentTerms
objects()

monogenic_variants_disease <- monogenic_variants_complete[monogenic_variants_complete$`Parent term` %in% diseaseParentTerms,]
nrow(monogenic_variants_complete) # 3950
nrow(monogenic_variants_disease) # 653
rm(monogenic_variants_disease)
# -------------------------------------------------------------------------

# after looking within the excel table I am seeing that using parent terms is actually discarding a lot of SNPs which are clearly of interest when looking at the disease / trait. Which tells me that I may just want to stick with my inital list and use the appropriate, more general termonology to refer to any findings.

# -------------------------------------------------------------------------

# object of interest
diseaseFst_popsVSall
allFst_popsVSall
allUniqueMonogenicVars

monogenic_Fst_popsVSall <- allFst_popsVSall[rownames(allFst_popsVSall)%in% allUniqueMonogenicVars,]
length(allUniqueMonogenicVars) # 1111
nrow(monogenic_Fst_popsVSall) # 835  .. Most SNPs made it through, now we look for those 5-10 highest in each population


######## ORGIN CODE

# split by column
list_diseaseFst_cols <- lapply(names(diseaseFst_popsVSall), function(x){return(data.frame(diseaseFst_popsVSall[, x, drop = F]))})

# grap top 20 per col ... preserve row names ... gpt4 helped here,
 <- lapply(list_diseaseFst_cols, function(x){
  sorted_indices <- order(x[, 1], decreasing = TRUE)
  sorted_data <- x[sorted_indices, , drop = FALSE]
  top_20_data <- sorted_data[1:20, , drop = FALSE]
  rownames(top_20_data) <- rownames(sorted_data)[1:20]
  return(top_20_data)
})
##############



list_monogenicFst_cols <- lapply(names(monogenic_Fst_popsVSall), function(x){ return(data.frame(monogenic_Fst_popsVSall[ , x, drop = F]))})

list_monogenicFst_cols_10 <- lapply(list_monogenicFst_cols, function(x){
  sorted_indices <- order(x[, 1], decreasing = TRUE)
  sorted_data <- x[sorted_indices, , drop = FALSE]
  top_10_data <- sorted_data[1:10, , drop = FALSE]
  rownames(top_10_data) <- rownames(sorted_data)[1:10]
  return(top_10_data)
})

class(list_monogenicFst_cols_10[[1]])
View(list_monogenicFst_cols_10[[1]]) # rows are named, good.
names(list_monogenicFst_cols_10) <- names(list_diseaseFst_cols_20)

# creating top monogenic SNPs data.frame
#
rowNames_10 <- as.character(sapply(list_monogenicFst_cols_10, function(x){ return(rownames(x))}))
rowNames_10_mat <- sapply(list_monogenicFst_cols_10, function(x){ return(rownames(x))})
View(rowNames_10_mat)

top10Monogenicvariants_perPopulation <- SNPanno_GWAS_ensembl[SNPanno_GWAS_ensembl$VariantID %in% unique(rowNames_10), ]
features_of_interest
FOI <- c("VariantID", "MAPPED_GENE", "EnsVar_most_severe_consequence", "DISEASE/TRAIT", "MAPPED_TRAIT",  "Parent term" )
top10Monogenicvariants_perPopulation <- as.data.frame(top10Monogenicvariants_perPopulation)[,FOI]



# Need to annotate with population labels
#
population_labels_mono <- sapply(unique(rowNames_10), function(row_name) {
  column_indices <- which(rowNames_10_mat == row_name, arr.ind = TRUE)[, 2]
  column_names <- unique(colnames(rowNames_10_mat)[column_indices])
  return(paste(column_names, collapse = "|"))
})

# need to bind the pop labels I made to these variants:

popLabels_10 <- as.data.frame(population_labels_mono)
popLabels_10['VariantID'] <- rownames(popLabels_10)
View(popLabels)
t10_mono <- merge(top10Monogenicvariants_perPopulation, popLabels_10, by = 'VariantID')

View(t10_mono)

write.csv(t10_mono, file = 'top_10_monogenic_variants_per_population.csv')

# want to see the Fst values within the table as well ideally... so we would need to grab Fst values and format them to match the pop labels.
#
# ^^ this seems pretty annoying to do.. need to digest the popLabels_10 obj against the top10Monogenicvariants_perpopulation object to create a data-munged col
# it is also just necessary for the process to be complete. So lets do it.

View(popLabels_10)

nrow(popLabels_10) # 93

Fst_per_population <- character()
for( row in 1:nrow(popLabels_10) ){
  id <- popLabels_10[row, 'VariantID']
  pops_with_id <- popLabels_10[row, 'population_labels_mono']
  pops_with_id <- strsplit(pops_with_id, '|', fixed = TRUE)[[1]]
  entry <- character()

  for(pop in pops_with_id){
    entry <- c(entry, round(list_monogenicFst_cols_10[[pop]][rownames(list_monogenicFst_cols_10[[pop]]) %in% id, ], 4) ) # grabs the Fst value

    if(length(entry) > 1){ #collapse vector into single string
      entry <- paste(entry, collapse = "|")
    }
  }

  Fst_per_population <- c(Fst_per_population, entry) # add Fst's which correspond to population_labels into a vector
}




tSplit <- strsplit(popLabels_10[1,1], '|', fixed = TRUE) # fixed necessary... as using regular expression for | breaks things.. as '|' is OR in regex. matching before or after the bar
tSplit
print(as.character(tSplit)[1])# ok something about converting to character is creating problems..  HERE WAS THE KEY ISSUE

# this is how we lookup the variant value per population:
list_monogenicFst_cols_10[['GWD']][rownames(list_monogenicFst_cols_10[['GWD']]) %in% 'rs920184', ]

# GOT IT
Fst_per_population[1:10]



# Bind new col to monogenic table -----------------------------------------
#
Fst_per_population_mono <- as.data.frame(Fst_per_population)
Fst_per_population_mono['VariantID'] <- rownames(popLabels_10)
View(Fst_per_population_mono)
t10_mono_with_fst <- merge(t10_mono, Fst_per_population_mono, by = 'VariantID')
View(t10_mono_with_fst)

save(t10_mono_with_fst, file = "./workingData/popsVSmetapop_topVariants/top10_mongenics_table.rds")

write.csv(t10_mono_with_fst, file = "top_10_mono_with_fst.csv")


# COMPLETE monogenic table ------------------------------------------------




# Update non monogenic table with Fst numbers ---------------------------------

list_diseaseFst_cols_20
popLabels
t20

View(list_diseaseFst_cols_20)
View(popLabels)
View(t20)

Fst_per_population <- character()
for( row in 1:nrow(popLabels) ){
  id <- popLabels[row, 'VariantID']
  pops_with_id <- popLabels[row, 'population_labels']
  pops_with_id <- strsplit(pops_with_id, '|', fixed = TRUE)[[1]]
  entry <- character()

  for(pop in pops_with_id){
    entry <- c(entry, round(list_diseaseFst_cols_20[[pop]][rownames(list_diseaseFst_cols_20[[pop]]) %in% id, ], 4) ) # grabs the Fst value

    if(length(entry) > 1){ #collapse vector into single string
      entry <- paste(entry, collapse = "|")
    }
  }

  Fst_per_population <- c(Fst_per_population, entry) # add Fst's which correspond to population_labels into a vector
}

Fst_per_population[1:10] # looks good



Fst_per_population_disease <- as.data.frame(Fst_per_population)
Fst_per_population_disease['VariantID'] <- rownames(popLabels)
View(Fst_per_population_disease)
t20_disease_with_fst <- merge(t20, Fst_per_population_disease, by = 'VariantID')
View(t20_disease_with_fst)

save(t20_disease_with_fst, file = "./workingData/popsVSmetapop_topVariants/top_20_fst_disease_variants.rds")

write.csv(t20_disease_with_fst, file = "top_20_disease_with_fst.csv")



# 4-23-23 Finding top traits to seek incidence for ------------------------

# gonna use a term-counter to see which terms are most highly represented:
#

t20_disease_with_fst
t10_mono_with_fst

mono_topDiseasesTable_diseaseTrait <- table(t10_mono_with_fst$`DISEASE/TRAIT`)
mono_topDiseasesTable_mappedTrait <- table(t10_mono_with_fst$MAPPED_TRAIT)

allDisease_topDiseaseTable_diseaseTrait <- table(t20_disease_with_fst$`DISEASE/TRAIT`)
allDisease_topDiseaseTable_mappedTrait <- table(t20_disease_with_fst$MAPPED_TRAIT)

View(mono_topDiseasesTable_diseaseTrait) # mostly unhelpful terms
View(mono_topDiseasesTable_mappedTrait) # Lets merge the diseases found from separate sources into here
View(allDisease_topDiseaseTable_diseaseTrait) # very clear top 10 coming out here. Wonderful
View(allDisease_topDiseaseTable_mappedTrait)

# merging monogenic tables together..

monogenic_variants_complete
colnames(monogenic_variants_complete)
colnames(t10_mono_with_fst)

thin_MVC <- monogenic_variants_complete[, c('VariantID', 'DISEASE/TRAIT', "MAPPED_TRAIT", 'monogenic_disease_name')]

nu_t10_mono_with_fst <- merge(t10_mono_with_fst, thin_MVC, by.x = c('VariantID', 'DISEASE/TRAIT', 'MAPPED_TRAIT'), by.y = c('VariantID', 'DISEASE/TRAIT', 'MAPPED_TRAIT') )
nu_t10_mono_with_fst <- unique(nu_t10_mono_with_fst) # because duplicates in the original DFs being merged in we have to remove them to produce a meaningful DF

dim(t10_mono_with_fst)
dim(nu_t10_mono_with_fst) # only 18 new rows.. reasonable

View(monogenic_variants_complete)
View(nu_t10_mono_with_fst)

# make table to count occurences of diseases:

mono_topDiseasesTable <- table(nu_t10_mono_with_fst$monogenic_disease_name)

View(mono_topDiseasesTable)

View(allDisease_topDiseaseTable_diseaseTrait)
########## COMPLETE ###############
########## Lists generated from tables, stored in my planning document in google #########



















































