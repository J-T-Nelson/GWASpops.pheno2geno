# Digesting tables for disease prevalence in 31 1kGenomes populations

prevalence_IHME <- read.csv("./disease_prevalence/IHME-GBD_2019_DATA-ba87a9da-1.csv")
colnames(prevalence_IHME)
IHME_countries <- unique(prevalence_IHME$location_name)
IHME_countries

thousGenomes_populations <- pops$Pop_Ancestry[2:32]
thousGenomes_populations

prevalence_PKU <- read.csv("./disease_prevalence/PKU_disease_prevalence_by_country.csv")
prevalence_CF <- read.csv("./disease_prevalence/worldwideRatesCysticFibrosis_J-Guo_A-Garratt.csv")
head(prevalence_PKU)
head(prevalence_CF)

PKU_countries <- unique(prevalence_PKU$World.region...Country)
CF_countries <- unique(prevalence_CF$Country)

PKU_countries
CF_countries


# Modify tables for analysis:
#
# PKU percentage column must be generated:

prevalence_PKU['prevalence_percentage'] <- sapply(prevalence_PKU$Prevalence..1...X.,
                                                  function(x){ return( 1/as.numeric(gsub(',', '', x)) )} )
prevalence_PKU$prevalence_percentage # numbers are quite small, may need some kind of transform to make them more visible/interesting. Perhaps a log transform
View(prevalence_PKU)
prevPKUnoNA <- prevalence_PKU[apply(prevalence_PKU, 1, function(row) any(!is.na(row))), ] # didn't quite do it.
class(prevPKUnoNA[1,1])

prevalence_PKU <- prevalence_PKU[ -c(1, 65, 66, 67, 68) ,] # removal by index is sufficient here. (negative sign!! forgot this for too long)

# CF percent col:

prevalence_CF['prevalence_percentage_estimate'] <- sapply(prevalence_CF$Estimated.CF.prevalence..per.10.000.,
                                                          function(x){ return( as.numeric(x)/10000 )} )
warnings()
View(prevalence_CF)


# -------------------------------------------------------------------------

# filter down IHME table to only relevant rows

unique(prevalence_IHME$metric_name)
prevalence_IHME <- prevalence_IHME[prevalence_IHME$metric_name == 'Percent', ] # from 18k to 6k rows
head(prevalence_IHME)

# Chat GPT code:

library(dplyr)

prevalence_IHME_summary <- prevalence_IHME %>%
  group_by(location_name, cause_name) %>%
  summarise(avg_val = mean(val, na.rm = TRUE)) %>%
  ungroup()

head(prevalence_IHME_summary) # good. Remarkably simple code. Suppose I do need to work on my dplyr skills here and there


# -------------------------------------------------------------------------

# Now we map onto prevalence_IHME_summary the average allele frequency associated with each disease for each population.
# then graph the results

#  XX 1. remove irrelevant populations, XX 2. associate mapped 1kGenomes populations labels with new row (entered as vectors) PER table
#     3. workout how to get allele frequencies for SNP sets per population 4. write function to go through table, average allele frequency for allele set and associate it per row 5. graph each disease showing population average allele frequency next to prevalance for that population (prevalence will need some transform to be compared on the same scale I think. Possibly just normalization?)

# first  CF list
pops$Pop_Ancestry[2:32]

bad_countries <-  c("Austria", "Belgium", "Bulgaria", "Cyprus", "Czechia", "Denmark", "France", "Germany", "Greece", "Hungary", "Ireland", "Latvia", "Lithuania", "Luxembourg", "Malta", "Netherlands", "Poland", "Portugal", "Romania", "Slovakia", "Slovenia", "Sweden", "Estonia", "Switzerland", "Norway", "Ukraine", "Russia", "Albania", "Armenia", "Croatia", "Serbia", "Moldova", "Belarus", "Georgia", "Bosnia and Herzegovina", "Iceland", "Canada", "Guatemala", "Cuba", "Haiti", "Dominican Republic", "Honduras", "Nicaragua", "El Salvador", "Costa Rica", "Panama", "Jamaica", "Trinidad", "Indonesia", "Philippines", "Thailand", "Myanmar", "South Korea", "North Korea", "Uzbekistan", "Malaysia", "Nepal", "Kazakhstan", "Cambodia", "Azerbaijan", "Tajikistan", "Laos", "Kyrgyzstan", "Turkmenistan", "Singapore", "Mongolia", "Israel", "Turkey", "Oman", "Qatar", "Bahrain", "Saudi Arabia", "Iran", "Jordan", "Palestine", "Lebanon", "Syria", "UAE", "Kuwait", "Iraq", "Afghanistan", "Yemen", "Australia", "New Zealand", "Papua New Guinea", "Brazil", "Chile", "Argentina", "Venezuela", "Ecuador", "Bolivia", "Paraguay", "Uruguay", "Egypt", "South Africa", "Algeria", "Libya", "Sudan", "Morocco")
# bad countries looks good now. ChatGPT didn't get things right at first.
sum(!(prevalence_CF$Country %in% bad_countries))

prevalence_CF <- prevalence_CF[!(prevalence_CF$Country %in% bad_countries), ]
prevalence_CF <- prevalence_CF[!(prevalence_CF$Country %in% "North Macedonia"), ]



# PKU bad row removal

bad_countries_pku <- c("Philippines", "Singapore", "Taiwan", "Thailand", "Bahrain", "Iran", "Iraq", "Israel", "Jordan", "Saudi Arabia", "Turkey", "UAE", "Egypt", "Tunisia", "Austria", "Belarus", "Bosnia-Herzegovina", "Bulgaria", "Croatia", "Cyprus", "Czech Republic", "Denmark", "Estonia", "France", "Germany", "Greece", "Hungary", "Iceland", "Ireland", "Latvia", "Lithuania", "Moldova", "Netherlands", "Norway", "Poland", "Portugal", "Romania", "Russia", "Serbia", "Slovak Republic", "Slovenia", "Sweden", "Switzerland", "Ukraine", "Canada", "Australia", "New Zealand", "Argentina", "Brazil", "Chile", "Costa Rica", "Cuba")

prevalence_PKU$World.region...Country <- trimws(prevalence_PKU$World.region...Country) # there was trailing white space messing things up

sum(!(prevalence_PKU$`World.region...Country` %in% bad_countries_pku))

prevalence_PKU <- prevalence_PKU[!(prevalence_PKU$`World.region...Country` %in% bad_countries_pku), ]
prevalence_PKU <- prevalence_PKU[!(prevalence_PKU$World.region...Country %in% 'Korea'), ]


# associating 1k genomes population titles with populations of prevalence data:

IHME_labels <- unique(prevalence_IHME_summary$location_name)
PKU_labels <- unique(prevalence_PKU$World.region...Country)
CF_labels <- unique(prevalence_CF$Country)
IHME_labels
PKU_labels
CF_labels

thousGenAbbreviations <- gsub("1000GENOMES:phase_3:","", pops$Population_Abbreviation[2:32])
thousGenAbbreviations
thousGenomes_populations

fcgpt <- sapply(1:length(thousGenomes_populations), function(x){ paste(thousGenomes_populations[x], " - ", thousGenAbbreviations[x])})
fcgpt


# Create a dictionary to map location names to their corresponding abbreviations
#
# First IHME:
location_to_abbrev <- list(
  Bangladesh = c("BEB"),
  Barbados = c("ACB"),
  China = c("CDX", "CHB", "CHS"),
  Colombia = c("CLM"),
  England = c("GBR"),
  Finland = c("FIN"),
  Gambia = c("GWD"),
  Italy = c("TSI"),
  Kenya = c("LWK"),
  `Mexico City` = c("MXL"),
  Nigeria = c("ESN", "YRI"),
  Peru = c("PEL"),
  `Puerto Rico` = c("PUR"),
  Punjab = c("PJL"),
  Scotland = c("GBR"),
  `Sierra Leone` = c("MSL"),
  Spain = c("IBS"),
  `Sri Lanka` = c("STU"),
  Tōkyō = c("JPT"),
  `Viet Nam` = c("KHV")
)

# Create a list of vectors with 3-letter abbreviations for each location in IHME_labels
IHME_abbrevs <- lapply(IHME_labels, function(x) location_to_abbrev[[x]])

# Create the data.frame
# df1 <- data.frame(labels = IHME_labels, abbrevs = I(unlist(IHME_abbrevs)))
df1 <- data.frame(labels = IHME_labels, abbrevs = I(IHME_abbrevs)) # I() specifies inhibiting coercion. WOW. that is actually amazing.
View(df1)
?I

dim(prevalence_IHME_summary) #200 3

prevalence_IHME_summary <- merge(prevalence_IHME_summary, df1, by.x='location_name', by.y='labels')
dim(prevalence_IHME_summary) # 200 4 ... good merge


# SECOND PKU:

# Create a dictionary to map location names to their corresponding abbreviations
location_to_abbrev2 <- list(
  China = c("CDX", "CHB", "CHS"),
  India = c("BEB", "GIH", "ITU", "PJL", "STU"),
  Japan = c("JPT"),
  Finland = c("FIN"),
  Italy = c("TSI"),
  Spain = c("IBS"),
  `United Kingdom` = c("GBR", "CEU"),
  `United States of America` = c("ASW", "MXL", "PUR"),
  Mexico = c("MXL"),
  Peru = c("PEL")
)

# Create a list of vectors with 3-letter abbreviations for each location in the second set of labels
second_abbrevs <- lapply(PKU_labels, function(x) location_to_abbrev2[[x]])

# Create the data.frame
df2 <- data.frame(labels = PKU_labels, abbrevs = I((second_abbrevs)))

dim(prevalence_PKU) # 10 7
prevalence_PKU <- merge(prevalence_PKU, df2, by.x="World.region...Country", by.y='labels')
dim(prevalence_PKU)# 10 8


# THIRD CF:

# Create a dictionary to map location names to their corresponding abbreviations
location_to_abbrev3 <- list(
  Italy = c("TSI"),
  Spain = c("IBS"),
  UK = c("GBR", "CEU"),
  Finland = c("FIN"),
  USA = c("ASW", "MXL", "PUR"),
  Mexico = c("MXL"),
  India = c("BEB", "GIH", "ITU", "PJL", "STU"),
  China = c("CDX", "CHB", "CHS"),
  Pakistan = c("PJL"),
  Bangladesh = c("BEB"),
  Japan = c("JPT"),
  Vietnam = c("KHV"),
  `Sri Lanka` = c("STU"),
  Colombia = c("CLM"),
  Peru = c("PEL"),
  Nigeria = c("ESN", "GWD", "LWK", "MSL", "YRI")
)

# Create a list of vectors with 3-letter abbreviations for each location in the last set of labels
last_abbrevs <- lapply(CF_labels, function(x) location_to_abbrev3[[x]])

# Create the data.frame
df3 <- data.frame(labels = CF_labels, abbrevs = I(last_abbrevs))

dim(prevalence_CF) # 16 8
prevalence_CF <- merge(prevalence_CF, df3, by.x='Country' , by.y='labels')
dim(prevalence_CF) # 16 9 ... Good


# -------------------------------------------------------------------------

# 3. workout how to get allele frequencies for SNP sets per population

# get SNP sets per disease, coming from monogenic and top20disease tables
#
# then get all 1k genomes allele frequencies on those alleles from stored data structure..

t20_disease_with_fst
nu_t10_mono_with_fst

View(nu_t10_mono_with_fst)

CF_snps <- c("rs975722", "rs6977665", "rs12188164", "rs2036100")

PKU_snps <- nu_t10_mono_with_fst$VariantID[nu_t10_mono_with_fst$monogenic_disease_name %in% c("Phenylketonuria", "Phenylketonuria (PKU)/Phenylalanine hydroxylase deficiency (PAH deficiency)")]
PKU_snps <- unique(PKU_snps)
PKU_snps

# getting the rest of the diease associated SNPs will be less nuanced, all coming from t20_disease..

dNames <- c('Type 2 diabetes', "Schizophrenia", "Hypertension", "Inflammatory bowel disease", "Breast cancer", "Colorectal cancer", "Multiple sclerosis", "Psoriasis", "Asthma", "Prostate cancer")

t20_prevalenceDiseases <- t20_disease_with_fst[t20_disease_with_fst$`DISEASE/TRAIT` %in% dNames, ]
unique(t20_prevalenceDiseases$`DISEASE/TRAIT`) # looks good.


# GPT4 Code:
#
# Split the data.frame by the 'DISEASE/TRAIT' column
split_data <- split(t20_prevalenceDiseases, t20_prevalenceDiseases$`DISEASE/TRAIT`)

# Initialize an empty list to store the named vectors
variant_id_vectors <- list()

# Loop through each split data.frame and extract unique values of the 'VariantID' column
for (disease_trait in names(split_data)) {
  unique_variant_ids <- unique(split_data[[disease_trait]]$VariantID)
  variant_id_vectors[[disease_trait]] <- unique_variant_ids
}

# Set the names of the list elements to be the same as the unique values in the 'DISEASE/TRAIT' column
names(variant_id_vectors) <- names(split_data)


## 4. write function to go through table, average allele frequency for allele set and associate it per row
## 5. graph each disease showing population average allele frequency next to prevalence for that population (prevalence will need some transform to be compared on the same scale I think. Possibly just normalization?)
##
##  ^^ I have the SNPs, I just need to load in the master data structure containing perSNP snp frequencies and average the per population frequency of all SNPs for each disease, then associate those values back to the tables I have been creating.

# then graphing. # and creating a table of correlations per disease per population between prevalence and fst frequency.


# Session 2 start 5-3 -----------------------------------------------------

ls1 <- ls()
load("./workingData/full_data_for_analysis/full_popAlleleFreq_PerAllele.rds")
ls2 <- ls()
ls2[!ls2%in%ls1]


length(full_popAlleleFreq_perAllele) #213,710
View(full_popAlleleFreq_perAllele)

# data we want is in here, will just have to figure out how to most effectively retrieve results.

View(variant_id_vectors)
View(PKU_snps)
View(CF_snps)

# combining all snps per disease

top_disease_snps <- variant_id_vectors
class(top_disease_snps)
top_disease_snps[["PKU"]] <- PKU_snps
top_disease_snps[["CF"]] <- CF_snps
View(top_disease_snps) # looks good.


# -------------------------------------------------------------------------

class(full_popAlleleFreq_perAllele[[1]])


# standardizing PKU and CF tables:

std_prevalence_PKU <- prevalence_PKU[,c(1,7,8)] #grab desired columns
std_prevalence_PKU['cause_name'] <- rep("Phenylketonuria", nrow(std_prevalence_PKU)) # add new column with disease name
colnames(std_prevalence_PKU) <- c("location_name", "avg_val", "abbrevs", "cause_name") #rename cols to standard names
std_prevalence_PKU <- std_prevalence_PKU[,c(1,4,2,3)] # readjust order


std_prevalence_CF <- prevalence_CF[, c(1,8,9)]
std_prevalence_CF['cause_name'] <- rep("Cystic Fibrosis", nrow(std_prevalence_CF))
colnames(std_prevalence_CF) <- c("location_name", "avg_val", 'abbrevs', 'cause_name')
std_prevalence_CF <- std_prevalence_CF[,c(1,4,2,3)]

#combine tables now
prevalence_all <- dplyr::bind_rows(prevalence_IHME_summary, std_prevalence_CF, std_prevalence_PKU)
View(prevalence_all) # great, working with this will be much easier now.

# table standardization complete ------------------------------------------

# per disease make DF where cols are pops, rows are snps,
# per col average and insert into new DF where cols are pops single row with averages of a diseases SNPs per population
# 3 tables get different treatment from here: add new data into each table.. ( would be easier if tables had standardized forms )

top_disease_snps

# making DF per disease

DF_list_ref <- list()
disease_names <- names(top_disease_snps)

for (i in seq_along(top_disease_snps)) { # GPT4 helped make this
  disease_name <- disease_names[i]
  disease <- top_disease_snps[[i]]

  DF_list_ref[[disease_name]] <- list()
  for (snpID in disease) {
    DF_list_ref[[disease_name]][[snpID]] <- full_popAlleleFreq_perAllele[[snpID]]

  }
}


names(top_disease_snps[1])


View(DF_list)
?grepl


tmod <- DF_list[[1]][[1]]
tmod <- as.data.frame(tmod[grepl("1000GENOMES", tmod$population), ])
tmod['population'] <- gsub("1000GENOMES:phase_3:",'', tmod$population)
attr(tmod, "Ancestral_Allele") <- calc_ancestralAllele(tmod)
attributes(tmod)
tmod <- tmod[!tmod$allele %in% attr(tmod, "Ancestral_Allele"),]
row.names(tmod) <- tmod$population
gudDF <- tmod[,"frequency", drop = FALSE]
class(gudDF)

View(gudDF)

is.na(attr(tmod, "Ancestral_Allele") |)


DF_list <- list()
disease_names <- names(top_disease_snps)
item_count = 0

for (i in seq_along(top_disease_snps)) { # GPT4 helped make this
  disease_name <- disease_names[i]
  disease <- top_disease_snps[[i]]

  DF_list[[disease_name]] <- list()
  for (snpID in disease) {
    item_count <- item_count + 1
    # print(item_count)
    # print(snpID)

    tempDF <- full_popAlleleFreq_perAllele[[snpID]]
    tempDF <- as.data.frame(tempDF[grepl("1000GENOMES", tempDF$population), ]) # selecting rows of interest.. 1k genomes rows
    tempDF['population'] <- gsub("1000GENOMES:phase_3:",'', tempDF$population) # renaming to only abbreviations
    # resetting ancestral alleles if needed
    if( is.na(attr(tempDF, "Ancestral_Allele")) | !(attr(tempDF, "Ancestral_Allele")%in% unique(tempDF$allele) )){
      attr(tempDF, "Ancestral_Allele") <- calc_ancestralAllele(tempDF)
    }

    tempDF <- tempDF[!tempDF$allele %in% attr(tempDF, "Ancestral_Allele"),] # removing ancestral alleles, as they're not desirable
    row.names(tempDF) <- tempDF$population # setting names for later DF stacking ... having errors too often
    tempDF <- tempDF[,"frequency", drop = FALSE] #grabbing only frequency col

    DF_list[[disease_name]][[snpID]] <- tempDF
  }
}

# FINALLY we have the desired result.. fucking impossible data to make assumptions about....

View(DF_list_ref)



# smash DF_list together by disease ---------------------------------------
# then take average per population
# then attach to other data frame
# then graph
# then correlation

?eval.parent
?substitute

DF_list_flat <- list()
for (i in seq_along(DF_list)) {
  disease_name <- disease_names[i]


  for(j in seq_along(DF_list[[i]])){

  }
}


combined_DFs <- lapply(DF_list, function(x){
  combined_df <- do.call(rbind, x)
  return(combined_df)
  })
View(combined_DFs) # no good.

browTest <- dplyr::bind_cols(DF_list[[1]])
View(browTest) # no good.


as.character(unique(purrr::flatten(purrr::flatten(unique(prevalence_all$abbrevs)))))
# [1] "BEB" "ACB" "CDX" "CHB" "CHS" "CLM" "GBR" "FIN" "GWD" "TSI" "LWK" "MXL" "ESN" "YRI" "PEL" "PUR" "PJL" "MSL" "IBS" "STU" "JPT" "KHV" "GIH" "ITU"
# [25] "CEU" "ASW"

pop_labels


# GPT4 code:
#
# Function to process each sublist
process_sublist <- function(sublist, pop_labels) {
  pop_means <- lapply(pop_labels, function(pop) {
    # Extract rows with matching population labels
    matching_rows <- lapply(sublist, function(df) {
      df[df$name == pop, ]
    })

    # Calculate mean value
    mean_value <- mean(unlist(lapply(matching_rows, function(row) row$frequency)), na.rm = TRUE)
    return(mean_value)
  })

  # Create a data.frame with the mean values
  result_df <- data.frame(name = pop_labels, mean_value = unlist(pop_means), row.names = pop_labels)
  return(result_df)
}

# Process the list of data.frames
output_list <- lapply(DF_list, process_sublist, pop_labels = pop_labels)


# attempt 2 from GPT4 ...

# Function to process each sublist
process_sublist <- function(sublist, pop_labels) {
  pop_means <- lapply(pop_labels, function(pop) {

    # Extract rows with matching population labels
    matching_rows <- lapply(sublist, function(df) {
      df[df$name == pop, ]
    })
    print(class(matching_rows))
    print(class(matching_rows[[1]]))
    print(matching_rows[[1]][1])
    # Calculate mean value
    mean_value <- mean(unlist( lapply(matching_rows, function(row) row)), na.rm = TRUE)
    return(mean_value)
  })

  # Create a data.frame with the mean values
  result_df <- data.frame(name = pop_labels, mean_value = unlist(pop_means), row.names = pop_labels)
  return(result_df)
}

# Process the list of data.frames
output_list <- lapply(DF_list, process_sublist, pop_labels = pop_labels)


# returning all NaN ... not sure where the code is fucking up. Getting too tired to keep this up


# Going to try my own method.
length(pop_labels)

fixed_DF_list <- lapply(DF_list, function(subList){
  for(i in seq_along(subList)){
    tempDF <- subList[[i]]
    tempDF <- tempDF[!(rownames(tempDF) == "ALL"), , drop=F] # remove all value
    print(class(tempDF))

    for(pop_lab in pop_labels){
      if( !(pop_lab %in% rownames(tempDF)) ){
        tempDF[pop_lab, ] <- 0
      }
    }
    subList[[i]] <- tempDF
  }
  return(subList)
})

# ok the data frames are now consistent ... we can just ignore the 0's when taking means I think

# -------- DATA MADE CONSISTENT FOR EASE OF PROCESSING -----------------------------------------------------------------

# Bind DF's together now.

bound_DF_list <- lapply(fixed_DF_list, function(x){ dplyr::bind_cols(x)})
View(bound_DF_list)


## Function to take mean:

mean_0 <- function(numeric_vec){
  #ignores 0 values for taking mean

  # Remove 0 values from the vector
  non_zero_values <- numeric_vec[numeric_vec != 0]
  mean_value <- mean(non_zero_values)

  return(mean_value)
}


# attempt at taking means...
#
mean_frequency_DF_list <- lapply(bound_DF_list, function(sublist){
  temp_means <- numeric()
  for(i in 1:length(nrow(sublist))){
    temp_means[i] <- mean_0(sublist[i, ])
  }
  names(temp_means) <- rownames(sublist)
  return(temp_means)
})


# GPT4 helping..
mean_frequency_DF_list <- lapply(bound_DF_list, function(sublist) {

  temp_means <- apply(sublist, 1, function(row) { # apply is good for applying functions to data frames by col or row..
    mean_0(row)
  })

  # Create a new data frame with the mean values and row names
  result_df <- data.frame(mean_value = temp_means, row.names = rownames(df))

  return(result_df)
})


# Attach means to other data frame now ------------------------------------

# first remove rows where avg_val is NA in prevalence_all
# use prevalence_all to query against mean_frequency_DF_list ... PER ROW: 2 processing paths, always returning a single number per query (a mean) get a mean of means when there is more than one population within 'abbrevs' col for a row.


# do this then graph. then correlation then paper writing hard mode.


# 5-6-23 session start  ---------------------------------------------------------------------


# Removing NA rows wrt avg_val: 2 NAs for cystic fibrosis nigeria

sum(is.na(prevalence_all$avg_val)) # 2
dim(prevalence_all) # 226 4
prevalence_all <- prevalence_all[-which(is.na(prevalence_all$avg_val)), ]
dim(prevalence_all) # 224 4


lapply(mean_frequency_DF_list, summary) # means range from .33-.55 or so.
mean(prevalence_all$avg_val) # 0.0113 .. very small.

# summary(purrr::flatten(purrr::flatten(mean_frequency_DF_list)))

# mean(mean_frequency_DF_list[[1]][,1]) # just learned you cannot coerce a 1 col data frame into a numeric with as.numeric... idk why.

mean(sapply(mean_frequency_DF_list, function(x){mean(x[,1])})) # done, OVERALL MEAN: 0.4268217

# The means are quite different in each set, so I think per disease, we will make two graphs, one graph for prevalence, and one graph for frequency, the two bar graphs can be compared side by side for each disease. (with ranges automatically adjusting I think this will be out best bet for visual comparison)
#

## Setup graphing function:

# split up prevalence DFs, per each one query respective mean freq DF for values, store both in variables generate two plots and save them both, name according to arguments provided.

prev_DFs <- split(prevalence_all, prevalence_all$cause_name)
View(prev_DFs)




graph_prevalence_vs_meanFreq <- function(prev_DF, mean_freq_list){
  # pass single prev_df in and the whole mean_freq_list

  current_disease_name <- unique(prev_DF$cause_name)

  # get list of pops
  pop_abbrevs <- prev_DF$abbrevs
  print(pop_abbrevs)

  # get frequency means for respective pops
  aggregate_allele_freqs <- numeric()
  for(list in pop_abbrevs){
    if(length(list) == 1){
      mean_freq <- mean_freq_list[[current_disease_name]][list,] # grabbing value by row name ref
      aggregate_allele_freqs <- append(aggregate_allele_freqs, mean_freq)
    }
    else{ # dealing with cases where multiple populations map onto a given diseases location category
      temp_mean <- numeric()
      for(label in list){
        mean_freq <- mean_freq_list[[current_disease_name]][label,]
        temp_mean <- append(temp_mean, mean_freq)
      }
      mean_freq <- mean(temp_mean)
      aggregate_allele_freqs <- append(aggregate_allele_freqs, mean_freq)
    }
  }

  #generate label vector for populations, (combining multiple populations where necessary)
  allele_freq_labels <- sapply(pop_abbrevs, function(labs){
    if(length(labs) == 1){ return(labs)}
    else{return(paste(labs, collapse = "|"))}
  })
  print(allele_freq_labels)

  # generate super pop DF
  super_pops <- bind_super_pops(pop_abbrevs, superPops)
  super_pops <- cbind(super_pops, allele_freq_labels)
  print(super_pops)

  prev_DF['super_population'] <- super_pops$super_population # adding super population column to be used for ordering of graphed samples

  # Generate and save two plots
  # 1. allele_freq_labels x aggregate_allele_freqs
  # 2. prev_DF$location_name x prev_DF$avg_val

  # Set up side-by-side layout
  par(mfrow = c(1, 2))

  # Plot 1: Allele frequency labels x Aggregate allele frequencies
  plot1 <- ggplot() +
    geom_point(aes(x = super_pops$allele_freq_labels, y = aggregate_allele_freqs), size = 4) +
    labs(title = paste("Mean Allele Frequency of", current_disease_name), x = "Population", y = "Aggregate Allele Frequency") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  print(plot1)

  # Plot 2: Location names x Prevalence
  plot2 <- ggplot() +
    geom_point(aes(x = prev_DF$location_name, y = prev_DF$avg_val), size = 4) +
    labs(title = paste("Prevalence of", current_disease_name), x = "Location", y = "Prevalence") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  print(plot2)

}

mean(c(mean_frequency_DF_list[['Asthma']]['CDX',], mean_frequency_DF_list[['Asthma']]['CHB',], mean_frequency_DF_list[['Asthma']]['CHS',]))

# helper func for graphing func, returns the superpopulations which would be associated with a provided abbrvs list
bind_super_pops <- function(pop_abbr, super_pops_DF){
  proxyLabs <- data.frame(sub_population = sapply(pop_abbr, function(x){return(x[1])}))

  # Add a temporary column with the original row numbers
  proxyLabs$row_num <- 1:nrow(proxyLabs)

  tempMerge <- merge(proxyLabs, super_pops_DF)

  # Sort the merged data frame based on the original row numbers
  tempMerge <- tempMerge[order(tempMerge$row_num), ]

  # Remove the temporary row number column
  tempMerge$row_num <- NULL

  # Reset row names to be the natural order of the data frame (now indices of the original DF are preserved wrt the values in the rows)
  rownames(tempMerge) <- 1:nrow(tempMerge)

  return(tempMerge)
}

bspTordered <- bind_super_pops(prev_DFs$Asthma$abbrevs, superPops)
View(bspTordered)
bspT2ordered <- cbind(bspTordered, prev_DFs$Asthma$abbrevs)
View(bspT2ordered)




?cbind

# manual confirmation of mean of means for combined population means ... looks good.


# Try out graphing function:
#



graph_prevalence_vs_meanFreq(prev_DF = prev_DFs$Asthma, mean_freq_list = mean_frequency_DF_list)


superPops
class(superPops)



# New graphing func from GPT4 ---------------------------------------------


graph_prevalence_vs_meanFreq_bar <- function(prev_DF, mean_freq_list, cur_disease_name = ""){
  # pass single prev_df in and the whole mean_freq_list

  current_disease_name <- cur_disease_name #unique(prev_DF$cause_name)

  # get list of pops
  pop_abbrevs <- prev_DF$abbrevs

  # get frequency means for respective pops
  aggregate_allele_freqs <- numeric()
  for(list in pop_abbrevs){
    if(length(list) == 1){
      mean_freq <- mean_freq_list[[current_disease_name]][list,] # grabbing value by row name ref
      aggregate_allele_freqs <- append(aggregate_allele_freqs, mean_freq)
    }
    else{ # dealing with cases where multiple populations map onto a given diseases location category
      temp_mean <- numeric()
      for(label in list){
        mean_freq <- mean_freq_list[[current_disease_name]][label,]
        temp_mean <- append(temp_mean, mean_freq)
      }
      mean_freq <- mean(temp_mean)
      aggregate_allele_freqs <- append(aggregate_allele_freqs, mean_freq)
    }
  }

  #generate label vector for populations, (combining multiple populations where necessary)
  allele_freq_labels <- sapply(pop_abbrevs, function(labs){
    if(length(labs) == 1){ return(labs)}
    else{return(paste(labs, collapse = "|"))}
  })

  # generate super pop DF
  super_pops <- bind_super_pops(pop_abbrevs, superPops)
  super_pops <- cbind(super_pops, allele_freq_labels, aggregate_allele_freqs)
  #super_pops <- super_pops[order(super_pops$super_population), ] # reorder for graphing
  #print(super_pops)

  prev_DF['super_population'] <- super_pops$super_population # adding super population column to be used for ordering of graphed samples
  #prev_DF <- prev_DF[order(prev_DF$super_population), ] # reorder for graphing

  # Generate and save two plots
  # 1. allele_freq_labels x aggregate_allele_freqs
  # 2. prev_DF$location_name x prev_DF$avg_val

  # Plot 1: Allele frequency labels x Aggregate allele frequencies
  plot1 <- ggplot(super_pops) +
    geom_bar(aes(x = allele_freq_labels, y = aggregate_allele_freqs,
                 fill = super_population), stat = "identity") +
    labs(title = paste("Mean Allele Frequency of", current_disease_name), x = "Population (1000Genomes)", y = "Aggregate Allele Frequency") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  print(plot1)

  plot1Name <- paste("meanAlleleFreq_", current_disease_name, ".png")
  ggsave(plot1Name, scale = 1.75, device = "png", bg = "white")

  # Plot 2: Location names x Prevalence
  plot2 <- ggplot(prev_DF) +
    geom_bar(aes(x = location_name, y = avg_val, fill = super_population), stat = "identity") +
    labs(title = paste("Prevalence of", current_disease_name), x = "Location of Sample Collection", y = "Prevalence") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  print(plot2)

  plot2Name <- paste("prevalence_", current_disease_name, ".png")
  ggsave(plot2Name, scale = 1.75, device = "png", bg = "white")

  return(cor(prev_DF$avg_val, aggregate_allele_freqs)) # returning correlation coefficient for numeric interpretation.

}



graph_prevalence_vs_meanFreq_bar(prev_DF = prev_DFs$Asthma, mean_freq_list = mean_frequency_DF_list)
debug(graph_prevalence_vs_meanFreq_bar)
getwd()

setwd("./GeneralPlots/final_report_plots/")


correlation_coefficients <- sapply(prev_DFs, function(x){graph_prevalence_vs_meanFreq_bar(prev_DF = x, mean_frequency_DF_list)})

# error for colon and rectum cancer.. going to remove from list

mod_prev_DFs <- prev_DFs[!(names(prev_DFs) %in% c("Colon and rectum cancer", "Cystic Fibrosis"))]

correlation_coefficients <- sapply(mod_prev_DFs, function(x){graph_prevalence_vs_meanFreq_bar(prev_DF = x, mean_frequency_DF_list)})

debug(graph_prevalence_vs_meanFreq_bar)
graph_prevalence_vs_meanFreq_bar(prev_DFs$`Colon and rectum cancer`, mean_frequency_DF_list)

# name inconsistencies creating issues.

names(mean_frequency_DF_list) <- c("Asthma", "Breast cancer", "Colorectal cancer", "Hypertension", "Inflammatory bowel disease", "Multiple sclerosis", "Prostate cancer", "Psoriasis", "Schizophrenia", "Type 2 diabetes","Phenylketonuria", "Cystic Fibrosis")

names(prev_DFs) <- c("Asthma", "Breast cancer", "Colorectal cancer", "Cystic Fibrosis", "Type 2 diabetes", "Hypertension", "Inflammatory bowel disease", "Multiple sclerosis", "Phenylketonuria", "Prostate cancer", "Psoriasis", "Schizophrenia")


undebug(graph_prevalence_vs_meanFreq_bar)
correlation_coefficients <- sapply(prev_DFs, function(x){graph_prevalence_vs_meanFreq_bar(prev_DF = x, mean_frequency_DF_list)})


debug(graph_prevalence_vs_meanFreq_bar)
graph_prevalence_vs_meanFreq_bar(prev_DFs$`Colorectal cancer`, mean_frequency_DF_list)

correlation_coefficients <- sapply(1:length(prev_DFs), function(i){
  graph_prevalence_vs_meanFreq_bar(prev_DF = prev_DFs[[i]], mean_frequency_DF_list, cur_disease_name = names(prev_DFs)[i])})

correlation_coefficients_DF <- data.frame(correlation_coefficients, disease = names(prev_DFs))
View(correlation_coefficients_DF)

#reorder and export correlation coefficient DF

correlation_coefficients_DF <- correlation_coefficients_DF[, c(2,1)]
correlation_coefficients_DF <- correlation_coefficients_DF[-c(11,12), ]

write.csv(correlation_coefficients_DF, file = "diseasePrevalence_vs_meanSNP_frequency.csv", row.names = F)

?write.csv
# COMPLETE ----------------------------------------------------------------













































































