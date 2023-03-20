# Fst Analysis 1: top 1%, 5%, 10% across all popPairs, and within individual popPairs.

all_popPairs_sum <- sapply(1:nrow(master_fst), function(i){
  return(sum(master_fst[i, ], na.rm = TRUE))
}) # this call is taking a long time to run... so going to have to run it later.

# I may have to do some kind of weighting involving counting up the # NA per row and multiplying by some factor related to that number... Otherwise rows with complete entry sets are automatically weighted more heavily... I could also instead fill NA slots with an average value... hmmm... will have to think about what is most well suited to this task.

thousGenPops

?grep


# create by population list of Fst data.frames
# then create separate list of data.frames which posseses the sum's of each row across those small DFs

popNames <- thousGenPops$Population_Abbreviation

full_fst_perPop <- lapply(popNames, function(x) {
  popFsts <- master_fst[ , grep(x, colnames(master_fst))]
})

names(full_fst_perPop) <- popNames

# save the novel data structure

getwd()
save(full_fst_perPop, file = 'full_fst_perPop.rds')


# filter by top 1, 5 and 10%, create DF's containing those values..

plot(all_popPairs_sum) # very basic plot actually looks interesting. Some clear stand outs.. though we will have to do some kind of normalization to figure out the true story here.

boxplot(all_popPairs_sum)


dens <- density(all_popPairs_sum)

plot(dens, main = "Density of population pair Fst values") # this is roughly what we're looking for. and honestly it looks more similar to expectation of a random distribution than I would have anticipated.


summary(all_popPairs_sum) # min looks good.


# Generating Better figures for reporting and presentation.
#
#   Starting with normalizing data

# finding NA rows to find rsIDs with missing 1k genomes data.. want to understand why there are missing values.

naVals <- is.na(master_fst)

rows_with_na <- apply(master_fst, 1, function(row) any(is.na(row)))
rows_with_na[1:10]
onlyPosRows <- rows_with_na[rows_with_na == TRUE] # 57575 rows.
57575/159559 # 0.3608383 ... so 36% of rows / rsIDs contain incomplete data. Meaning the normalization should have a strong affect on what comes to the top.

onlyPosRows[1:20]
# rs35472707  rs78444298   rs3850625 rs111827672   rs3790604  rs78793397  rs78018265  rs17046344 rs138088171   rs4722675    rs664223 rs193123434 rs139235514
# TRUE        TRUE        TRUE        TRUE        TRUE        TRUE        TRUE        TRUE        TRUE        TRUE        TRUE        TRUE        TRUE
# rs139941554  rs12686364  rs72754495 rs185373991 rs181468389 rs147833373 rs148767741
# TRUE        TRUE        TRUE        TRUE        TRUE        TRUE        TRUE

alleles_to_check <- names(onlyPosRows)[1:10]
check_alleles <- full_popAlleleFreq_perAllele[alleles_to_check]

num_missing_pops <- sapply(check_alleles, function(x) {
  thousGenPops$Population_Abbreviation %in% x$population
})
# alright, we actually see that ALL 1k genomes pops exist in all 10 of these.. I know for a fact I have seen some examples where that is not the case, however it seems the majority of cases (alleles / variantIDs) have all 1k genomes pops.
# What we ARE seeing however is that when an allele is totally fixed, (which according to the limited number of samples checked) we don't see reporting for both alleles, as the rest of the data is necessarily inferred. ... if we wanted to get this analysis more precisely correct, we would be filling in those missing values before calculating Fst.. as I believe that missing minor allele frequency in the data is not being properly dealt with by my functions doing calculations. (LETS CHECK THIS ASSUMPTION ABOUT MY CALCULATION)

# OK figured something out: we are only discarding major alleles... so when the major is 1.0, we get no minor allele, and pops get dropped... (COULD AND SHOULD CHECK TO SEE IF ANY ROWS WITH NA VALUES ARE MISSING AND OF THE 1kGenomes POPS... this would give us more info on whats going on in our data), If we assume all NA are from having 1.0 in the major allele however I don't think there is a meaningful way to fill the NA values.. as the FST will come from a comparison which was never made... making it the average would not be good and making it 0 wouldn't be either, as those are just inaccruate, as it could be basically anything I belive. (anything from 0-1) ... so I think the best way to work with the data now would be doing a weighting on the current values by giving relative amplification for the amount of values missing.. just to get a better picture of the relative impacts for now...
#
# WE WILL WANT TO RECALCULATE THESE VALUES IN THE FUTURE BY DOING A BIT OF PREPROCESSING WITHIN THE CALCULATE FST FUNC TO GET BETTER / MORE ACCURATE RESULTS THOUGH... NOT TIME FOR THAT TODAY.

# NOW FOR NORMALIZATION:
#
# count # NA in each row...
# create factor to multiple sum'd value by... do this formula: (total # pop pairs)/((total # pop pairs) - (# missing pairs) )
#  (# missing pairs) == num NA

num_na_perRow <- apply(naVals, 1, sum)
num_na_perRow[1:100]
summary(num_na_perRow) # max and min are what we expect to see.

# for now we should remove any data that are too muddled. i.e. > 50% na
# 496/2 == 248

bad_rows = num_na_perRow >247
sum(bad_rows) # 28518
28518/159559 # 17.9 % rows are bad ... meaning our current reports are somewhat subject to change, but we know why and what to do to make corrections.


# REMOVING BAD ROWS:

bad_rows[1:100]

good_fst_rows <- master_fst[!bad_rows, ] # 131041 remaining observations.

# now to normalize

num_na_perRow_goodRows <- num_na_perRow[!bad_rows]

all_popPairsSum_normalized <- sapply(1:nrow(good_fst_rows), function(i){
  return(sum(good_fst_rows[i, ], na.rm = TRUE))
})

# naming 'fixed' because its a sort of inverse normalization.. and not a standard type of normalization
denom = (all_popPairsSum_normalized - num_na_perRow_goodRows)
fixing_factor = 496/(496 - num_na_perRow_goodRows)

all_popPairsSum_fixed <- all_popPairsSum_normalized/denom

summary(denom) # this doesn't look right. FIGURE THIS OUT
summary(num_na_perRow_goodRows) # looks fine. .. Their lengths are looking good too. Maybe the order is somehow scrambled in one? That would explain seeing the min where its at ... NO we should have

# think I just did the calculation wrong outright
summary(fixing_factor) # yeah this looks good. Mostly 1 values.. with max of 1.96.
all_popPairsSum_fixed <- all_popPairsSum_normalized*fixing_factor

summary(all_popPairsSum_fixed) # max is actually the same as before!
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
# 0.00   13.20   24.69   31.12   42.84  197.01

# the before fixing results:
summary(all_popPairs_sum)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
# 0.000   6.346  18.376  24.943  36.467 197.015




# with values fixed, we can now generate figures and tables
# first show a density chart.
# Then show the scatter displaying the standout values
# then show the top 10, 5 , and 1% of SNPs... (how to display this though? )

# MAKE PLOTS FIRST, BEAUTIFY LATER

dens <- density(all_popPairsSum_fixed)

plot(dens, main = "Density of pop-pair Fst - Fixed") # Tail is a little more fat and curve is a bit more smooth than before

plot(all_popPairsSum_fixed) # we want to see some quantile lines.

# getting quantiles of interest 1, 5, 10


q1 <- quantile(all_popPairsSum_fixed, 0.01) # at 2.1 I think I am getting the smallest 1% not the largest
q5 <- quantile(all_popPairsSum_fixed, 0.05)
q10 <- quantile(all_popPairsSum_fixed, 0.1)
q1
q10


q99 <- quantile(all_popPairsSum_fixed, 0.99)
q99 # 109.3 ... Good
q95 <- quantile(all_popPairsSum_fixed, 0.95)
q90 <- quantile(all_popPairsSum_fixed, 0.90)
q95 # 79.57039
q90 # 64.51429

?abline

# adding title's labels

plot(all_popPairsSum_fixed, main = "pop-pair Fst Sums - Fixed", ylab = "Sum of Variant Fst", xlab = "Index of Sample", pch = 20, cex = .7)

abline(h = q99, col = 'green', lwd = 2)
abline(h = q95, col = 'blue', lwd = 2)
abline(h = q90, col = 'red', lwd = 2)

text(-4000, y = q99, "1%", pos = 3, offset = .5, col = 'darkgreen')
text(-4000, y = q95, "5%", pos = 3, offset = .5, col = 'darkblue')
text(-4000, y = q90, "10%", pos = 3, offset = .5, col = 'darkred')

?text
?par
colors()

# looks great now. Lots of work to get it done... but maybe faster in the future with new knowledge.


plot(dens, main = "Density of pop-pair Fst - Fixed")

abline(v = q99, col = 'green', lwd = 2)
abline(v = q95, col = 'blue', lwd = 2)
abline(v = q90, col = 'red', lwd = 2)
legend("topright", inset = 0.02, legend = sumSts, bg = "white", cex = .80, text.font = 4, box.lwd = 0)


summary_stats <- c("Mean" = mean(all_popPairsSum_fixed),
                   "Median" = median(all_popPairsSum_fixed),
                   "Standard Deviation" = sd(all_popPairsSum_fixed))

sumSts <- c("Mean: 3.12", "Median: 24.69", "Standard Deviation: 23.91", "Min: 0.00", "Max: 197.01")
summary_stats

# How to proceed? What is most important to report on from here? Find which of these 10% are disease related? (such filtering may be quite difficult, as there is no column specifying the nature of the trait in this categorical way yet.)



# seeking col that may offer insight into which traits are diseases and which are not.
uni_clin_sig <- unique(full_SNP_Annotations_GWASc_Ensembl$EnsVar_clinical_significance)
uni_clin_sig # 139 results, not going to be useful for capturing all disease traits... wonder how to do such a complex filtering

















































