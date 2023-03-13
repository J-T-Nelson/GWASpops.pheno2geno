# testing hudson Fst estimator to understand if order outputs different values for same pairs


# TESTING THIS ESTIMATOR
HudsonFst <- function(n1, n2, p1, p2){
  return( ((p1-p2)**2 - (p1*(1-p1))/(n1-1) - (p2*(1-p2))/(n2-1)) / (p1*(1-p2) + p2*(1-p1)) )
}


pop1 <- .56
pop2 <- .84
pop3 <- .25
pop4 <- .01

n1 <- 534
n2 <- 233
n3 <- 28
n4 <- 4000

testPops <- c(pop1, pop2, pop3, pop4)
testSampleSizes <- c(n1, n2, n3, n4)

popPairs <- combn(testPops, 2)
popPairs
#       [,1] [,2] [,3] [,4] [,5] [,6]
# [1,] 0.56 0.56 0.56 0.84 0.84 0.25
# [2,] 0.84 0.25 0.01 0.25 0.01 0.01


ssPairs <- combn(testSampleSizes,2)
ssPair
#       [,1] [,2] [,3] [,4] [,5] [,6]
# [1,]  534  534  534  233  233   28
# [2,]  233   28 4000   28 4000 4000

choose(4,2) # 6 ... 6 unique combinations here

# reverse rows to reverse order of entry into function

popPairs_rev <- popPairs[c(2,1), ]
popPairs_rev
      # [,1] [,2] [,3] [,4] [,5] [,6]
# [1,] 0.84 0.25 0.01 0.25 0.01 0.01
# [2,] 0.56 0.56 0.56 0.84 0.84 0.25

ssPairs_rev <- ssPairs[c(2,1), ]


for(i in 1:ncol(popPairs)){
  HudsonFst(ssPairs[1, i], ssPairs[2, i], popPairs[1, i], popPairs[1, i])
}

# didn't know I could use lapply this way... thanks chatGPT
results <- lapply(1:ncol(popPairs), function(i) {
  HudsonFst(ssPairs[1, i], ssPairs[2, i], popPairs[1, i], popPairs[2, i])
})

results_rev <- lapply(1:ncol(popPairs), function(i) {
  HudsonFst(ssPairs_rev[1, i], ssPairs_rev[2, i], popPairs_rev[1, i], popPairs_rev[2, i])
})


captureList <- list()
for(i in 1:length(comparisonNames)){
  captureList[[ comparisonNames[i] ]] <- HudsonFst(ssPairs[1, i], ssPairs[2, i], popPairs[1, i], popPairs[1, i])
}


comparisonNames <- c('pop1_pop2', 'pop1_pop3', 'pop1_pop4', 'pop2_pop3', 'pop2_pop4', 'pop3_pop4')
comparisonNames_rev <- c('pop2_pop1', 'pop3_pop1', 'pop4_pop1', 'pop3_pop2', 'pop4_pop2', 'pop4_pop3')
names(results) <- comparisonNames
names(results_rev) <- comparisonNames_rev

results
# $pop1_pop2
# [1] 0.1684634
#
# $pop1_pop3
# [1] 0.1673458
#
# $pop1_pop4
# [1] 0.5405069
#
# $pop2_pop3
# [1] 0.5083228
#
# $pop2_pop4
# [1] 0.826114
#
# $pop3_pop4
# [1] 0.1986395


results_rev
# $pop2_pop1
# [1] 0.1684634
#
# $pop3_pop1
# [1] 0.1673458
#
# $pop4_pop1
# [1] 0.5405069
#
# $pop3_pop2
# [1] 0.5083228
#
# $pop4_pop2
# [1] 0.826114
#
# $pop4_pop3
# [1] 0.1986395


# we can see clearly here that order doesn't matter. The same result has been produced regardless of the order of processing for each pair, even with substantially different values.
#
# Sanity test complete.








# Comparing lapply to for loop  -------------------------------------------


results_rev <- lapply(1:ncol(popPairs), function(i) {
  HudsonFst(ssPairs_rev[1, i], ssPairs_rev[2, i], popPairs_rev[1, i], popPairs_rev[2, i])
})


captureList <- list()
for(i in 1:length(comparisonNames)){
  captureList[[ comparisonNames[i] ]] <- HudsonFst(ssPairs[1, i], ssPairs[2, i], popPairs[1, i], popPairs[1, i])
}

# lapply is cleaner, but overall they are quite similar in the amount of effort it takes to compose each... so I will use my new powers with care

