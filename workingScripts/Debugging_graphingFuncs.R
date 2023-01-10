
# Not sure if bugs exist, but this will test for them firstly  ------------
# secondly I will debug any issues that come up in this process!
#



# I think I will create my own testing functions.. ------------------------
  # and then after work is done, I will look at testThat package and see if it offers any advantages on my own testing framework
  # data in is always populations data..


# the goal of the testing functions is to attempt all graphing func option-permutations s.t. errors that may occur for users are made obvious
#
test_varPopFreq <- function(dataIn, tableNum = 1) {

  dataIn <- dataIn[[2]][[tableNum]]

  #graph_varPopFrequencies(dataIn)
  graph_varPopFrequencies(dataIn, graph_style = 'bar', orderless = F, ascending = F, graph_scale = 3, sqrtYscale = F, yUpperLim = NA) # DEFAULT
  graph_varPopFrequencies(dataIn, graph_style = 'bar', orderless = T, ascending = F, graph_scale = 3, sqrtYscale = F, yUpperLim = NA)
  graph_varPopFrequencies(dataIn, graph_style = 'bar', orderless = T, ascending = T, graph_scale = 3, sqrtYscale = F, yUpperLim = NA)
  graph_varPopFrequencies(dataIn, graph_style = 'bar', orderless = T, ascending = T, graph_scale = 4, sqrtYscale = T, yUpperLim = NA)
  graph_varPopFrequencies(dataIn, graph_style = 'bar', orderless = T, ascending = T, graph_scale = 4, sqrtYscale = T, yUpperLim = .5)
  graph_varPopFrequencies(dataIn, graph_style = 'bar', orderless = T, ascending = T, graph_scale = 4, sqrtYscale = F, yUpperLim = .5)
  graph_varPopFrequencies(dataIn, graph_style = 'bar', orderless = T, ascending = F, graph_scale = 4, sqrtYscale = F, yUpperLim = .5)
  graph_varPopFrequencies(dataIn, graph_style = 'bar', orderless = F, ascending = F, graph_scale = 4, sqrtYscale = F, yUpperLim = .5)
  graph_varPopFrequencies(dataIn, graph_style = 'bar', orderless = F, ascending = T, graph_scale = 4, sqrtYscale = F, yUpperLim = .5)
  graph_varPopFrequencies(dataIn, numPopulations = 20, graph_style = 'bar', orderless = F, ascending = T, graph_scale = 4, sqrtYscale = F, yUpperLim = .5)
  graph_varPopFrequencies(dataIn, numPopulations = 20, graph_style = 'bar', orderless = F, ascending = T, graph_scale = 4, sqrtYscale = F, yUpperLim = .5, export = T, file_name = 'barVarPop_fullOptions_testExport')
  graph_varPopFrequencies(dataIn, graph_style = 'bar', orderless = F, ascending = T, graph_scale = 4, sqrtYscale = F, yUpperLim = 1, export = T, file_name = 'barVarPop_fullOptions_testExport2')


  graph_varPopFrequencies(dataIn, graph_style = 'point', orderless = F, ascending = F, graph_scale = 3, sqrtYscale = F, yUpperLim = NA) # DEFAULT
  graph_varPopFrequencies(dataIn, graph_style = 'point', orderless = T, ascending = F, graph_scale = 3, sqrtYscale = F, yUpperLim = NA)
  graph_varPopFrequencies(dataIn, graph_style = 'point', orderless = T, ascending = T, graph_scale = 3, sqrtYscale = F, yUpperLim = NA)
  graph_varPopFrequencies(dataIn, graph_style = 'point', orderless = T, ascending = T, graph_scale = 4, sqrtYscale = T, yUpperLim = NA)
  graph_varPopFrequencies(dataIn, graph_style = 'point', orderless = T, ascending = T, graph_scale = 4, sqrtYscale = T, yUpperLim = .5)
  graph_varPopFrequencies(dataIn, graph_style = 'point', orderless = T, ascending = T, graph_scale = 4, sqrtYscale = F, yUpperLim = .5)
  graph_varPopFrequencies(dataIn, graph_style = 'point', orderless = T, ascending = F, graph_scale = 4, sqrtYscale = F, yUpperLim = .5)
  graph_varPopFrequencies(dataIn, graph_style = 'point', orderless = F, ascending = F, graph_scale = 4, sqrtYscale = F, yUpperLim = .5)
  graph_varPopFrequencies(dataIn, graph_style = 'point', orderless = F, ascending = T, graph_scale = 4, sqrtYscale = F, yUpperLim = .5)
  graph_varPopFrequencies(dataIn, numPopulations = 20, graph_style = 'point', orderless = F, ascending = T, graph_scale = 4, sqrtYscale = F, yUpperLim = .5)
  graph_varPopFrequencies(dataIn, numPopulations = 20, graph_style = 'point', orderless = F, ascending = T, graph_scale = 4, sqrtYscale = F, yUpperLim = .5, export = T, file_name = 'pointVarPop_fullOptions_testExport')
  graph_varPopFrequencies(dataIn, graph_style = 'point', orderless = F, ascending = T, graph_scale = 4, sqrtYscale = F, yUpperLim = 1, export = T, file_name = 'pointVarPop_fullOptions_testExport2')
}




test_singlePopAlleleFreq <- function(dataIn, popTableNum = 1){

  masterTable <- dataIn[[1]]
  SinglePopTable <- dataIn[[3]][[popTableNum]]

  #graph_singlePopAlleleFreq(dataIn)
  graph_singlePopAlleleFreq(masterTable, SinglePopTable, graph_style = 'point', facet_graph = F, orderless = F, ascending = F, show_genes =  F, pValAsSize = F, graph_scale = 3, sqrtYscale = F, yUpperLim = NA) #default
  graph_singlePopAlleleFreq(masterTable, SinglePopTable, graph_style = 'point', facet_graph = T, orderless = F, ascending = F, show_genes =  F, pValAsSize = F, graph_scale = 3, sqrtYscale = F, yUpperLim = NA)
  graph_singlePopAlleleFreq(masterTable, SinglePopTable, graph_style = 'point', facet_graph = T, orderless = T, ascending = F, show_genes =  F, pValAsSize = F, graph_scale = 3, sqrtYscale = F, yUpperLim = NA)
  graph_singlePopAlleleFreq(masterTable, SinglePopTable, graph_style = 'point', facet_graph = T, orderless = T, ascending = T, show_genes =  F, pValAsSize = F, graph_scale = 3, sqrtYscale = F, yUpperLim = NA)
  graph_singlePopAlleleFreq(masterTable, SinglePopTable, graph_style = 'point', facet_graph = T, orderless = T, ascending = T, show_genes =  T, pValAsSize = F, graph_scale = 3, sqrtYscale = F, yUpperLim = NA)
  graph_singlePopAlleleFreq(masterTable, SinglePopTable, graph_style = 'point', facet_graph = T, orderless = T, ascending = T, show_genes =  T, pValAsSize = T, graph_scale = 3, sqrtYscale = F, yUpperLim = NA)
  graph_singlePopAlleleFreq(masterTable, SinglePopTable, graph_style = 'point', facet_graph = T, orderless = T, ascending = T, show_genes =  T, pValAsSize = T, graph_scale = 4, sqrtYscale = F, yUpperLim = NA)
  graph_singlePopAlleleFreq(masterTable, SinglePopTable, graph_style = 'point', facet_graph = T, orderless = T, ascending = T, show_genes =  T, pValAsSize = T, graph_scale = 4, sqrtYscale = T, yUpperLim = NA)
  graph_singlePopAlleleFreq(masterTable, SinglePopTable, graph_style = 'point', facet_graph = T, orderless = T, ascending = T, show_genes =  T, pValAsSize = T, graph_scale = 4, sqrtYscale = T, yUpperLim = .75)
  graph_singlePopAlleleFreq(masterTable, SinglePopTable, numVariants = 15, graph_style = 'point', facet_graph = T, orderless = T, ascending = T, show_genes =  T, pValAsSize = T, graph_scale = 4, sqrtYscale = T, yUpperLim = .75, export = T, file_name = "pointSinglePop_full_options_testExport")
  graph_singlePopAlleleFreq(masterTable, SinglePopTable, graph_style = 'point', facet_graph = T, orderless = T, ascending = T, show_genes =  T, pValAsSize = T, graph_scale = 4, sqrtYscale = T, yUpperLim = .75, export = T, file_name = "pointSinglePop_full_options_testExport1")
  graph_singlePopAlleleFreq(masterTable, SinglePopTable, numVariants = 15, graph_style = 'point', facet_graph = T, orderless = T, ascending = T, show_genes =  T, pValAsSize = T, graph_scale = 4, sqrtYscale = T, yUpperLim = 1, export = T, file_name = "pointSinglePop_full_options_testExport2")



  graph_singlePopAlleleFreq(masterTable, SinglePopTable, graph_style = 'bar', facet_graph = F, orderless = F, ascending = F, show_genes =  F, pValAsSize = F, graph_scale = 3, sqrtYscale = F, yUpperLim = NA) #default
  graph_singlePopAlleleFreq(masterTable, SinglePopTable, graph_style = 'bar', facet_graph = T, orderless = F, ascending = F, show_genes =  F, pValAsSize = F, graph_scale = 3, sqrtYscale = F, yUpperLim = NA)
  graph_singlePopAlleleFreq(masterTable, SinglePopTable, graph_style = 'bar', facet_graph = T, orderless = T, ascending = F, show_genes =  F, pValAsSize = F, graph_scale = 3, sqrtYscale = F, yUpperLim = NA)
  graph_singlePopAlleleFreq(masterTable, SinglePopTable, graph_style = 'bar', facet_graph = T, orderless = T, ascending = T, show_genes =  F, pValAsSize = F, graph_scale = 3, sqrtYscale = F, yUpperLim = NA)
  graph_singlePopAlleleFreq(masterTable, SinglePopTable, graph_style = 'bar', facet_graph = T, orderless = T, ascending = T, show_genes =  T, pValAsSize = F, graph_scale = 3, sqrtYscale = F, yUpperLim = NA)
  graph_singlePopAlleleFreq(masterTable, SinglePopTable, graph_style = 'bar', facet_graph = T, orderless = T, ascending = T, show_genes =  T, pValAsSize = T, graph_scale = 3, sqrtYscale = F, yUpperLim = NA)
  graph_singlePopAlleleFreq(masterTable, SinglePopTable, graph_style = 'bar', facet_graph = T, orderless = T, ascending = T, show_genes =  T, pValAsSize = T, graph_scale = 4, sqrtYscale = F, yUpperLim = NA)
  graph_singlePopAlleleFreq(masterTable, SinglePopTable, graph_style = 'bar', facet_graph = T, orderless = T, ascending = T, show_genes =  T, pValAsSize = T, graph_scale = 4, sqrtYscale = T, yUpperLim = NA)
  graph_singlePopAlleleFreq(masterTable, SinglePopTable, graph_style = 'bar', facet_graph = T, orderless = T, ascending = T, show_genes =  T, pValAsSize = T, graph_scale = 4, sqrtYscale = T, yUpperLim = .75)
  graph_singlePopAlleleFreq(masterTable, SinglePopTable, numVariants = 15, graph_style = 'bar', facet_graph = T, orderless = T, ascending = T, show_genes =  T, pValAsSize = T, graph_scale = 4, sqrtYscale = T, yUpperLim = .75, export = T, file_name = "pointSinglePop_full_options_testExport")
  graph_singlePopAlleleFreq(masterTable, SinglePopTable, graph_style = 'bar', facet_graph = T, orderless = T, ascending = T, show_genes =  T, pValAsSize = T, graph_scale = 4, sqrtYscale = T, yUpperLim = .75, export = T, file_name = "barSinglePop_full_options_testExport1")
  graph_singlePopAlleleFreq(masterTable, SinglePopTable, numVariants = 15, graph_style = 'bar', facet_graph = T, orderless = T, ascending = T, show_genes =  T, pValAsSize = T, graph_scale = 4, sqrtYscale = T, yUpperLim = 1, export = T, file_name = "barSinglePop_full_options_testExport2")
}





# Test functions exist for each graph type,--------------------------------
# now we run each data set through the test funcs.

?load_all

?load
?list.files

files <- list.files(path = './data/transformed_data_for_graphing')
for(f in files){
  load(f)
}


# this func should load all data from the transformed data set, then should return all that data in a list s.t. it can easily be iterated through for graph testing.
loadGraphDataBetter <- function() {
  files <- list.files(path = './data/transformed_data_for_graphing')
  emptyL <- vector('list', length(files))

  for(f in 1:length(files)){
    emptyL[f] <- load(files[f]) # not sure if this is valid syntax.
  }
  return(emptyL)
}

#simpler version of loading all data.
loadGraphData <- function() {
  files <- list.files(path = './data/transformed_data_for_graphing')
  for(f in files){
    load(f)
  }
}

# run all premade data through graphing funcs for error detection...
# This current function is quite nested in its structure for efficient coding... but it may not be very easy to debug.. I also have several concerns about whether graphs can be generated in the way I have architected.. we will have to do some testing for the test funcs before actual use of them.
testAllData_graph_varPop <- function(variables) {

  graphData <- loadGraphDataBetter()

  for(d in graphData){
    test_varPopFreq(d)
  }

  for(d in graphData){
    test_singlePopAlleleFreq(d)
  }
}




# hardcoding data loading and testing -------------------------------------

loadTransformedData <- function() {
  load('./data/transformed_data_for_graphing/airPollutionPops.rds')

  load('./data/transformed_data_for_graphing/alcConsumpPops.rds')

  load('./data/transformed_data_for_graphing/bCarcinomaPops.rds')

  load('./data/transformed_data_for_graphing/colorectalCancerPops.rds')

  load('./data/transformed_data_for_graphing/IBFPops.rds')

  load('./data/transformed_data_for_graphing/IntPops.rds')

  load('./data/transformed_data_for_graphing/malabsorptionSyndPops.rds')

  load('./data/transformed_data_for_graphing/lungCancerPops.rds')

  load('./data/transformed_data_for_graphing/neuroticismPops.rds')

  load('./data/transformed_data_for_graphing/prostateCancerPops.rds')

  load('./data/transformed_data_for_graphing/substanceAbusePops.rds')

  load('./data/transformed_data_for_graphing/malabsorptionSyndPops.rds')
}

loadTransformedData() # doesn't fucking work because it all loads within the local environment.... maybe a good github question. I can just make a script and run source() on the script to do the same shit. aight.



test_varPopFreq(AirPollutionAlldata)
test_singlePopAlleleFreq(AirPollutionAlldata)



# First want to just test the graphing functions without any nesting in other structures --------
AirPol_popFreqTable <- AirPollutionAlldata[[2]][[1]]
graph_varPopFrequencies(AirPol_popFreqTable) # good

graph_singlePopAlleleFreq(AirPollutionAlldata[[1]], AirPollutionAlldata[[3]][[1]]) # good

test_varPopFreq(AirPol_popFreqTable)

graph_varPopFrequencies(AirPol_popFreqTable, graph_style = 'bar', orderless = F, ascending = F, graph_scale = 3, sqrtYscale = F, yUpperLim = NA) #works
graph_varPopFrequencies(AirPol_popFreqTable, graph_style = 'bar', orderless = T, ascending = T, graph_scale = 4, sqrtYscale = T, yUpperLim = .5) #works


# Debugging begins... not sure which call is failing...
# ERROR: Error in if (yUpperLim) { : missing value where TRUE/FALSE needed
# ^^^^ the code was just wrong... I think I fixed it.



test_varPopFreq(AirPol_popFreqTable) # success. Everything executed.. many messages, only a single graph is returned.... not sure how to get around that if I want to actually view them all as output... but we can think about that later. For right now, just verifying that the graphs execute is enough. We can worry about aesthetics and visual output later.
# one obvious way to view ouput is to just save the damn things as they get generated.

# TESTING SINGLE POP GRAPHING NOW


test_singlePopAlleleFreq(AirPollutionAlldata) # success! Bunch of warnings and such but nothing broke! ... original data set is validated!



# Write up tests for all other data sets now.  ----------------------------


debug(test_varPopFreq)
debug(test_singlePopAlleleFreq)
undebug(test_varPopFreq)
undebug(test_singlePopAlleleFreq)

graph_varPopFrequencies(alcConsumpAlldata, graph_style = 'bar', orderless = F, ascending = F, graph_scale = 3, sqrtYscale = F, yUpperLim = NA)

test_varPopFreq(alcConsumpAlldata) # success -  data removal (What is happening with this?)
test_singlePopAlleleFreq(alcConsumpAlldata) # ERROR Error in seq_len(n) : argument must be coercible to non-negative integer

test_varPopFreq(bCarcinomaAlldata) # success - data removal...
test_singlePopAlleleFreq(bCarcinomaAlldata) # success ... data removal... some graphs are appearing empty - (need to investigate why, as this should only be happening when graphing options are set incorrectly.. just need to varify good data isn't being misrepresented)

test_varPopFreq(colorectalCancerAlldata) # data removal... apparently empty graphs again..
test_singlePopAlleleFreq(colorectalCancerAlldata) # data removal (DR)

test_varPopFreq(IBFAlldata) # minor DR
test_singlePopAlleleFreq(IBFAlldata) # Error in seq_len(n) : argument must be coercible to non-negative integer

test_varPopFreq(IntAlldata)
test_singlePopAlleleFreq(IntAlldata) # Error in seq_len(n) : argument must be coercible to non-negative integer

test_varPopFreq(lungCancerAlldata)
test_singlePopAlleleFreq(lungCancerAlldata) # no errors

test_varPopFreq(malabsorptionSyndAlldata) # DR
test_singlePopAlleleFreq(malabsorptionSyndAlldata) # DR

test_varPopFreq(neuroticismAlldata)
test_singlePopAlleleFreq(neuroticismAlldata) # Error in seq_len(n) : argument must be coercible to non-negative integer

test_varPopFreq(prostateCancerAlldata)
test_singlePopAlleleFreq(prostateCancerAlldata) # SOME DATA REMOVAL IS HAPPENEING AS A RESULT OF THE yUpperLim parameter!

test_varPopFreq(substanceAbuseAlldata)
test_singlePopAlleleFreq(substanceAbuseAlldata)






# Since I cant do what I need to within a function,  ----------------------
# nor can I easily work out how to pass args to a script, I will be hardcoding the testing now..

debug(graph_varPopFrequencies)
undebug(graph_varPopFrequencies)

debug(graph_singlePopAlleleFreq)
undebug(graph_singlePopAlleleFreq)


malabsorptionSyndAlldata
AirPollutionAlldata
prostateCancerAlldata
colorectalCancerAlldata
substanceAbuseAlldata
lungCancerAlldata
bCarcinomaAlldata
IBFAlldata
alcConsumpAlldata
neuroticismAlldata
IntAlldata

# with the above names I can now just swap things out ad hoc in order to quickly test the different data sets and actually see all the graphs

dataStart <- IntAlldata
tableNum = 2
popTableNum = 2
dataIn <- dataStart[[2]][[tableNum]]

graph_varPopFrequencies(dataIn, graph_style = 'bar', orderless = F, ascending = F, graph_scale = 3, sqrtYscale = F, yUpperLim = NA) # DEFAULT
graph_varPopFrequencies(dataIn, graph_style = 'bar', orderless = T, ascending = F, graph_scale = 3, sqrtYscale = F, yUpperLim = NA)
graph_varPopFrequencies(dataIn, graph_style = 'bar', orderless = T, ascending = T, graph_scale = 3, sqrtYscale = F, yUpperLim = NA)
graph_varPopFrequencies(dataIn, graph_style = 'bar', orderless = T, ascending = T, graph_scale = 4, sqrtYscale = T, yUpperLim = NA)
graph_varPopFrequencies(dataIn, graph_style = 'bar', orderless = T, ascending = T, graph_scale = 4, sqrtYscale = T, yUpperLim = .5)
graph_varPopFrequencies(dataIn, graph_style = 'bar', orderless = T, ascending = T, graph_scale = 4, sqrtYscale = F, yUpperLim = .5)
graph_varPopFrequencies(dataIn, graph_style = 'bar', orderless = T, ascending = F, graph_scale = 4, sqrtYscale = F, yUpperLim = .5)
graph_varPopFrequencies(dataIn, graph_style = 'bar', orderless = F, ascending = F, graph_scale = 4, sqrtYscale = F, yUpperLim = .5)
graph_varPopFrequencies(dataIn, graph_style = 'bar', orderless = F, ascending = T, graph_scale = 4, sqrtYscale = F, yUpperLim = .5)
graph_varPopFrequencies(dataIn, numPopulations = 20, graph_style = 'bar', orderless = F, ascending = T, graph_scale = 4, sqrtYscale = F, yUpperLim = .5)
graph_varPopFrequencies(dataIn, numPopulations = 20, graph_style = 'bar', orderless = F, ascending = T, graph_scale = 4, sqrtYscale = F, yUpperLim = .5, export = T, file_name = 'barVarPop_fullOptions_testExport')
graph_varPopFrequencies(dataIn, graph_style = 'bar', orderless = F, ascending = T, graph_scale = 4, sqrtYscale = F, yUpperLim = .5, export = T, file_name = 'barVarPop_fullOptions_testExport2')

graph_varPopFrequencies(dataIn, graph_style = 'point', orderless = F, ascending = F, graph_scale = 3, sqrtYscale = F, yUpperLim = NA) # DEFAULT
graph_varPopFrequencies(dataIn, graph_style = 'point', orderless = T, ascending = F, graph_scale = 3, sqrtYscale = F, yUpperLim = NA)
graph_varPopFrequencies(dataIn, graph_style = 'point', orderless = T, ascending = T, graph_scale = 3, sqrtYscale = F, yUpperLim = NA)
graph_varPopFrequencies(dataIn, graph_style = 'point', orderless = T, ascending = T, graph_scale = 4, sqrtYscale = T, yUpperLim = NA)
graph_varPopFrequencies(dataIn, graph_style = 'point', orderless = T, ascending = T, graph_scale = 4, sqrtYscale = T, yUpperLim = .5)
graph_varPopFrequencies(dataIn, graph_style = 'point', orderless = T, ascending = T, graph_scale = 4, sqrtYscale = F, yUpperLim = .5)
graph_varPopFrequencies(dataIn, graph_style = 'point', orderless = T, ascending = F, graph_scale = 4, sqrtYscale = F, yUpperLim = .5)
graph_varPopFrequencies(dataIn, graph_style = 'point', orderless = F, ascending = F, graph_scale = 4, sqrtYscale = F, yUpperLim = .5)
graph_varPopFrequencies(dataIn, graph_style = 'point', orderless = F, ascending = T, graph_scale = 4, sqrtYscale = F, yUpperLim = .5)
graph_varPopFrequencies(dataIn, numPopulations = 20, graph_style = 'point', orderless = F, ascending = T, graph_scale = 4, sqrtYscale = F, yUpperLim = .5)
graph_varPopFrequencies(dataIn, numPopulations = 20, graph_style = 'point', orderless = F, ascending = T, graph_scale = 4, sqrtYscale = F, yUpperLim = .5, export = T, file_name = 'pointVarPop_fullOptions_testExport')
graph_varPopFrequencies(dataIn, graph_style = 'point', orderless = F, ascending = T, graph_scale = 4, sqrtYscale = F, yUpperLim = .5, export = T, file_name = 'pointVarPop_fullOptions_testExport2')


masterTable <- dataStart[[1]]
SinglePopTable <- dataStart[[3]][[popTableNum]]

graph_singlePopAlleleFreq(masterTable, SinglePopTable, graph_style = 'point', facet_graph = F, orderless = F, ascending = F, show_genes =  F, pValAsSize = F, graph_scale = 3, sqrtYscale = F, yUpperLim = NA) #default
graph_singlePopAlleleFreq(masterTable, SinglePopTable, graph_style = 'point', facet_graph = T, orderless = F, ascending = F, show_genes =  F, pValAsSize = F, graph_scale = 3, sqrtYscale = F, yUpperLim = NA)
graph_singlePopAlleleFreq(masterTable, SinglePopTable, graph_style = 'point', facet_graph = T, orderless = T, ascending = F, show_genes =  F, pValAsSize = F, graph_scale = 3, sqrtYscale = F, yUpperLim = NA)
graph_singlePopAlleleFreq(masterTable, SinglePopTable, graph_style = 'point', facet_graph = T, orderless = T, ascending = T, show_genes =  F, pValAsSize = F, graph_scale = 3, sqrtYscale = F, yUpperLim = NA)
graph_singlePopAlleleFreq(masterTable, SinglePopTable, graph_style = 'point', facet_graph = T, orderless = T, ascending = T, show_genes =  T, pValAsSize = F, graph_scale = 3, sqrtYscale = F, yUpperLim = NA)
graph_singlePopAlleleFreq(masterTable, SinglePopTable, graph_style = 'point', facet_graph = T, orderless = T, ascending = T, show_genes =  T, pValAsSize = T, graph_scale = 3, sqrtYscale = F, yUpperLim = NA)
graph_singlePopAlleleFreq(masterTable, SinglePopTable, graph_style = 'point', facet_graph = T, orderless = T, ascending = T, show_genes =  T, pValAsSize = T, graph_scale = 4, sqrtYscale = F, yUpperLim = NA)
graph_singlePopAlleleFreq(masterTable, SinglePopTable, graph_style = 'point', facet_graph = T, orderless = T, ascending = T, show_genes =  T, pValAsSize = T, graph_scale = 4, sqrtYscale = T, yUpperLim = NA)
graph_singlePopAlleleFreq(masterTable, SinglePopTable, graph_style = 'point', facet_graph = T, orderless = T, ascending = T, show_genes =  T, pValAsSize = T, graph_scale = 4, sqrtYscale = T, yUpperLim = .75)
graph_singlePopAlleleFreq(masterTable, SinglePopTable, numVariants = 15, graph_style = 'point', facet_graph = T, orderless = T, ascending = T, show_genes =  T, pValAsSize = T, graph_scale = 4, sqrtYscale = T, yUpperLim = .75, export = T, file_name = "pointSinglePop_full_options_testExport")
graph_singlePopAlleleFreq(masterTable, SinglePopTable, numVariants = 15, graph_style = 'point', facet_graph = T, orderless = T, ascending = T, show_genes =  T, pValAsSize = T, graph_scale = 1, sqrtYscale = T, yUpperLim = NA, export = T, file_name = "pointSinglePop_full_options_testExport")
graph_singlePopAlleleFreq(masterTable, SinglePopTable, graph_style = 'point', facet_graph = T, orderless = T, ascending = T, show_genes =  T, pValAsSize = T, graph_scale = 1, sqrtYscale = T, yUpperLim = .75, export = T, file_name = "pointSinglePop_full_options_testExport1")
graph_singlePopAlleleFreq(masterTable, SinglePopTable, numVariants = 15, graph_style = 'point', facet_graph = T, orderless = T, ascending = T, show_genes =  T, pValAsSize = T, graph_scale = 1, sqrtYscale = T, yUpperLim = .75, export = T, file_name = "pointSinglePop_full_options_testExport2")

graph_singlePopAlleleFreq(masterTable, SinglePopTable, graph_style = 'bar', facet_graph = F, orderless = F, ascending = F, show_genes =  F, pValAsSize = F, graph_scale = 3, sqrtYscale = F, yUpperLim = NA) #default
graph_singlePopAlleleFreq(masterTable, SinglePopTable, graph_style = 'bar', facet_graph = T, orderless = F, ascending = F, show_genes =  F, pValAsSize = F, graph_scale = 3, sqrtYscale = F, yUpperLim = NA)
graph_singlePopAlleleFreq(masterTable, SinglePopTable, graph_style = 'bar', facet_graph = T, orderless = T, ascending = F, show_genes =  F, pValAsSize = F, graph_scale = 3, sqrtYscale = F, yUpperLim = NA)
graph_singlePopAlleleFreq(masterTable, SinglePopTable, graph_style = 'bar', facet_graph = T, orderless = T, ascending = T, show_genes =  F, pValAsSize = F, graph_scale = 3, sqrtYscale = F, yUpperLim = NA)
graph_singlePopAlleleFreq(masterTable, SinglePopTable, graph_style = 'bar', facet_graph = T, orderless = T, ascending = T, show_genes =  T, pValAsSize = F, graph_scale = 3, sqrtYscale = F, yUpperLim = NA)
graph_singlePopAlleleFreq(masterTable, SinglePopTable, graph_style = 'bar', facet_graph = T, orderless = T, ascending = T, show_genes =  T, pValAsSize = T, graph_scale = 3, sqrtYscale = F, yUpperLim = NA)
graph_singlePopAlleleFreq(masterTable, SinglePopTable, graph_style = 'bar', facet_graph = T, orderless = T, ascending = T, show_genes =  T, pValAsSize = T, graph_scale = 4, sqrtYscale = F, yUpperLim = NA)
graph_singlePopAlleleFreq(masterTable, SinglePopTable, graph_style = 'bar', facet_graph = T, orderless = T, ascending = T, show_genes =  T, pValAsSize = T, graph_scale = 4, sqrtYscale = T, yUpperLim = NA)
graph_singlePopAlleleFreq(masterTable, SinglePopTable, graph_style = 'bar', facet_graph = T, orderless = T, ascending = T, show_genes =  T, pValAsSize = T, graph_scale = 4, sqrtYscale = T, yUpperLim = .75)

graph_singlePopAlleleFreq(masterTable, SinglePopTable, numVariants = 15, graph_style = 'bar', facet_graph = T, orderless = T, ascending = T, show_genes =  T, pValAsSize = T, graph_scale = 1, sqrtYscale = T, yUpperLim = .75, export = T, file_name = "pointSinglePop_full_options_testExport") # error

graph_singlePopAlleleFreq(masterTable, SinglePopTable, graph_style = 'bar', facet_graph = T, orderless = T, ascending = T, show_genes =  T, pValAsSize = T, graph_scale = 1, sqrtYscale = T, yUpperLim = .75, export = T, file_name = "barSinglePop_full_options_testExport1")

graph_singlePopAlleleFreq(masterTable, SinglePopTable, numVariants = 15, graph_style = 'bar', facet_graph = T, orderless = T, ascending = T, show_genes =  T, pValAsSize = T, graph_scale = 1, sqrtYscale = T, yUpperLim = .75, export = T, file_name = "barSinglePop_full_options_testExport2") # error





































