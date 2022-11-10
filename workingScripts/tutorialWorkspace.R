
All1kGenomes <- masterList[[3]][[1]]
graph_singlePopAlleleFreq(MT, All1kGenomes)
graph_singlePopAlleleFreq(MT, All1kGenomes, numVariants = 10, graph_style = 'bar', orderless = T)
graph_singlePopAlleleFreq(MT, All1kGenomes, facet_graph = T, pValAsSize = T, orderless = T)

interestingSNP <- masterList[[2]][[10]]
interestingSNP2 <- masterList[[2]][[23]]

graph_varPopFrequencies(interestingSNP)
graph_varPopFrequencies(interestingSNP, ascending = T , yUpperLim = .25)
graph_varPopFrequencies(interestingSNP2, ascending = T , yUpperLim = .25, graph_style = 'point')
?graph_varPopFrequencies

attr(interestingSNP, 'VariantName')
names(interestingSNP)
name(intersetingSNP)
attributes(interestingSNP)

?Populations
?GWAS.asso.study.data
