# AUTHOR: Jon Tanner Nelson
# LAST UPDATED: 10-4-2022
# CONTENTS: graph_varPopFrequencies(), graph_singlePopAlleleFreq(), filterMerge_masterTable_singlePopTable()
# PURPOSE:  Plotting, graphing, and associated helper functions



#' graph_varPopFrequencies
#'
#' graphing function which takes in a population allele frequency table for a single variant (rs00000000), and outputs a customizable graph of either 'bar' or 'point' geometry-styles.
#'
#' See 'Arguments' for details on graph customization / options. It is strongly suggested to play around with the various optional arguments in order to discover the range of possible outputs when seeking any specific look or style of graph.
#'
#' popFreqTable's are the data frames which are generated by createMT(), get_ensVariants(). The tables are named after the Variant which they contain data for. They contain 'allele', 'population', 'allele_count' and 'frequency' variables
#'
#' @param popFreqTable Data in. - data.table of population data for a single variant (SNP), popFreqTables are output from get_ensVariants() which is called by createMT().
#' @param numPopulations is the number of populations shown on the graph. If it is not adjusted by the user, the graph is unsorted by default, if the user changes the value of numPopulations the data graphed will be sorted descending by default, but can be adjusted to ascending by using the 'ascending' option.
#' @param graph_style can be either 'bar' or 'point' .. changes the geometry of the graph, the two different styles convey some different information as well. (Try them both!)
#' @param ascending if TRUE the graphs will display in acsending order rather than descending order (which is the default behavior)
#' @param orderless if TRUE, will overide 'ascending' and enforce orderlessness upon the x axis elements
#' @param export determines whether or not to export the graph generated
#' @param file_name allows the user to specify the name of the file output if 'export' is set to TRUE
#' @param graph_scale determines the scale of the graph for exported files, larger values generate higher resolution graphs
#' @param sqrtYscale When TRUE, very small frequency values are made visible by elongating the lower range of the y-axis. sqrtYscale option will over-ride the 'yUpperLim' option
#' @param yUpperLim sets the y-axis to a custom value, by default the graphing function ggplot() automatically determines the y-axis limit within this function call.
#'
#' @return ggplot2 graph object
#'
#' @examples
#' # using package data object 'testMasterList' for this example
#'    populationFrequencyTable_1 <- testMasterList[[2]][[1]]
#'
#' #using all default options for the graph, if you wish to learn more about what graphs this function can produce, check the vignettes on https://github.com/J-T-Nelson/GWASpops.pheno2geno
#'    graph_varPopFrequencies(populationFrequencyTable_1)
#'
#'
#' @import ggplot2
#' @importFrom dplyr filter
#' @importFrom dplyr select
#' @importFrom dplyr arrange
#' @importFrom dplyr desc
#' @importFrom dplyr distinct
#' @importFrom dplyr inner_join
#' @importFrom stats reorder
#'
#' @export
graph_varPopFrequencies <- function(popFreqTable,
                                        numPopulations = nrow(preparedData),
                                        graph_style = 'bar',
                                        orderless = FALSE,
                                        ascending = FALSE,
                                        export = FALSE,
                                        file_name = 'PopFrequencyGraph',
                                        graph_scale = 3,
                                        sqrtYscale = FALSE,
                                        yUpperLim = NA ) {


  preparedData <- popFreqTable %>% filter(popFreqTable$allele != attr(popFreqTable, 'Ancestral_Allele'))

  if(numPopulations > nrow(preparedData)){ #dealing with invalid user input for numPopulations
    warning('numPopulations specified is greater than total populations in table supplied. terminating function', '\nnumPopulations in table supplied = ', nrow(preparedData))
    stop()
  }

  if(numPopulations < nrow(preparedData)){
    preparedData <- arrange(preparedData, desc(frequency))
  }

  pointSize <- expression()
  if(graph_style == 'point'){
    graphGeometry <- expression(geom_point(shape = 'circle', aes(colour = allele, size = allele_count,  eval(ordering))))
    pointSize <- expression(scale_size_binned(n.breaks = 10))
  } else{ graphGeometry <- expression(geom_col(aes(fill = allele, eval(ordering)), width = .8, alpha = .8))}

  ordering <- expression(population, population)
  if(!orderless){
    ordering <-  expression(reorder(population, -frequency))
    if(ascending){
      ordering <- expression(reorder(population, frequency))
    }
  }

  sqYscale <- expression() # blank expression() will return NULL which don't affect the graph, upon activation these can fill with code-expressions which will modify the graph
  yUp_lim <- expression()

  if(sqrtYscale){
    sqYscale <- expression(scale_y_sqrt())
    if(!is.na(yUpperLim)){
      cat('WARNING: both yUpperLim and sqrtYscale are TRUE, yUpperLim will override sqrtYscale\n')
    }
  }
  if(!is.na(yUpperLim)){
    yUp_lim <- expression(ylim(0, yUpperLim))
    cat('WARNING: If y-axis limit value is below max values being graphed, the values will not be included in the graph!!\n')
  }


  graph <- preparedData[1:numPopulations, ] %>% ggplot(aes(population, frequency)) +
    eval(graphGeometry) +
    eval(pointSize) +
    theme_bw(base_size = 15) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    labs(title = paste0('Variant ID: ', attr(popFreqTable, 'VariantID')),
         subtitle = paste0('Ancestral Allele: ', attr(popFreqTable, 'Ancestral_Allele')),
         size = 'Individuals Sampled') +
    xlab('Population Ancestry')+ylab('Allele Frequency') +
    eval(sqYscale) + eval(yUp_lim)


  if(export){
    cat('graph_scale is 3 by default, ratio is 13:8\n')
    file_name <- paste0(file_name, '.png')
    ggsave(filename = file_name, scale = graph_scale)
  }
  return(graph)

}


#' graph_singlePopAlleleFreq
#'
#' Graphing function which takes in a master table and a single-population allele frequency table and outputs a customizable graph which is either 'point' or 'bar' geometry / style.
#'
#'
#' See 'Arguments' for details on graph customization / options. It is strongly suggested to play around with the various optional arguments in order to discover the range of possible outputs when seeking any specific look or style of graph.
#'
#' @param masterTable Data in. - Table of SNP data from GWAS catalog and Ensembl-variants endpoint. Output by get_ensVariants(). Contains 40 variables.
#' @param singlePopTable Data in - A population allele / SNP frequency table which is produced by singlePopTransform() function. Contains 'VariantID' 'allele_count' 'frequency' 'allele' variables, as well as a population attribute
#' @param numVariants specifys the number of variants to be displayed by the graph.
#' @param graph_style can be either 'bar' / 'col' or 'point' .. changes the geometry of the graph, the two different styles convey some different information as well. (Try them both!) (bar and col are synonymous and produce the same graph)
#' @param facet_graph when TRUE, graph will be faceted according to 'Reported_trait'.
#' @param orderless when TRUE, the graph will be unordered and thus bars / points will be randomly distributed instead of sorted
#' @param ascending if TRUE the graphs will display in acsending order rather than descending order (which is the default behavior)
#' @param show_genes displays gene names of each variants above their respective data point on the graph.
#' @param pValAsSize changes point sizes according to p-value of each SNP, ONLY FUNCTIONS for 'point' graph_style
#' @param export determines whether or not to export the graph generated
#' @param file_name allows the user to specify the name of the file output if 'export' is set to TRUE
#' @param graph_scale determines the scale of the graph for exported files, larger values generate higher resolution graphs
#' @param sqrtYscale When TRUE, very small frequency values are made visible by elongating the lower range of the y-axis. sqrtYscale option will over-ride the 'yUpperLim' option
#' @param yUpperLim sets the y-axis to a custom value, by default graphing function ggplot() automatically determines the y-axis limit; yUpperLim will override sqrtYscale if both are active
#'
#' @return ggplot2 graph object
#'
#' @examples
#' # using package data object 'testMasterList' for this example
#'    testMasterTable <- testMasterList[[1]]
#'    singlePopulation <- testMasterList[[3]][[1]]
#'
#' # graphing with default options, see vignettes at https://github.com/J-T-Nelson/GWASpops.pheno2geno for more details on graphs which this function can produce
#'    graph_singlePopAlleleFreq(testMasterTable, singlePopulation)
#'
#' @export
graph_singlePopAlleleFreq <- function(
                                      masterTable,
                                      singlePopTable,
                                      numVariants = nrow(preparedData),
                                      graph_style = 'point',
                                      facet_graph = FALSE,
                                      orderless = FALSE,
                                      ascending = FALSE,
                                      show_genes = FALSE,
                                      pValAsSize = FALSE,
                                      export = FALSE,
                                      file_name = 'SinglePopFreqGraph',
                                      graph_scale = 3,
                                      sqrtYscale = FALSE,
                                      yUpperLim = NA){

  #Below function fully fixes data for this graph type; merges necessary rows and cols of masterTable with singlePopTable
  preparedData <- filterMerge_masterTable_singlePopTable(masterTable, singlePopTable)

  if(numVariants > nrow(preparedData)){ #dealing with invalid user input for numVariants
    warning('numVariants specified is greater than total variants in table supplied. terminating function', '\nnumVariants in table supplied = ', nrow(preparedData))
    stop()
  }

  if(numVariants < nrow(preparedData)){ #if less than max rows are selected, the highest frequency rows will be preferentially graphed
    preparedData <- arrange(preparedData, desc(frequency))
  }

  ordering <- expression(VariantID)
  if(!orderless){# if orderless is TRUE, skip ordering statements, place an empty expression instead
    ordering <-  expression(reorder(VariantID, -frequency))
    if(ascending){
      ordering <- expression(reorder(VariantID, frequency))
    }
  }

  sqYscale <- expression() # blank expression() will return NULL which don't affect the graph, upon activation these can fill with code-expressions which will modify the graph
  yUp_lim <- expression()
  if(sqrtYscale){
    sqYscale <- expression(scale_y_sqrt())
  }
  if(!is.na(yUpperLim)){
    yUp_lim <- expression(ylim(0, yUpperLim))
    cat('WARNING: If y-axis limit value is below max values being graphed, the values will not be included in the graph!!\n')
  }
  if(!is.na(yUpperLim) && sqrtYscale){
    cat('WARNING: if yUpperLim and sqrtYscale are both activated, yUpperLim will take priority and override sqrtYscale')
  }

  faceting <- expression()
  if(facet_graph){
    faceting <- expression(facet_wrap(vars(`Reported trait`)))
  }

  gene_display <- expression()
  if(show_genes){
    gene_display <- expression(geom_text(aes(x = VariantID, y = frequency, label = `Mapped gene`, angle = 20,
                                             vjust = -1., hjust = 0), size = 3.5 , color = "#0a0a0a"))
  }

  if(pValAsSize){
    pVal <- -log10(GWAScat_sciConvert(preparedData$`P-value`))
    preparedData$`P-value` <- pVal
  }

  # below if/else statements are unavoidable due to option incompatibility in the function. They simply permit different critical options to function correctly
  if(graph_style == 'bar' || graph_style == 'col'){

    if(orderless){

      graph <- preparedData[1:numVariants, ] %>%
        ggplot() +
        geom_col(aes(x = VariantID, y = frequency, group = `Reported trait`, fill = EnsVar_most_severe_consequence)) +
        eval(gene_display) +
        theme_light(base_size = 14) + theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
        ylab('Allele Frequency') + xlab('Variant ID') +
        labs(title = paste0('Population: ', attr(singlePopTable, 'population')), color = 'Variant Consequence') +
        eval(sqYscale) + eval(yUp_lim) + eval(faceting)

    } else {

        graph <- preparedData[1:numVariants, ] %>%
          ggplot(aes(x = VariantID, y = frequency )) +
          geom_col(aes(group = `Reported trait`, fill = EnsVar_most_severe_consequence, eval(ordering))) +
          theme_light(base_size = 14) + theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
          ylab('Allele Frequency') + xlab('Variant ID') +
          labs(title = paste0('Population: ', attr(singlePopTable, 'population')), color = 'Variant Consequence') +
          eval(gene_display) + eval(sqYscale) + eval(yUp_lim) + eval(faceting)
      }
  } else {

    if(pValAsSize){ # adjusting size according to P-value ... this option was not code-able in a modular fashion thus the repetitious code
                   # x = eval(ordering) because this is how you tell ggplot to order the x axis a specified way
      graph <- preparedData[1:numVariants, ] %>%
        ggplot(aes(x = eval(ordering), y = frequency, size = `P-value`)) +
        geom_point(aes(x = eval(ordering), y = frequency, group = `Reported trait`, color = EnsVar_most_severe_consequence), shape = "circle") +
        theme_light(base_size = 14) + theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
        ylab('Allele Frequency') + xlab('Variant ID') +
        labs(title = paste0('Population: ', attr(singlePopTable, 'population')), subtitle = "Point Area = -log(P Value)", color = 'Variant Consequence') +
        eval(gene_display) + eval(sqYscale) + eval(yUp_lim) + eval(faceting) +
        guides(size = "none")
    } else {

      graph <- preparedData[1:numVariants, ] %>%
        ggplot(aes(x = eval(ordering), y = frequency)) +
        geom_point(aes(x = eval(ordering), y = frequency, group = `Reported trait`, color = EnsVar_most_severe_consequence), shape = "circle", size = 3) +
        theme_light(base_size = 14) + theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
        ylab('Allele Frequency') + xlab('Variant ID') +
        labs(title = paste0('Population: ', attr(singlePopTable, 'population')), color = 'Variant Consequence') +
        eval(gene_display) + eval(sqYscale) + eval(yUp_lim) + eval(faceting) +
        guides(size = "none")
    }
  }

  if(export){ #exporting code block
    cat('graph_scale is 3 by default, ratio is 13:8\n')
    file_name <- paste0(file_name, '.png')
    ggsave(filename = file_name, scale = graph_scale )
  }
  return(graph)
}


# HELPER FUNCS: DATA TRANSFORMATION FOR GRAPHING  ---------------------------------------

#' filterMerge_masterTable_singlePopTable
#'
#' Merges useful rows and cols of masterTable with singlePopTable. helper func for graph_singlePopAlleleFreq(). Transforms data for graphing purposes.
#'
#' Takes data from masterTable (GWAS and Ensembl Variants endpoint data) and merges it with the sinlge-population allele frequency data. Thus opening up graphing options by combining variables from 2 distinct tables into 1 new table.
#'
#' @param masterTable object produced by createMT(..., varAnnotations = TRUE, population_data = TRUE). Contains ~40 variables to start, only useful variables are merged into the new table used for graphing.
#' @param singlePopTable a population allele / SNP frequency table which is produced by singlePopTransform() function. Contains 'VariantID' 'allele_count' 'frequency' 'allele' variables, as well as a population attribute
#'
#' @return data.frame / tibble
#'
#' @example N/A
#'
#' @noRd
filterMerge_masterTable_singlePopTable <- function(masterTable, singlePopTable) {

  #filtering out undesired cols and rows of master table; Also removing duplicate rows which are created from removal of distinguishing data
  MT <- select(masterTable, VariantID, 'Reported trait', 'P-value',
               'P-value annotation', 'Discovery sample number and ancestry', OR,
               CI, 'Mapped gene', Location, 'Study accession',
               EnsVar_most_severe_consequence) %>%
    distinct(.keep_all = FALSE)


  MT <- filter(MT, MT$EnsVar_most_severe_consequence %in% c('missense_variant',
                                                            'regulatory_region_variant',
                                                            'TF_binding_site_variant',
                                                            'intergenic_variant'))

  #merging the two tables .. inner join is used in order to not include rows from either table which have no match in the other (and to keep all cols from each)
  popFreq_MTmerger <- inner_join(singlePopTable, MT)

  return(popFreq_MTmerger)
}
































