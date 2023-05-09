# This work has been moved from "disease_prevalence_vs_snp_frequency.R"
# my failed graph work, attempting to get order to work:


?reorder


# going to try working things more manually. do all preprocessing in one func.. then graphing manually first.



preprocess_prev_alleleFreq_graph <- function(prev_DF, mean_freq_list){
  # pass single prev_df in and the whole mean_freq_list

  current_disease_name <- unique(prev_DF$cause_name)

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
  return(list(prev_DF, super_pops))
}


prev_meanFreq_plot <- function(prev_DF, super_pops){
  ## fill with function graphing code
}



graph_prev_vs_meanFreq <- function(prev_DF, mean_freq_list){
  preprocessed_data <- preprocess_prev_alleleFreq_graph(prev_DF, mean_freq_list)
  prev_DF <- preprocessed_data[[1]]
  super_pops <- preprocessed_data[[2]]
  prev_meanFreq_plot(prev_DF, super_pops)
}



# use preprocessing code then manual plotting funcs to get the plots functioning


processed_data <- preprocess_prev_alleleFreq_graph(prev_DF = prev_DFs$Asthma, mean_freq_list = mean_frequency_DF_list)
full_plot_data <- data.frame(processed_data[[1]], allele_freq_labels = processed_data[[2]]$allele_freq_labels,
                             aggregate_allele_freqs = processed_data[[2]]$aggregate_allele_freqs)
View(full_plot_data)


# Plot 1: Allele frequency labels x Aggregate allele frequencies
plot1 <- ggplot(full_plot_data) +
  geom_bar(aes(x = reorder(allele_freq_labels, super_population), y = aggregate_allele_freqs,
               fill = super_population), stat = "identity") +
  labs(title = paste("Mean Allele Frequency of", "Asthma"), x = "Population", y = "Aggregate Allele Frequency") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
print(plot1)

# Plot 2: Location names x Prevalence
plot2 <- ggplot(full_plot_data) +
  geom_bar(aes(x = reorder(location_name, super_population), y = avg_val, fill = super_population), stat = "identity") +
  labs(title = paste("Prevalence of", "Asthma"), x = "Location", y = "Prevalence") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
print(plot2)



library(esquisse)


full_plot_data_ordered <- full_plot_data[order(full_plot_data$super_population), ]
View(full_plot_data_ordered)



plot1 <- ggplot(full_plot_data_ordered) +
  geom_bar(aes(x = allele_freq_labels, y = aggregate_allele_freqs,
               fill = super_population), stat = "identity") +
  labs(title = paste("Mean Allele Frequency of", "Asthma"), x = "Population", y = "Aggregate Allele Frequency") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
print(plot1)


plot2 <- ggplot(full_plot_data_ordered) +
  geom_bar(aes(x = location_name, y = avg_val, fill = super_population), stat = "identity") +
  labs(title = paste("Prevalence of", "Asthma"), x = "Location", y = "Prevalence") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
print(plot2)




?factor

factor(full_plot_data_ordered$super_population)

str(full_plot_data)

# Define the desired order of super_population levels
super_population_order <- c("SouthAsian", "African", "EastAsian", "AdmixedAmerican", "European")

# Convert the super_population column to a factor and specify the levels
full_plot_data$super_population <- factor(full_plot_data$super_population, levels = super_population_order)
View(full_plot_data)
class(full_plot_data$super_population)

# Plot 1: Allele frequency labels x Aggregate allele frequencies
plot1 <- ggplot(full_plot_data) +
  geom_bar(aes(x = reorder(allele_freq_labels, barOrder), y = aggregate_allele_freqs,
               fill = super_population), stat = "identity") +
  labs(title = paste("Mean Allele Frequency of", "Asthma"), x = "Population", y = "Aggregate Allele Frequency") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
print(plot1)

# Plot 2: Location names x Prevalence
plot2 <- ggplot(full_plot_data) +
  geom_bar(aes(x = reorder(location_name, barOrder), y = avg_val, fill = super_population), stat = "identity") +
  labs(title = paste("Prevalence of", "Asthma"), x = "Location", y = "Prevalence") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
print(plot2)




rownames(full_plot_data) <- 1:nrow(full_plot_data)
barOrder <- as.numeric(rownames(full_plot_data[order(full_plot_data$super_population), ]))
barOrder
full_plot_data$location_name



# Cannot reorder. Have tried many different things. Nothing is doi --------







