# Clustering Analysis:
#
# Strategy:
#   1. For all SNPs above .4 fst against the meta population create table with values of interest from masterTable ✔
#   2. Generate table for all parent terms first: all pops with binary values for having or not having a given SNP of fst .4 or more✔
#   3. generate an alternate table where the Fst value for each pop wrt each SNP is the value (so no binary) ✔
#   4. preprocess in any way necessary for PCA ?? ✔ (R PCA built in func does centering for you and standardization optionally)
#   5. PCA down to two components ✔
#   6. Check the variance captured.. if above 75% then we can move onto graphing, if not we must rework data with some preprocessing or reconsider our approach. ✔ (not quite 75%, but results look good. )
#   7. graph data and see if natural clusters form. Label all points. Possibly use colors to indicate other interesting qualities about each population ✔


# saving important data struc
getwd()
save(above.4_names, file = './workingData/popsVSmetapop_topVariants/above4_names.rds')


uni.4_names <- unique(as.character(purrr::flatten(above.4_names)))

colnames(SNPanno_GWAS_ensembl)

above.4_Table <- SNPanno_GWAS_ensembl[SNPanno_GWAS_ensembl$VariantID %in% uni.4_names, ]
# 3219 obs .. from 477 unique variant IDs

#generate binary table
binPopAbove.4_perVar <- lapply(above.4_names, function(x){as.numeric(uni.4_names %in% x)})
binPopAbove.4_perVar <- dplyr::bind_cols(binPopAbove.4_perVar) # not transposed correctly, but otherwise looks good.

# generate Fst based table
fstPopAbove.4_perVar <- diseaseFst_popsVSall[rownames(diseaseFst_popsVSall) %in% uni.4_names, ]
colnames(fstPopAbove.4_perVar) <- gsub('1000GENOMES:phase_3:','',gsub('1000GENOMES:phase_3:ALL','',gsub('-X-','',colnames(fstPopAbove.4_perVar))))
# looks good! seems much more informative than the binary method.. but may be more noisy too.

# referencing video to learn her method.. taking some notes in my project doc to capture important details.
# Then going to find the right R libraries for PCA and graphing.

?prcomp

pca_1 <- prcomp(fstPopAbove.4_perVar)
#proportion of variance explained
PVE <- (pca_1$sdev**2) / sum(pca_1$sdev**2)
PVE # ok we have [1] 0.3324426994 0.2158026732 0.1119567747 0.0606991383 0.0476742640 0.0366536310 0.0267510183 0.0229971601
# lots of variance explained in the first 2 considering what I saw in my reference video not above 75% though
PVE[1]+PVE[2] # 0.5482454 <- the cumulative variance explained by component 1 and 2

plot(pca_1$x[,1],pca_1$x[,2]) # oh shit... looks like we may need to transpose before doing PCA.. maybe we know nothing yet

# so we are seeing each SNP of our 477 uniques as a point. The problem is labeling them to be honest... there isn't a clear way to say which point belongs to which population unless we use some coarse measure like highest Fst for that SNP gets the label. There would also be 31 colors.

# going to transpose the input and see what we get:
pca_2 <- prcomp(t(fstPopAbove.4_perVar))
PVE_2 <- (pca_2$sdev**2) / sum(pca_2$sdev**2)
PVE[1]+PVE[2] # 0.5482454 ... exactly the same, BUT most of it is contained within the first component now. I am sure Dr. Scheild would like to explain why... I would like to know.
plot(pca_2$x[,1],pca_2$x[,2])


library(ggplot2)

p1 <- ggplot(mtcars, aes(x = mpg, y = wt)) + geom_point()
print(p1)

pop_labels <- colnames(fstPopAbove.4_perVar)
pca_2gg <- ggplot(as.data.frame(pca_2$x), aes(x = pca_2$x[,1], y = pca_2$x[,2]))+geom_point()+labs(pop_labels)+xlab("First Principle Component")+
  ylab("Second Principle Component")+ggtitle("PCA on Variants above 0.4 Fst")+geom_text_repel(aes(label=pop_labels), vjust=1, hjust=1, size=3.5)+
  legend(x = legend_labs$Population_Abbreviation, y = legend_labs$PopAnces_Graph_Labels)
print(pca_2gg)

?geom_text

# want to fix overlapping:
install.packages("ggrepel")
library(ggrepel)
# works beautifully and simply. Nice.

# now to get a helpful key and possibly some color associating sub-populations with super-populations


head(pops)

superPops <- data.frame(c("ACB", "AFR", "AMR", "ASW", "BEB", "CDX", "CEU", "CHB", "CHS", "CLM", "EAS", "ESN", "EUR", "FIN", "GBR", "GIH", "GWD", "IBS", "ITU", "JPT", "KHV", "LWK", "MSL", "MXL", "PEL", "PJL", "PUR", "SAS", "STU", "TSI", "YRI"), c("African", "African", "AdmixedAmerican","African","SouthAsian","EastAsian","European", "EastAsian", "EastAsian", "AdmixedAmerican","EastAsian", "African", "European", "European", "European", "SouthAsian", "African", "European", "SouthAsian", "EastAsian", "EastAsian","African", "African", 'AdmixedAmerican','AdmixedAmerican', "SouthAsian", 'AdmixedAmerican', "SouthAsian","SouthAsian", "European", "African"))
pop_labels

# Citation for super populations
# Gaspar HA, Breen G. Probabilistic ancestry maps: a method to assess and visualize population substructures in genetics. BMC Bioinformatics. 2019 Mar 7;20(1):116. doi: 10.1186/s12859-019-2680-1. PMID: 30845922; PMCID: PMC6407257.
#

legend_labs <- pops[ 2:32,c(1,4)]
legend_labs$Population_Abbreviation <- gsub("1000GENOMES:phase_3:",'',legend_labs$Population_Abbreviation)
legend_labs


# I think I need to make a data.frame with properly associated values to generate a better plot with more color now.

pca_2$x[,1][1:2] # nice they're labeled
length(pca_2$x[,1])

pca_2_DF <- data.frame(pca_2$x[,1], pca_2$x[,2], superPops)
colnames(pca_2_DF) <- c("First Principal Component", "Second Principal Component", "Population Abbrevation", "Super Population")


pca_2DFgg <- ggplot(pca_2_DF, aes(x = `First Principal Component`, y = `Second Principal Component`, colour=`Super Population`))+
  geom_point(shape = "circle", size = 2)+
  labs(pca_2_DF$`Population Abbrevation`)+
  xlab("First Principle Component")+
  ylab("Second Principle Component")+
  ggtitle("PCA on Variants above 0.4 Fst")+
  geom_text_repel(aes(label=pop_labels), vjust=1, hjust=1, size=3.5)+
  scale_color_hue(direction = 1)+
  theme_minimal(base_size = 14)

print(pca_2DFgg)
ggsave(filename = 'PCA_allSubPopsVSmetaPop_withSuperPopulationColor.PNG', scale = 3)
?ggsave

library(esquisse)
esquisser(pca_2_DF) # great for getting the basics correct!

# legendLabs is the abbrevation legend, somehow get it posted within the image of the PCA plot.

# generate sub-parent-term plots:

# 1. extract each parent term as a table
# 2. generate plots


colnames(above.4_Table)
unique(above.4_Table$`Parent term`) # seeing all 17 parent terms

# refiltering according to parental term before splitting off.

diseaseParentTerms <- c('Cancer', 'Other disease', 'Digestive system disorder', 'Cardiovascular disease', 'Neurological disorder', 'Immune system disorder', 'Metabolic disorder')

SNPannoTable_disease <- full_SNP_Anno_withParentalTerms[full_SNP_Anno_withParentalTerms$`Parent term` %in% diseaseParentTerms ,]

above.4_T_diseaseOnly <- above.4_Table[above.4_Table$`Parent term` %in% diseaseParentTerms,]
# down from 3219 obs to 843 obs.

# check that all SNPs we expect to see remain:
uni.4_names[1:10] # bunch of VariantIDs
length(uni.4_names) #477
sum(uni.4_names %in% above.4_T_diseaseOnly$VariantID) #477 Good. We arent losing any.

# split according to parent term, then redo PCA calculation and plotting of that calculation

parentTermsTable <- split(above.4_T_diseaseOnly, above.4_T_diseaseOnly$`Parent term`)

sapply(parentTermsTable, nrow) # number of row per, however names aren't necessarily unique so we will see less unique SNPs to use per disease category
# Cancer    Cardiovascular disease Digestive system disorder    Immune system disorder        Metabolic disorder     Neurological disorder
# 118                        73                       125                       129                        20                       225
# Other disease
# 153

# num unique SNPS per term
sapply(parentTermsTable, \(x){length(unique(x$VariantID))})
# Cancer    Cardiovascular disease Digestive system disorder    Immune system disorder        Metabolic disorder     Neurological disorder
# 66                        51                        75                        74                        18                       162
# Other disease    Total
# 105              477

nrow(fstPopAbove.4_perVar) #477 <- number of unique SNP names in the non-split PCA


# Do first analysis manually, then use as template


cancerFst <- fstPopAbove.4_perVar[rownames(fstPopAbove.4_perVar) %in% parentTermsTable$Cancer$VariantID ,]

pca_Cancer <- prcomp(t(cancerFst))
PVE_Cancer <- (pca_Cancer$sdev**2) / sum(pca_Cancer$sdev**2)
PVE_Cancer[1]+PVE_Cancer[2] # 0.6789227
plot(pca_Cancer$x[,1],pca_Cancer$x[,2])

pca_DF_cancer <- data.frame(pca_Cancer$x[,1], pca_Cancer$x[,2], superPops)
colnames(pca_DF_cancer) <- c("First Principal Component", "Second Principal Component", "Population Abbrevation", "Super Population")

pca_CancerDFgg <- ggplot(pca_DF_cancer, aes(x = `First Principal Component`, y = `Second Principal Component`, colour=`Super Population`))+
  geom_point(shape = "circle", size = 2)+
  labs(pca_DF_cancer$`Population Abbrevation`)+
  xlab("First Principle Component")+
  ylab("Second Principle Component")+
  ggtitle("Cancer - PCA on Variants above 0.4 Fst")+
  geom_text_repel(aes(label=pop_labels), vjust=1, hjust=1, size=3.5)+
  scale_color_hue(direction = 1)+
  theme_minimal(base_size = 14)

print(pca_CancerDFgg)


####################

cardiovascular_diseaseFst <- fstPopAbove.4_perVar[rownames(fstPopAbove.4_perVar) %in% parentTermsTable$`Cardiovascular disease`$VariantID ,]

pca_cardiovascular_disease <- prcomp(t(cardiovascular_diseaseFst))
PVE_cardiovascular_disease <- (pca_cardiovascular_disease$sdev**2) / sum(pca_cardiovascular_disease$sdev**2)
PVE_cardiovascular_disease[1]+PVE_cardiovascular_disease[2] # 0.7501739
plot(pca_cardiovascular_disease$x[,1],pca_cardiovascular_disease$x[,2])

pca_DF_cardiovascular_disease <- data.frame(pca_cardiovascular_disease$x[,1], pca_cardiovascular_disease$x[,2], superPops)
colnames(pca_DF_cardiovascular_disease) <- c("First Principal Component", "Second Principal Component", "Population Abbrevation", "Super Population")

pca_cardiovascular_diseaseDFgg <- ggplot(pca_DF_cardiovascular_disease, aes(x = `First Principal Component`, y = `Second Principal Component`, colour=`Super Population`))+
  geom_point(shape = "circle", size = 2)+
  labs(pca_DF_cardiovascular_disease$`Population Abbrevation`)+
  xlab("First Principle Component")+
  ylab("Second Principle Component")+
  ggtitle("Cardiovascular Disease - PCA on Variants above 0.4 Fst")+
  geom_text_repel(aes(label=pop_labels), vjust=1, hjust=1, size=3.5)+
  scale_color_hue(direction = 1)+
  theme_minimal(base_size = 14)

print(pca_cardiovascular_diseaseDFgg)


####################


digestive_disorderFst <- fstPopAbove.4_perVar[rownames(fstPopAbove.4_perVar) %in% parentTermsTable$`Digestive system disorder`$VariantID ,]

pca_digestive_disorder <- prcomp(t(digestive_disorderFst))
PVE_digestive_disorder <- (pca_digestive_disorder$sdev**2) / sum(pca_digestive_disorder$sdev**2)
PVE_digestive_disorder[1]+PVE_digestive_disorder[2] # 0.7359454
plot(pca_digestive_disorder$x[,1],pca_digestive_disorder$x[,2])

pca_DF_digestive_disorder <- data.frame(pca_digestive_disorder$x[,1], pca_digestive_disorder$x[,2], superPops)
colnames(pca_DF_digestive_disorder) <- c("First Principal Component", "Second Principal Component", "Population Abbrevation", "Super Population")

pca_digestive_disorderDFgg <- ggplot(pca_DF_digestive_disorder, aes(x = `First Principal Component`, y = `Second Principal Component`, colour=`Super Population`))+
  geom_point(shape = "circle", size = 2)+
  labs(pca_DF_digestive_disorder$`Population Abbrevation`)+
  xlab("First Principle Component")+
  ylab("Second Principle Component")+
  ggtitle("Digestive System Disorder - PCA on Variants above 0.4 Fst")+
  geom_text_repel(aes(label=pop_labels), vjust=1, hjust=1, size=3.5)+
  scale_color_hue(direction = 1)+
  theme_minimal(base_size = 14)

print(pca_digestive_disorderDFgg)


####################


immune_disorderFst <- fstPopAbove.4_perVar[rownames(fstPopAbove.4_perVar) %in% parentTermsTable$`Immune system disorder`$VariantID ,]

pca_immune_disorder <- prcomp(t(immune_disorderFst))
PVE_immune_disorder <- (pca_immune_disorder$sdev**2) / sum(pca_immune_disorder$sdev**2)
PVE_immune_disorder[1]+PVE_immune_disorder[2] # 0.69826
plot(pca_immune_disorder$x[,1],pca_immune_disorder$x[,2])

pca_DF_immune_disorder <- data.frame(pca_immune_disorder$x[,1], pca_immune_disorder$x[,2], superPops)
colnames(pca_DF_immune_disorder) <- c("First Principal Component", "Second Principal Component", "Population Abbrevation", "Super Population")

pca_immune_disorderDFgg <- ggplot(pca_DF_immune_disorder, aes(x = `First Principal Component`, y = `Second Principal Component`, colour=`Super Population`))+
  geom_point(shape = "circle", size = 2)+
  labs(pca_DF_immune_disorder$`Population Abbrevation`)+
  xlab("First Principle Component")+
  ylab("Second Principle Component")+
  ggtitle("Immune System Disorder - PCA on Variants above 0.4 Fst")+
  geom_text_repel(aes(label=pop_labels), vjust=1, hjust=1, size=3.5)+
  scale_color_hue(direction = 1)+
  theme_minimal(base_size = 14)

print(pca_immune_disorderDFgg)


####################


metabolic_disorderFst <- fstPopAbove.4_perVar[rownames(fstPopAbove.4_perVar) %in% parentTermsTable$`Metabolic disorder`$VariantID ,]

pca_metabolic_disorder <- prcomp(t(metabolic_disorderFst))
PVE_metabolic_disorder <- (pca_metabolic_disorder$sdev**2) / sum(pca_metabolic_disorder$sdev**2)
PVE_metabolic_disorder[1]+PVE_metabolic_disorder[2] # 0.7501739
plot(pca_metabolic_disorder$x[,1],pca_metabolic_disorder$x[,2])

pca_DF_metabolic_disorder <- data.frame(pca_metabolic_disorder$x[,1], pca_metabolic_disorder$x[,2], superPops)
colnames(pca_DF_metabolic_disorder) <- c("First Principal Component", "Second Principal Component", "Population Abbrevation", "Super Population")

pca_metabolic_disorderDFgg <- ggplot(pca_DF_metabolic_disorder, aes(x = `First Principal Component`, y = `Second Principal Component`, colour=`Super Population`))+
  geom_point(shape = "circle", size = 2)+
  labs(pca_DF_metabolic_disorder$`Population Abbrevation`)+
  xlab("First Principle Component")+
  ylab("Second Principle Component")+
  ggtitle("Metabolic Disorder - PCA on Variants above 0.4 Fst")+
  geom_text_repel(aes(label=pop_labels), vjust=1, hjust=1, size=3.5)+
  scale_color_hue(direction = 1)+
  theme_minimal(base_size = 14)

print(pca_metabolic_disorderDFgg)


####################


neurological_disorderFst <- fstPopAbove.4_perVar[rownames(fstPopAbove.4_perVar) %in% parentTermsTable$`Neurological disorder`$VariantID ,]

pca_neurological_disorder <- prcomp(t(neurological_disorderFst))
PVE_neurological_disorder <- (pca_neurological_disorder$sdev**2) / sum(pca_neurological_disorder$sdev**2)
PVE_neurological_disorder[1]+PVE_neurological_disorder[2] # 0.7752358
plot(pca_neurological_disorder$x[,1],pca_neurological_disorder$x[,2])

pca_DF_neurological_disorder <- data.frame(pca_neurological_disorder$x[,1], pca_neurological_disorder$x[,2], superPops)
colnames(pca_DF_neurological_disorder) <- c("First Principal Component", "Second Principal Component", "Population Abbrevation", "Super Population")

pca_neurological_disorderDFgg <- ggplot(pca_DF_neurological_disorder, aes(x = `First Principal Component`, y = `Second Principal Component`, colour=`Super Population`))+
  geom_point(shape = "circle", size = 2)+
  labs(pca_DF_neurological_disorder$`Population Abbrevation`)+
  xlab("First Principle Component")+
  ylab("Second Principle Component")+
  ggtitle("Neurological Disorder - PCA on Variants above 0.4 Fst")+
  geom_text_repel(aes(label=pop_labels), vjust=1, hjust=1, size=3.5)+
  scale_color_hue(direction = 1)+
  theme_minimal(base_size = 14)

print(pca_neurological_disorderDFgg)


####################


other_diseaseFst <- fstPopAbove.4_perVar[rownames(fstPopAbove.4_perVar) %in% parentTermsTable$`Other disease`$VariantID ,]

pca_other_disease <- prcomp(t(other_diseaseFst))
PVE_other_disease <- (pca_other_disease$sdev**2) / sum(pca_other_disease$sdev**2)
PVE_other_disease[1]+PVE_other_disease[2] # 0.718618
plot(pca_other_disease$x[,1],pca_other_disease$x[,2])

pca_DF_other_disease <- data.frame(pca_other_disease$x[,1], pca_other_disease$x[,2], superPops)
colnames(pca_DF_other_disease) <- c("First Principal Component", "Second Principal Component", "Population Abbrevation", "Super Population")

pca_other_diseaseDFgg <- ggplot(pca_DF_other_disease, aes(x = `First Principal Component`, y = `Second Principal Component`, colour=`Super Population`))+
  geom_point(shape = "circle", size = 2)+
  labs(pca_DF_other_disease$`Population Abbrevation`)+
  xlab("First Principle Component")+
  ylab("Second Principle Component")+
  ggtitle("Other Disease - PCA on Variants above 0.4 Fst")+
  geom_text_repel(aes(label=pop_labels), vjust=1, hjust=1, size=3.5)+
  scale_color_hue(direction = 1)+
  theme_minimal(base_size = 14)

print(pca_other_diseaseDFgg)


####################


colnames(diseaseFst_popsVSall) <- gsub('1000GENOMES:phase_3:','',gsub('1000GENOMES:phase_3:ALL','',gsub('-X-','',colnames(diseaseFst_popsVSall))))

pca_all_variants <- prcomp(t(diseaseFst_popsVSall))
PVE_all_variants <- (pca_all_variants$sdev**2) / sum(pca_all_variants$sdev**2)
PVE_all_variants[1]+PVE_all_variants[2] # 0.5696698 .. as expected, less variance is explained when there is more data.. tougher to reduce the data meaningfully
plot(pca_all_variants$x[,1],pca_all_variants$x[,2])

pca_DF_all_variants <- data.frame(pca_all_variants$x[,1], pca_all_variants$x[,2], superPops)
colnames(pca_DF_all_variants) <- c("First Principal Component", "Second Principal Component", "Population Abbrevation", "Super Population")

pca_all_variantsDFgg <- ggplot(pca_DF_all_variants, aes(x = `First Principal Component`, y = `Second Principal Component`, colour=`Super Population`))+
  geom_point(shape = "circle", size = 2)+
  labs(pca_DF_all_variants$`Population Abbrevation`)+
  xlab("First Principle Component")+
  ylab("Second Principle Component")+
  ggtitle("PCA on All Disease Variants")+
  geom_text_repel(aes(label=pop_labels), vjust=1, hjust=1, size=3.5)+
  scale_color_hue(direction = 1)+
  theme_minimal(base_size = 14)

print(pca_all_variantsDFgg)



####################


colnames(allFst_popsVSall) <- gsub('1000GENOMES:phase_3:','',gsub('1000GENOMES:phase_3:ALL','',gsub('-X-','',colnames(allFst_popsVSall))))

pca_all_variants_ALL <- prcomp(t(allFst_popsVSall))
PVE_all_variants_ALL <- (pca_all_variants_ALL$sdev**2) / sum(pca_all_variants_ALL$sdev**2)
PVE_all_variants_ALL[1]+PVE_all_variants_ALL[2] # 0.5797212
plot(pca_all_variants_ALL$x[,1],pca_all_variants_ALL$x[,2])

pca_DF_all_variants_ALL <- data.frame(pca_all_variants_ALL$x[,1], pca_all_variants_ALL$x[,2], superPops)
colnames(pca_DF_all_variants_ALL) <- c("First Principal Component", "Second Principal Component", "Population Abbrevation", "Super Population")

pca_all_variants_ALLDFgg <- ggplot(pca_DF_all_variants_ALL, aes(x = `First Principal Component`, y = `Second Principal Component`, colour=`Super Population`))+
  geom_point(shape = "circle", size = 2)+
  labs(pca_DF_all_variants_ALL$`Population Abbrevation`)+
  xlab("First Principle Component")+
  ylab("Second Principle Component")+
  ggtitle("PCA on All Variants")+
  geom_text_repel(aes(label=pop_labels), vjust=1, hjust=1, size=3.5)+
  scale_color_hue(direction = 1)+
  theme_minimal(base_size = 14)

print(pca_all_variants_ALLDFgg)


#################### gonna cluster the other way around


colnames(allFst_popsVSall) <- gsub('1000GENOMES:phase_3:','',gsub('1000GENOMES:phase_3:ALL','',gsub('-X-','',colnames(allFst_popsVSall))))

pca_all_variants_ALL <- prcomp(allFst_popsVSall)
PVE_all_variants_ALL <- (pca_all_variants_ALL$sdev**2) / sum(pca_all_variants_ALL$sdev**2)
PVE_all_variants_ALL[1]+PVE_all_variants_ALL[2] # 0.5797212
plot(pca_all_variants_ALL$x[,1],pca_all_variants_ALL$x[,2])

pca_DF_all_variants_ALL <- data.frame(pca_all_variants_ALL$x[,1], pca_all_variants_ALL$x[,2], superPops)
colnames(pca_DF_all_variants_ALL) <- c("First Principal Component", "Second Principal Component", "Population Abbrevation", "Super Population")

pca_all_variants_ALLDFgg <- ggplot(pca_DF_all_variants_ALL, aes(x = `First Principal Component`, y = `Second Principal Component`, colour=`Super Population`))+
  geom_point(shape = "circle", size = 2)+
  labs(pca_DF_all_variants_ALL$`Population Abbrevation`)+
  xlab("First Principle Component")+
  ylab("Second Principle Component")+
  ggtitle("PCA on All Variants")+
  geom_text_repel(aes(label=pop_labels), vjust=1, hjust=1, size=3.5)+
  scale_color_hue(direction = 1)+
  theme_minimal(base_size = 14)

print(pca_all_variants_ALLDFgg)







































