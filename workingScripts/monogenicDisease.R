# import GWAS data, filter down for reporting on relevant values to monogenic disease questions
#


#"C:\Users\mrtne\Desktop\Kulathinal\GWASc_All_Data_20230114\gwas_catalog_v1.0.2-studies_r2023-01-14.tsv"
#"C:\Users\mrtne\Desktop\Kulathinal\GWASc_All_Data_20230114\gwas_catalog-ancestry_r2023-01-14.tsv"
Associations <- data.table::fread("C:\\Users\\mrtne\\Desktop\\Kulathinal\\GWASc_All_Data_20230114\\gwas_catalog_v1.0.2-associations_e108_r2023-01-14.tsv")
Studies <- data.table::fread("C:\\Users\\mrtne\\Desktop\\Kulathinal\\GWASc_All_Data_20230114\\gwas_catalog_v1.0.2-studies_r2023-01-14.tsv")
Ancestry <- data.table::fread("C:\\Users\\mrtne\\Desktop\\Kulathinal\\GWASc_All_Data_20230114\\gwas_catalog-ancestry_r2023-01-14.tsv")


# DESKTOP FILE LOAD -------------------------------------------------------

"D:\\Kulathinal_files\\GWASc_All_Data_20230114\\"

Associations <- data.table::fread("D:\\Kulathinal_files\\GWASc_All_Data_20230114\\gwas_catalog_v1.0.2-associations_e108_r2023-01-14.tsv")
Studies <- data.table::fread("D:\\Kulathinal_files\\GWASc_All_Data_20230114\\gwas_catalog_v1.0.2-studies_r2023-01-14.tsv")
Ancestry <- data.table::fread("D:\\Kulathinal_files\\GWASc_All_Data_20230114\\gwas_catalog-ancestry_r2023-01-14.tsv")
# vcf files are the standard ... very interoperable... would be useful for the fst stuff.
#

names(Associations[1,])
names(Ancestry[1,])
names(Studies[1,])
colnames(Associations)
#STUDY ACCESSION is a master key across all tables.

# STATS TO REPORT ON:
#  num associations per trait / phenotype
#  num diff traits...
#  traits themselves
#  replication sample size?
#  num individuals - ancestral
#  broad ancestral category?
#  replication sample description?
#
# num associations. num diff traits. MAPPED_TRAIT DISEASE/TRAIT

?merge

# don't think I need to merge. Grabbing from the associations themselves is sufficient
#

# naming monogenic disorders for lookup against associations traits...
# problem is I don't have intelligent searching and thus I think it'd be easy to just use the wrong words wrt the ontology used by GWASc... so manual search and data collection is probably best here.
# The best solution I can think of is to simply turn each separate word into its own pattern then to search through the set of terms within the associations table that match any and to add all of those to a respective list / df as rows then we can report on the summary stats...
#

monoGene1 <- c("Cystic Fibrosis", "Duchenne Muscular Dystrophy", "Friedreich's Ataxia",  "Frataxin", "ALS", "Genetic Amyotrophic Lateral Sclerosis", "Amyotrophic Lateral Sclerosis", "Hemophilia", "Huntington’s Disease", "IRDs", "Inherited Retinal Disorders", "Inherited Retinal Dystrophies", "Rett Syndrome", "Sickle Cell Disease", "Malaria", "Spinal Muscular Atrophy", "SMA")


# LIST 2
#
# Cystic fibrosis.
# Deafness that’s present at birth (congenital).
# Duchenne muscular dystrophy.
# Familial hypercholesterolemia, a type of high cholesterol disease.
# Hemochromatosis (iron overload).
# Neurofibromatosis type 1 (NF1).
# Sickle cell disease.
# Tay-Sachs disease.

# LIST 3
#
# Phenylketonuria (PKU)
# Cystic fibrosis
# Sickle-cell anemia
# Albinism, oculocutaneous, type II
# Huntington's disease
# Myotonic dystrophy type 1
# Hypercholesterolemia, autosomal dominant, type B
# Neurofibromatosis, type 1
# Polycystic kidney disease 1 and 2
# Hemophilia A
# Muscular dystrophy, Duchenne type
# Hypophosphatemic rickets
# Rett's syndrome

# LIST 4:
#
# Adrenoleukodystrophy (ALD)
# Agammaglobulinemia non-Bruton type	IGHM
# Alport syndrome	COL4A5
# Amyloid neuropathy – Andrade disease
# Angioneurotic oedema	C1NH
# Bartter syndrome type 4	BSND
# Blepharophimosis - ptosis - epicanthus inversus syndrome (BEPS)	FOXL2
# Brugada sindrome - Long QT syndrome-3	SCN5A
# Bruton agammaglobulinemia tyrosine kinase	BTK
# Ceroid lipofuscinosis neuronal type 2	CLN2
# Charcot Marie Tooth type 1A (CMT1A)	PMP22
# Charcot Marie Tooth type X (CMTX)	CMTX
# Chronic granulomatous disease (CGD)	CYBB
# Cystic Fibrosis  (CF)	CFTR
# Congenital adrenal hyperplasia (
#   Congenital disorder of glycosylation type Ia (CDG Ia)	PMM2
#   Congenital fibrosis of extraocular muscles 1 (CFEOM1)	KIF21A
#   Crigler-Najjar syndrome	UGT1A1
#   Deafness, autosomal recessive	CX26
#   Diamond-Blackfan anemia (DBA)	RPS19
#   Duchenne-Becker muscular dystrophy (DMD/DMB)	DMD
#   Duncan disease - X-linked lymphoproliferative syndrome (XLPD)	SH2D1A
#   Ectrodactyly ectodermal dysplasia and cleft lip/palate syndrome (EEC)
#   Epidermolysis bullosa dystrophica/pruriginosa	COL7A1
#   Exostoses multiple type I (EXT1)	EXT1
#   Exostoses multiple type II (EXT2)	EXT2
#   Facioscapulohumeral muscular dystrophy	FRG1
#   Factor VII
#   Familial Mediterranean Fever (FMF)	MEFV
#   Fanconi anemia A	FANCA
#   Fanconi anemia G	FANCG
#   Fragile-
#     Gangliosidosis (
#       Gaucher disease (GD)	GBA
#       Glanzmann thrombasthenia	ITGA2B
#       Glucose-6-phosphate dehydrogenase
#       Glutaric acidemia I	GCDH
#       Haemophilia
#       Haemophilia
#       Hand-foot-uterus syndrome	HOXD13
#       Hemophagocytic lymphohistiocytosis familial, type 2 (FHL2)	PRF1
#       Hypomagnesaemia primary	CLDN16
#       HYPOPHOSPHATASIA	ALPL
#       Holt-Oram Sindrome (HOS)	TBX5
#       Homocystinuria	MTHFR
#       Incontinentia pigmenti	NEMO
#       Lesch-Nyhan syndrome	HPRT
#       Limb-girdle muscular dystrophy type 2C (LGMD2C)	SGCG
#       Long QT syndrome-1	KCNQ1
#       Mannosidosis Alpha	MAN2B1
#       Marfan syndrome	FBN1
#       Methacrylic Aciduria, deficiency of beta-hydroxyisobutyryl-CoA deacylase	HIBCH
#       Mevalonic
#       Myotonic dystrophy (DM)
#       Myotonic dystrophy type 2 (DM2)	ZNF9
#       Mucopolysaccharidosis Type I - Hurler syndrome	IDUA
#       Mucopolysaccharidosis Type IIIA - Sanfilippo sindrome A (MPS3A)	SGSH
#       Mucopolysaccharidosis Type IIIB - Sanfilippo sindrome B (MPS3B)
#       Mucopolysaccharidosis Type VI (MPS VI) - Maroteaux-Lamy
#       Neuronal ceroid lipofuscinosis 1 - Batten's disease (CLN1)	PPT1
# Niemann-Pick disease
# Noonan sindrome	PTPN11
# Pancreatitis, hereditary (PCTT)	PRSS1
# Paramyotonia congenita (PMC)	SCN4A
# Phenylketonuria	PAH
# Polycystic kidney disease type 1 (
# Polycystic kidney disease type 2 (PKD2)	PKD2
# Polycystic kidney and hepatic disease-1 (ARPKD)	PKHD1
# Schwartz-Jampel/Stuve-Wiedemann syndrome	LIFR
# Sickle cell anemia
# Synpolydactyly (SPD1)	HOXA13
# Smith-Lemli-Opitz syndrome	DHCR7
# Spastic paraplegia type 3	SPG3A
# Spinal Muscular Atrophy (SMA)	SMN
# Spinocerebellar ataxia 3 (SCA3)	ATXN3
# Spinocerebellar ataxia 7 (SCA7)	ATXN7
# Stargardt disease	ABCA4
# Tay Sachs (TSD)	HEXA
# Thalassemia-α mental retardation syndrome	ATRX
# Thalassemia-
# Torsion dystonia, early onset (EOTD)	DYT1
# Tyrosinaemia type 1	FAH
# Tuberosclerosis 1	TSC1
# Tuberosclerosis 2	TSC2
# Wiskott-Aldrich Sindrome (WAS)	WAS

# PAPER SOURCE... A.NesterovaMonogenicRareDiseasesinBiomedicalDatabasesandTextMining.pdf
# This source is useful because its already reporting on the diseases which are most found in literature and thus are most likely to have representation on the GWAS catalog
#
# "Figure 2. Top ranked monogenic rare diseases by combination of point prevalence score and literature coverage."
#
# Cystic Fibrosis
# Hemophilia A
# Retinitis pigmentosa
# Hemophilia B
# Hereditary breast and ovarian cancer syndrome
# Huntington disease
# Neurofibromatosis type 1
# Duchenne muscular dystrophy
# Sickle cell anemia
# Tuberous sclerosis complex
# Charcot-Marie-Tooth disease / Hereditary motor...
# Marfan syndrome
# Fanconi anemia
# Titin-related limb-girdle muscular dystrophy R10
# ataxia-telangiectasia
# Severe combined immunodeficiency
# Lynch syndrome
# Congenital hypothyroidism
# Familial adenomatous polyposis
# Phenylketonuria
# Rett syndrome
# Microphthalmia-anophthalmia-coloboma
# Joubert syndrome
# Primary ciliary dyskinesia
# Autosomal dominant polycystic kidney disease
# Wilson disease
# Limb-girdlge muscular dystrophy
# Arryhthmogenic right entricular cardiomyopathy
# Chronic granulomatous disease
# Hereditary spastic paraplegia
# Familial thoracic aortic aneurysm and...
# Early infantile epileptic encephalopathy
# Alport Syndrome
#
#
#




# data inspection ---------------------------------------------------------
library(tidyverse)
numSnps <- length(unique(A$SNPS)) # 249,792

library(GWASpops.pheno2geno)
?createMT
ls(envir= as.environment("package:GWASpops.pheno2geno"))

library(utils)
package_native_routine_registration("GWASpops.pheno2geno")
ls(envir= as.environment("package:utils"))


curData = list.files('./GWAS.data')

mt1 <- createMT('./GWAS.data/air_pollution', varAnnotations = F, processData = F)

l = vector('list')
ll = vector('list')
for(f in curData){
  path <- paste0('./GWAS.data/', f)
  print(c(f,path))
  ll[f] = createMT(path, varAnnotations = F, processData = F) # worked the second time.. not sure what went wrong the first time. good to know
}

g = "dogman"
l2 = vector('list')
l2[g] = "hamshand"  # this works
l2

mt1 <- createMT('./GWAS.data/air_pollution', varAnnotations = F, processData = F)
mt2 <- createMT('./GWAS.data/alcohol_consumption', varAnnotations = F, processData = F)
mt3 <- createMT('./GWAS.data/breast_carcinoma', varAnnotations = F, processData = F)
mt4 <- createMT('./GWAS.data/colorectal_cancer', varAnnotations = F, processData = F)
mt5 <- createMT('./GWAS.data/inflammatory_bowel_disease', varAnnotations = F, processData = F)
mt6 <- createMT('./GWAS.data/Intelligence', varAnnotations = F, processData = F)
mt7 <- createMT('./GWAS.data/lung_cancer', varAnnotations = F, processData = F)
mt8 <- createMT('./GWAS.data/malabsorption_syndrome', varAnnotations = F, processData = F)
mt9 <- createMT('./GWAS.data/neuroticism_measurement', varAnnotations = F, processData = F)
mt10 <- createMT('./GWAS.data/prostate_cancer', varAnnotations = F, processData = F)
mt11 <- createMT('./GWAS.data/substance_abuse', varAnnotations = F, processData = F)

bindthem <- c(mt2,mt3,mt4,mt5,mt6,mt7,mt8,mt9,mt10,mt11)

for(t in bindthem){
  mt1 <- bind_rows(mt1, t)
}

for(i in 1:length(bindthem)){ # for some reason this produced 348k rows.. should be like 10k .. no idea why this might happen
  mt1 <- bind_rows(mt1, bindthem[i+1])
}

numSNPS2 <- length(unique(mt1$VariantID)) # 11013

ratioSNPS <- numSNPS2/numSnps

c(3,4,5) *(1/.044) # 68.18182  90.90909 113.63636 ~estimate of num hours to just call api for all data between 70-120 hrs = 3-5 days

?bind_rows



# NUM UNIQUE TRAITS

colnames(A)
# MAPPED_TRAIT_URI - URI of the EFO (experimental factor ontology) trait
# MAPPED_TRAIT - Mapped Experimental Factor Ontology trait for this study
# DISEASE/TRAIT - Disease or trait examined in the GWAS

MtURI <- A$MAPPED_TRAIT_URI
mapTrait <- A$MAPPED_TRAIT
diseaseTrait <- A$`DISEASE/TRAIT`

numMtURI <- length(unique(MtURI)) # 7479
numMT <- length(unique(mapTrait)) # 7479
numDT <- length(unique(diseaseTrait))# 21086

# why so many more unique disease/traits then mapped_traits? what is the difference?
# ANSWER: so the URI is some encoded number which references some term in the EFO framework. And thus the disease/trait is the more nuanced subdivision of those EFO terms. This explains why there are so many more disease traits then traits. Also why the num traitURI's matches the numTraits.
# I think it would make sense to search based on mapped traits then.. as that will return the sub-traits found within disease/trait col...

MtURI[1] #"http://www.ebi.ac.uk/efo/EFO_0000305"
mapTrait[1] # "breast carcinoma"

uniCF_MapTrait <- unique(CFrows$MAPPED_TRAIT)
uniCF_DisTrait <- unique(CFrows$`DISEASE/TRAIT`)

uniCF_MapTrait #"cystic fibrosis associated meconium ileus"                   "obsolete_cystic fibrosis, type 2 diabetes mellitus"
# [3] "lung disease severity measurement, obsolete_cystic fibrosis" "obsolete_cystic fibrosis"
# [5] "CFTR mutation carrier status"                                "cystic fibrosis-related diabetes"
#
uniCF_DisTrait #"Meconium ileus in cystic fibrosis"
# [2] "Cystic fibrosis-related diabetes"
# [3] "Lung disease severity in cystic fibrosis"
# [4] "Cystic fibrosis severity"
# [5] "CFTR mutation F508del heterozygosity in cystic fibrosis"
# [6] "CFTR mutation F508del heterozygosity in cystic fibrosis (PC-adjusted model)"
# [7] "Meconium Ileus in Cystic Fibrosis (adjusted for Cystic Fibrosis-Related Diabetes)"
# [8] "Cystic fibrosis-related diabetes in cystic fibrosis"
# [9] "Cystic Fibrosis-Related Diabetes in Cystic Fibrosis (adjusted for Meconium Ileus)"

matchesCF <- tolower(uniCF_DisTrait) %in% tolower(uniCF_MapTrait) # 1 match
matchesCF

# going to try and search the GWAS catalog to discover which search terms work well.. I want to know which col to search my monogenic traits by, and imagine whatever is mapped to the searching of GWASc itself would be better.. likely the less nuanced, as we can later cut out any traits which aren't interesting. ....
# One novel and major realization about the design of the capstone is that:
#   1. we need to choose a subset of traits/diseases to report upon to fulfill the second major presentation objective of the project, which would be comparing incidence of said trait in populations of interest to the SNP incidence. .. there is simply no way to fulfill this expansive degree of statistical literature search without narrowing down on a reasonably sized subset of traits we are interested in as subjects.
#   2. Some traits wouldn't make good candidates, such as intelligence.. as I just don't think that really suits the angle of the research questions as well as disease traits would.
#   3. Some traits within studies may be too nuanced to be compared to incidence statistics in the first place. Just looking at the reported disease/traits for CF, we can see that few can be said to report on incidence of CF itself, as much as they report on incidence of some associated trait that comes AFTER incidence is already confirmed. (like CF related diabetes for example.) Thus, I believe only by aiming at well suited disease/traits will we be able to make any sort of a logically interesting comparison between SNP-population incidence and population-to-population disease incidence. .... because of these difficulties, the process of identifying suitable targets for the project apears to be of paramount importance if we are to maintain the direction wrt the later presentation goals.
#   .... all of this needs to be discussed with Kulathinal, I could really use some guidance on how to approach these challenges.


# PLAN FOR SEEKING MONOGENIC D DATA:
#  search both MAPPED_TRAIT and DISEASE/TRAIT for our search terms... I think each search term may need its own regex written for optimal searching.
#  Write regex patterns for each. Load into array. hit associations for all matches for all regex patterns. split results based on level-1 = regex pattern / broad trait category, level-2 DISEASE/TRAIT col ..
#  final product should be a list with 2 levels thus which can be manually curated adhoc. Well call it the monogenicGWASlist
#
#
#  IMPORTANT: Func must be able to adjust caps sensitivity
#  function should be run once for caps sensitive terms and once for caps insensitive terms
#  FUNC SHOULD ALSO process sublists which may contain abbreviations which are always run w/ caps insensitivity

# No way to search background traits.
# cannot download such data.. so even when conditions such as CMT (charcot-Marie-Tooth disease) are the binding factor between a bunch of traits theres no way to actually search for it.




# OK we now have some data to report on for the monogenic disease investigation... We should breifly describe our search for this information and cite sources used... (that 1 paper I found) Describe the data processing and the results in an Rmd format. Also talk about the insights gained wrt other goals. Then shoot it over to Kulathinal in w/e format makes the most sense.. maybe a PDF


































































































































































