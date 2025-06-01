setwd(r4projects::get_project_wd())

setwd("3_data_analysis/for_ariel_20220714/")

metabolite_data <- 
  readr::read_csv("processed_metabolites_data.csv")

metabolite_annotation1 <- 
  readxl::read_xlsx("Metabolomics_Annotated_SignificantHits.xlsx", sheet = 1)

metabolite_annotation2 <- 
  readxl::read_xlsx("Metabolomics_Annotated_SignificantHits.xlsx", sheet = 2)

metabolite_annotation3 <- 
  readxl::read_xlsx("Metabolomics_Annotated_SignificantHits.xlsx", sheet = 3)

metabolite_annotation4 <- 
  readxl::read_xlsx("Metabolomics_Annotated_SignificantHits.xlsx", sheet = 4)

metabolite_annotation5 <- 
  readxl::read_xlsx("Metabolomics_Annotated_SignificantHits.xlsx", sheet = 5)

###Eicosapentaenoic acid (EPA)
epa_id <- 
  metabolite_annotation3$Variables[which(metabolite_annotation3$Compounds_ID == "301.217_11.4")]

##Docosahexaenoic acid (DHA)
dha_id <-
metabolite_annotation1$Variables[which(metabolite_annotation1$Compounds_ID...7 == "327.2327_11.9")]

###Arachidonic acid
aa_id <-
metabolite_annotation1$Variables[which(metabolite_annotation1$Compounds_ID...7 == "303.2327_11.9")]

library(tidyverse)

colnames(metabolite_data) <-
  colnames(metabolite_data) %>% 
  stringr::str_replace("\\/", "\\.")

idx1 <- 
match(epa_id, colnames(metabolite_data))
idx2 <-
match(dha_id, colnames(metabolite_data))
idx3 <- 
match(aa_id, colnames(metabolite_data))

temp_data <-
  metabolite_data[,c(1:6, c(idx1, idx2, idx3))]
  


