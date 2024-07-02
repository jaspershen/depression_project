###
no_source()

setwd(r4projects::get_project_wd())
setwd("data_analysis/metabolomics_data/")
library(tidyverse)
library(data.table)

data = readr::read_csv("processed_metabolites_data.csv")

colnames(data)

sample_info = data[,c(1:7)]
expression_data = data[-c(1:7)]

sample_info = 
  sample_info %>% 
  dplyr::select(-c("X1")) %>% 
  dplyr::rename(subject_id = id) %>% 
  dplyr::mutate(sample_id = paste(subject_id, Time, sep = "_")) %>% 
  dplyr::select(sample_id, subject_id, dplyr::everything())

variable_info = 
  data.frame(variable_id = colnames(expression_data))

annotation = readr::read_csv("210225 School_study_metabolomics_final_annotation_v2.csv")

raw_annotation = readr::read_csv("MSMS Annotation/annotation_result.csv") %>% 
  dplyr::select(name, mz, rt)


name = raw_annotation$name
name = 
  name %>% 
  stringr::str_replace("rplc_pos_|rplc_neg_|hilic_pos_|hilic_neg_", "")

idx = 
variable_info$variable_id %>% 
  stringr::str_replace("M", "") %>% 
  stringr::str_replace("m.z", "m/z") %>% 
  purrr::map(.f = function(x){
    # cat(x, " ")
    temp_idx = which(x == name)
    if(length(temp_idx) == 0){
      temp_idx = NA
    }
    
    if(length(temp_idx) > 1){
      temp_idx = temp_idx[1]
    }
    temp_idx
  }) %>% 
  unlist()

variable_info$mz = raw_annotation$mz[idx]
variable_info$rt = raw_annotation$rt[idx]

# annotation = 
# annotation %>% 
#   dplyr::mutate(new_id = 
#                   case_when(
#                     Mode == "pHILIC" ~ paste("hilic_pos", Compound, sep = "_"),
#                     Mode == "nHILIC" ~ paste("hilic_neg", Compound, sep = "_"),
#                     Mode == "pRPLC" ~ paste("rplc_pos", Compound, sep = "_"),
#                     Mode == "nRPLC" ~ paste("rplc_neg",Compound, sep = "_")
#                   ))
# 
# annotation = 
# annotation %>% 
#   dplyr::left_join(raw_annotation[,c("name", "mz", "rt")], by = c("new_id" = "name"))

annotation$Compound = 
  annotation$Compound %>% 
  paste("M", ., sep = "") %>% 
  stringr::str_replace("/", ".")

idx = which(is.na(match(annotation$Compound, variable_info$variable_id)))

variable_info = 
variable_info %>% 
  dplyr::left_join(annotation, by = c("variable_id" = "Compound"))

dir.create("data_preparation")
dir.create("data_preparation/peaks")
dir.create("data_preparation/metabolites")

expression_data = 
expression_data %>% 
  t() %>% 
  as.data.frame()

colnames(expression_data) = sample_info$sample_id
rownames(expression_data) == variable_info$variable_id

head(variable_info)

save(sample_info, file = "data_preparation/peaks/sample_info")
save(variable_info, file = "data_preparation/peaks/variable_info")
save(expression_data, file = "data_preparation/peaks/expression_data")

##only remain metabolites
variable_info = 
variable_info %>% 
  dplyr::filter(!is.na(Metabolite))

expression_data = 
  expression_data[variable_info$variable_id,]

save(variable_info, file = "data_preparation/metabolites/variable_info")
save(expression_data, file = "data_preparation/metabolites/expression_data")
save(sample_info, file = "data_preparation/metabolites/sample_info")
