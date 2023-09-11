no_function()
library(plyr)
library(tidyverse)

setwd(masstools::get_project_wd())
rm(list = ls())

setwd("data_analysis/Lipids/")

data = readr::read_csv("processed_lipid_data.csv")

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

dir.create("data_preparation")

expression_data = 
  t(expression_data) %>% 
  as.data.frame()

colnames(expression_data) = sample_info$sample_id
rownames(expression_data) == variable_info$variable_id


save(sample_info, file = "data_preparation/sample_info")
save(variable_info, file = "data_preparation/variable_info")
save(expression_data, file = "data_preparation/expression_data")
