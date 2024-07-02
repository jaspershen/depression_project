no_function()
library(plyr)
library(tidyverse)

setwd(r4projects::get_project_wd())
rm(list = ls())

setwd("3-data_analysis/Psychometrics/")

data = readr::read_csv("Psychometric Data.csv")

sample_info = data[,c(1:8)]
expression_data = data[-c(1:8)]

sample_info = 
  sample_info %>% 
  dplyr::select(-c("X1")) %>% 
  dplyr::rename(subject_id = id, sample_id = Sample) %>% 
  dplyr::select(sample_id, subject_id, dplyr::everything())

variable_info = 
  data.frame(variable_id = colnames(expression_data))

dir.create("data_preparation")

save(sample_info, file = "data_preparation/sample_info")
save(variable_info, file = "data_preparation/variable_info")
save(expression_data, file = "data_preparation/expression_data")
