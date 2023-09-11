no_function()
library(plyr)
library(tidyverse)

setwd(masstools::get_project_wd())
rm(list = ls())

setwd("data_analysis/Cytokines/")

data = readr::read_csv("corrected_cytokines_data.csv")

colnames(data)

sample_info = data[,c(1:8)]
expression_data = data[-c(1:8)]

sample_info = 
  sample_info %>% 
  dplyr::select(-c("X1")) %>% 
  dplyr::rename(subject_id = id, sample_id = Sample) %>% 
  dplyr::select(sample_id, subject_id, dplyr::everything())

sample_info$sample_id_old = sample_info$sample_id
sample_info$sample_id = paste(sample_info$subject_id, 
                              sample_info$Time, 
                              sep = "_")

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
