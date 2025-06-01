no_function()
library(plyr)
library(tidyverse)

setwd(r4projects::get_project_wd())
rm(list = ls())

data <- readr::read_csv("3_data_analysis/Psychometrics/Psychometric Data.csv")

dir.create("3_data_analysis/sample_info/", recursive = TRUE)
setwd("3_data_analysis/sample_info/")


sample_info = data

sample_info$sample_id <- 
  paste(sample_info$id, sample_info$Time, sep = "_")

sample_info <- 
sample_info %>% 
  dplyr::rename(subject_id = id) %>% 
  dplyr::select(sample_id, subject_id, record_id, dplyr::everything())


save(sample_info, file = "sample_info")
