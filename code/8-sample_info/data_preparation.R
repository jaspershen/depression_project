no_function()
library(plyr)
library(tidyverse)

setwd(masstools::get_project_wd())
rm(list = ls())

data <- readr::read_csv("data_analysis/Psychometrics/Psychometric Data.csv")

dir.create("data_analysis/sample_info/", recursive = TRUE)
setwd("data_analysis/sample_info/")


sample_info = data

sample_info$sample_id <- 
  paste(sample_info$id, sample_info$Time, sep = "_")

sample_info <- 
sample_info %>% 
  dplyr::rename(subject_id = id) %>% 
  dplyr::select(sample_id, subject_id, record_id, dplyr::everything())


save(sample_info, file = "sample_info")
