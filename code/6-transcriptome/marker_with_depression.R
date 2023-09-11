###
no_source()

setwd(masstools::get_project_wd())
rm(list = ls())

load("data_analysis/transcriptomics/data_preparation/transcriptomics_data")

dir.create("data_analysis/transcriptomics/marker_depression")
setwd("data_analysis/transcriptomics/marker_depression")

library(tidyverse)
library(data.table)

dim(transcriptomics_data)

###only remain protein-coding
library(tidymass)
transcriptomics_data <-
  transcriptomics_data %>%
  activate_mass_dataset(what = "variable_info") %>%
  dplyr::filter(!is.na(GENETYPE)) %>%
  dplyr::filter(GENETYPE == "protein-coding")

colnames(transcriptomics_data)
