###
no_source()
setwd(r4projects::get_project_wd())

rm(list = ls())
library(tidyverse)

load("3-data_analysis/metabolomics_data/data_preparation/sample_info")

subject_info <-
  sample_info %>%
  dplyr::distinct(subject_id, .keep_all = TRUE) %>%
  dplyr::select(subject_id, age, sex, bdi_total, safe) %>%
  dplyr::mutate(subject_id = as.character(subject_id))

library(tidyverse)
library(data.table)
library(tidymass)

load("3-data_analysis/transcriptomics/data_preparation/transcriptomics_data")

library(RNAAgeCalc)

dir.create("3-data_analysis/transcriptomics/age_prediction")
setwd("3-data_analysis/transcriptomics/age_prediction")

#####signature parameter

#####T1 point
data <-
  transcriptomics_data %>%
  activate_mass_dataset(what = "sample_info") %>%
  dplyr::filter(stringr::str_detect(sample_id, "T1")) %>%
  extract_expression_data()

chronage <-
  data.frame(sampleid = colnames(data)) %>%
  dplyr::mutate(subject_id = stringr::str_extract(sampleid, "[0-9]{1,2}") %>% as.character()) %>%
  dplyr::left_join(subject_info, by = "subject_id") %>%
  dplyr::select(sampleid, age)

t1_optimization_result <-
  purrr::map(c(
    "DESeq2",
    "Pearson",
    "Dev",
    "deMagalhaes",
    "GenAge",
    "GTExAge",
    "Peters",
    "all"
  ), function(x) {
    temp <-
      predict_age(
        exprdata = data,
        tissue = "blood",
        exprtype = "FPKM",
        chronage = chronage,
        signature = x
      )
    data.frame(signature = x,
               mean = mean(abs(temp$RNAAge - temp$ChronAge)),
               sd = sd(abs(temp$RNAAge - temp$ChronAge)))
  }) %>%
  do.call(rbind, .) %>%
  as.data.frame()

t1_optimization_result


#####T2 point
data <-
  transcriptomics_data %>%
  activate_mass_dataset(what = "sample_info") %>%
  dplyr::filter(stringr::str_detect(sample_id, "T2")) %>%
  extract_expression_data()

chronage <-
  data.frame(sampleid = colnames(data)) %>%
  dplyr::mutate(subject_id = stringr::str_extract(sampleid, "[0-9]{1,2}") %>% as.character()) %>%
  dplyr::left_join(subject_info, by = "subject_id") %>%
  dplyr::select(sampleid, age)

t2_optimization_result <-
  purrr::map(c(
    "DESeq2",
    "Pearson",
    "Dev",
    "deMagalhaes",
    "GenAge",
    "GTExAge",
    "Peters",
    "all"
  ), function(x) {
    temp <-
      predict_age(
        exprdata = data,
        tissue = "blood",
        exprtype = "FPKM",
        chronage = chronage,
        signature = x
      )
    data.frame(signature = x,
               mean = mean(abs(temp$RNAAge - temp$ChronAge)),
               sd = sd(abs(temp$RNAAge - temp$ChronAge)))
  }) %>%
  do.call(rbind, .) %>%
  as.data.frame()

#####T3 point
data <-
  transcriptomics_data %>%
  activate_mass_dataset(what = "sample_info") %>%
  dplyr::filter(stringr::str_detect(sample_id, "T3")) %>%
  extract_expression_data()

chronage <-
  data.frame(sampleid = colnames(data)) %>%
  dplyr::mutate(subject_id = stringr::str_extract(sampleid, "[0-9]{1,2}") %>% as.character()) %>%
  dplyr::left_join(subject_info, by = "subject_id") %>%
  dplyr::select(sampleid, age)

t3_optimization_result <-
  purrr::map(c(
    "DESeq2",
    "Pearson",
    "Dev",
    "deMagalhaes",
    "GenAge",
    "GTExAge",
    "Peters",
    "all"
  ), function(x) {
    temp <-
      predict_age(
        exprdata = data,
        tissue = "blood",
        exprtype = "FPKM",
        chronage = chronage,
        signature = x
      )
    data.frame(signature = x,
               mean = mean(abs(temp$RNAAge - temp$ChronAge)),
               sd = sd(abs(temp$RNAAge - temp$ChronAge)))
  }) %>%
  do.call(rbind, .) %>%
  as.data.frame()

t1_optimization_result
t2_optimization_result
t3_optimization_result


####signature GenAge is best


#####stype parameter

#####T1 point
data <-
  transcriptomics_data %>%
  activate_mass_dataset(what = "sample_info") %>%
  dplyr::filter(stringr::str_detect(sample_id, "T1")) %>%
  extract_expression_data()

chronage <-
  data.frame(sampleid = colnames(data)) %>%
  dplyr::mutate(subject_id = stringr::str_extract(sampleid, "[0-9]{1,2}") %>% as.character()) %>%
  dplyr::left_join(subject_info, by = "subject_id") %>%
  dplyr::select(sampleid, age)

t1_optimization_result <-
  purrr::map(c("all", "caucasian"), function(x) {
    temp <-
      predict_age(
        exprdata = data,
        tissue = "blood",
        exprtype = "FPKM",
        chronage = chronage,
        stype = x
      )
    data.frame(stype = x,
               mean = mean(abs(temp$RNAAge - temp$ChronAge)),
               sd = sd(abs(temp$RNAAge - temp$ChronAge)))
  }) %>%
  do.call(rbind, .) %>%
  as.data.frame()

#####T2 point
data <-
  transcriptomics_data %>%
  activate_mass_dataset(what = "sample_info") %>%
  dplyr::filter(stringr::str_detect(sample_id, "T2")) %>%
  extract_expression_data()

chronage <-
  data.frame(sampleid = colnames(data)) %>%
  dplyr::mutate(subject_id = stringr::str_extract(sampleid, "[0-9]{1,2}") %>% as.character()) %>%
  dplyr::left_join(subject_info, by = "subject_id") %>%
  dplyr::select(sampleid, age)

t2_optimization_result <-
  purrr::map(c("all", "caucasian"), function(x) {
    temp <-
      predict_age(
        exprdata = data,
        tissue = "blood",
        exprtype = "FPKM",
        chronage = chronage,
        stype = x
      )
    data.frame(stype = x,
               mean = mean(abs(temp$RNAAge - temp$ChronAge)),
               sd = sd(abs(temp$RNAAge - temp$ChronAge)))
  }) %>%
  do.call(rbind, .) %>%
  as.data.frame()

#####T3 point
data <-
  transcriptomics_data %>%
  activate_mass_dataset(what = "sample_info") %>%
  dplyr::filter(stringr::str_detect(sample_id, "T3")) %>%
  extract_expression_data()

chronage <-
  data.frame(sampleid = colnames(data)) %>%
  dplyr::mutate(subject_id = stringr::str_extract(sampleid, "[0-9]{1,2}") %>% as.character()) %>%
  dplyr::left_join(subject_info, by = "subject_id") %>%
  dplyr::select(sampleid, age)

t3_optimization_result <-
  purrr::map(c("all", "caucasian"), function(x) {
    temp <-
      predict_age(
        exprdata = data,
        tissue = "blood",
        exprtype = "FPKM",
        chronage = chronage,
        stype = x
      )
    data.frame(stype = x,
               mean = mean(abs(temp$RNAAge - temp$ChronAge)),
               sd = sd(abs(temp$RNAAge - temp$ChronAge)))
  }) %>%
  do.call(rbind, .) %>%
  as.data.frame()


t1_optimization_result
t2_optimization_result
t3_optimization_result

####stype caucasian is best
