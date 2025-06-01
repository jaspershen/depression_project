###
no_source()
library(tidyverse)
setwd(r4projects::get_project_wd())

rm(list = ls())

load("3_data_analysis/metabolomics_data/data_preparation/sample_info")

subject_info <-
  sample_info %>%
  dplyr::distinct(subject_id, .keep_all = TRUE) %>%
  dplyr::select(subject_id, age, sex, bdi_total, safe) %>%
  dplyr::mutate(subject_id = as.character(subject_id))

load("3_data_analysis/sample_info/sample_info")

library(tidyverse)
library(data.table)
library(tidymass)

load("3_data_analysis/transcriptomics/data_preparation/transcriptomics_data")

library(RNAAgeCalc)

dir.create("3_data_analysis/transcriptomics/age_prediction")
setwd("3_data_analysis/transcriptomics/age_prediction")

# #####T1 point
data <-
  transcriptomics_data %>%
  activate_mass_dataset(what = "sample_info") %>%
  dplyr::filter(stringr::str_detect(sample_id, "T1")) %>%
  extract_expression_data()

massdataset::export_mass_dataset(object = transcriptomics_data, file_type = "csv")

chronage <-
  data.frame(sampleid = colnames(data)) %>%
  dplyr::mutate(subject_id = stringr::str_extract(sampleid, "[0-9]{1,2}") %>% as.character()) %>%
  dplyr::left_join(subject_info, by = "subject_id") %>%
  dplyr::select(sampleid, age)
#
# t1_result <-
#   predict_age(
#     exprdata = data,
#     tissue = "blood",
#     exprtype = "FPKM",
#     chronage = chronage,
#     signature = "GenAge",
#     stype = "all"
#   )
# 
# t1_result <-
#   t1_result %>%
#   dplyr::mutate(sample_id = colnames(data))
# 
# write.csv(t1_result, file = "t1_result.csv", row.names = FALSE)
# save(t1_result, file = "t1_result")
# 
# #####T2 point
# data <-
#   transcriptomics_data %>%
#   activate_mass_dataset(what = "sample_info") %>%
#   dplyr::filter(stringr::str_detect(sample_id, "T2")) %>%
#   extract_expression_data()
# 
# chronage <-
#   data.frame(sampleid = colnames(data)) %>%
#   dplyr::mutate(subject_id = stringr::str_extract(sampleid, "[0-9]{1,2}") %>% as.character()) %>%
#   dplyr::left_join(subject_info, by = "subject_id") %>%
#   dplyr::select(sampleid, age)
# 
# t2_result <-
#   predict_age(
#     exprdata = data,
#     tissue = "blood",
#     exprtype = "FPKM",
#     chronage = chronage,
#     signature = "GenAge",
#     stype = "all"
#   )
# 
# t2_result <-
#   t2_result %>%
#   dplyr::mutate(sample_id = colnames(data))
# 
# write.csv(t2_result, file = "t2_result.csv", row.names = FALSE)
# save(t2_result, file = "t2_result")
# 
# #####T3 point
# data <-
#   transcriptomics_data %>%
#   activate_mass_dataset(what = "sample_info") %>%
#   dplyr::filter(stringr::str_detect(sample_id, "T3")) %>%
#   extract_expression_data()
# 
# chronage <-
#   data.frame(sampleid = colnames(data)) %>%
#   dplyr::mutate(subject_id = stringr::str_extract(sampleid, "[0-9]{1,2}") %>% as.character()) %>%
#   dplyr::left_join(subject_info, by = "subject_id") %>%
#   dplyr::select(sampleid, age)
# 
# t3_result <-
#   predict_age(
#     exprdata = data,
#     tissue = "blood",
#     exprtype = "FPKM",
#     chronage = chronage,
#     signature = "GenAge",
#     stype = "all"
#   )
# 
# t3_result <-
#   t3_result %>%
#   dplyr::mutate(sample_id = colnames(data))
# 
# write.csv(t3_result, file = "t3_result.csv", row.names = FALSE)
# save(t3_result, file = "t3_result")
# 
# #####T4 point
# data <-
#   transcriptomics_data %>%
#   activate_mass_dataset(what = "sample_info") %>%
#   dplyr::filter(stringr::str_detect(sample_id, "T4")) %>%
#   extract_expression_data()
# 
# chronage <-
#   data.frame(sampleid = colnames(data)) %>%
#   dplyr::mutate(subject_id = stringr::str_extract(sampleid, "[0-9]{1,2}") %>% as.character()) %>%
#   dplyr::left_join(subject_info, by = "subject_id") %>%
#   dplyr::select(sampleid, age)
# 
# t4_result <-
#   predict_age(
#     exprdata = data,
#     tissue = "blood",
#     exprtype = "FPKM",
#     chronage = chronage,
#     signature = "GenAge",
#     stype = "all"
#   )
# 
# t4_result <-
#   t4_result %>%
#   dplyr::mutate(sample_id = colnames(data))
# 
# write.csv(t4_result, file = "t4_result.csv", row.names = FALSE)
# save(t4_result, file = "t4_result")
# 
# 
# #####T5 point
# data <-
#   transcriptomics_data %>%
#   activate_mass_dataset(what = "sample_info") %>%
#   dplyr::filter(stringr::str_detect(sample_id, "T5")) %>%
#   extract_expression_data()
# 
# chronage <-
#   data.frame(sampleid = colnames(data)) %>%
#   dplyr::mutate(subject_id = stringr::str_extract(sampleid, "[0-9]{1,2}") %>% as.character()) %>%
#   dplyr::left_join(subject_info, by = "subject_id") %>%
#   dplyr::select(sampleid, age)
# 
# t5_result <-
#   predict_age(
#     exprdata = data,
#     tissue = "blood",
#     exprtype = "FPKM",
#     chronage = chronage,
#     signature = "GenAge",
#     stype = "all"
#   )
# 
# t5_result <-
#   t5_result %>%
#   dplyr::mutate(sample_id = colnames(data))
# 
# write.csv(t5_result, file = "t5_result.csv", row.names = FALSE)
# save(t5_result, file = "t5_result")
# 
# #####T6 point
# data <-
#   transcriptomics_data %>%
#   activate_mass_dataset(what = "sample_info") %>%
#   dplyr::filter(stringr::str_detect(sample_id, "T6")) %>%
#   extract_expression_data()
# 
# chronage <-
#   data.frame(sampleid = colnames(data)) %>%
#   dplyr::mutate(subject_id = stringr::str_extract(sampleid, "[0-9]{1,2}") %>% as.character()) %>%
#   dplyr::left_join(subject_info, by = "subject_id") %>%
#   dplyr::select(sampleid, age)
# 
# t6_result <-
#   predict_age(
#     exprdata = data,
#     tissue = "blood",
#     exprtype = "FPKM",
#     chronage = chronage,
#     signature = "GenAge",
#     stype = "all"
#   )
# 
# t6_result <-
#   t6_result %>%
#   dplyr::mutate(sample_id = colnames(data))
# 
# write.csv(t6_result, file = "t6_result.csv", row.names = FALSE)
# save(t6_result, file = "t6_result")


###all samples
data <-
  transcriptomics_data %>%
  activate_mass_dataset(what = "sample_info") %>%
  extract_expression_data()

chronage <-
  data.frame(sampleid = colnames(data)) %>%
  dplyr::mutate(subject_id = stringr::str_extract(sampleid, "[0-9]{1,2}") %>% as.character()) %>%
  dplyr::left_join(subject_info, by = "subject_id") %>%
  dplyr::select(sampleid, age)

result <-
  predict_age(
    exprdata = data,
    tissue = "blood",
    exprtype = "FPKM",
    chronage = chronage,
    signature = "GenAge",
    stype = "all"
  )

write.csv(result, file = "result.csv", row.names = FALSE)
save(result, file = "result")
