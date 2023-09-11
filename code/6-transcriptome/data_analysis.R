###
no_source()

setwd(masstools::get_project_wd())

rm(list = ls())

load("data_analysis/metabolomics_data/data_preparation/sample_info")

subject_info <-
  sample_info %>%
  dplyr::distinct(subject_id, .keep_all = TRUE) %>%
  dplyr::select(subject_id, age, sex, bdi_total, safe) %>%
  dplyr::mutate(subject_id = as.character(subject_id))

library(tidyverse)
library(data.table)
library(tidymass)

load("data_analysis/transcriptomics/data_preparation/transcriptomics_data")

library(RNAAgeCalc)

dir.create("data_analysis/transcriptomics/age_prediction")
setwd("data_analysis/transcriptomics/age_prediction")

###only remain the protein-code gene
transcriptomics_data <-
  transcriptomics_data %>%
  activate_mass_dataset(what = "variable_info") %>%
  dplyr::filter(!is.na(GENETYPE)) %>%
  dplyr::filter(GENETYPE == "protein-coding")

# #####T1 point
# data <-
#   transcriptomics_data %>%
#   activate_mass_dataset(what = "sample_info") %>%
#   dplyr::filter(stringr::str_detect(sample_id, "T1")) %>%
#   extract_expression_data()
#
# chronage <-
#   data.frame(sampleid = colnames(data)) %>%
#   dplyr::mutate(subject_id = stringr::str_extract(sampleid, "[0-9]{1,2}") %>% as.character()) %>%
#   dplyr::left_join(subject_info, by = "subject_id") %>%
#   dplyr::select(sampleid, age)
#
# t1_result <-
#   predict_age(
#     exprdata = data,
#     tissue = "blood",
#     exprtype = "FPKM",
#     chronage = chronage
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
#     chronage = chronage
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
#     chronage = chronage
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
#     chronage = chronage
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
#     chronage = chronage
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
#     chronage = chronage
#   )
#
# t6_result <-
#   t6_result %>%
#   dplyr::mutate(sample_id = colnames(data))
#
# write.csv(t6_result, file = "t6_result.csv", row.names = FALSE)
# save(t6_result, file = "t6_result")

load("t1_result")
load("t2_result")
load("t3_result")
load("t4_result")
load("t5_result")
load("t6_result")

###summary
t1_result <-
  t1_result %>%
  dplyr::mutate(subject_id = stringr::str_extract(sample_id, "[0-9]{1,2}")) %>%
  dplyr::arrange(subject_id)

t2_result <-
  t2_result %>%
  dplyr::mutate(subject_id = stringr::str_extract(sample_id, "[0-9]{1,2}")) %>%
  dplyr::arrange(subject_id)

t3_result <-
  t3_result %>%
  dplyr::mutate(subject_id = stringr::str_extract(sample_id, "[0-9]{1,2}")) %>%
  dplyr::arrange(subject_id)

t4_result <-
  t4_result %>%
  dplyr::mutate(subject_id = stringr::str_extract(sample_id, "[0-9]{1,2}")) %>%
  dplyr::arrange(subject_id)

t5_result <-
  t5_result %>%
  dplyr::mutate(subject_id = stringr::str_extract(sample_id, "[0-9]{1,2}")) %>%
  dplyr::arrange(subject_id)

t6_result <-
  t6_result %>%
  dplyr::mutate(subject_id = stringr::str_extract(sample_id, "[0-9]{1,2}")) %>%
  dplyr::arrange(subject_id)

t1_result <-
  t1_result %>%
  dplyr::mutate(age_diff = RNAAge - ChronAge,
                time_point = "T1") %>%
  dplyr::select(subject_id, age_diff, time_point)

t2_result <-
  t2_result %>%
  dplyr::mutate(age_diff = RNAAge - ChronAge,
                time_point = "T2") %>%
  dplyr::select(subject_id, age_diff, time_point)


t3_result <-
  t3_result %>%
  dplyr::mutate(age_diff = RNAAge - ChronAge,
                time_point = "T3") %>%
  dplyr::select(subject_id, age_diff, time_point)

t4_result <-
  t4_result %>%
  dplyr::mutate(age_diff = RNAAge - ChronAge,
                time_point = "T4") %>%
  dplyr::select(subject_id, age_diff, time_point)

t5_result <-
  t5_result %>%
  dplyr::mutate(age_diff = RNAAge - ChronAge,
                time_point = "T5") %>%
  dplyr::select(subject_id, age_diff, time_point)

t6_result <-
  t6_result %>%
  dplyr::mutate(age_diff = RNAAge - ChronAge,
                time_point = "T6") %>%
  dplyr::select(subject_id, age_diff, time_point)

temp_data <-
  rbind(t1_result,
        t2_result,
        t3_result,
        t4_result,
        t5_result,
        t6_result)

dim(t1_result)
dim(t2_result)
dim(t3_result)
dim(t4_result)
dim(t5_result)
dim(t6_result)

library(plyr)

###remove the subjects without T1
temp_data <-
  temp_data %>%
  dplyr::filter(subject_id %in% t1_result$subject_id)

temp_data <-
  temp_data %>%
  plyr::dlply(.variables = .(subject_id)) %>%
  purrr::map(function(x) {
    x <-
      x %>%
      dplyr::arrange(time_point)
    x$age_diff <-
      x$age_diff - x$age_diff[x$time_point == "T1"]
    x
  }) %>%
  do.call(rbind, .) %>%
  as.data.frame()

temp_data %>%
  ggplot(aes(time_point, age_diff)) +
  geom_boxplot() +
  geom_jitter()

temp_data %>%
  ggplot(aes(time_point, age_diff)) +
  geom_line(aes(group = subject_id,
                color = subject_id)) +
  theme_bw()
