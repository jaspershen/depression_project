###
no_source()
library(tidyverse)
setwd(r4projects::get_project_wd())

rm(list = ls())

load("data_analysis/metabolomics_data/data_preparation/sample_info")

subject_info <-
  sample_info %>%
  dplyr::distinct(subject_id, .keep_all = TRUE) %>%
  dplyr::select(subject_id, age, sex, bdi_total, safe) %>%
  dplyr::mutate(subject_id = as.character(subject_id))

load("data_analysis/sample_info/sample_info")

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
  activate_mass_dataset(what = "variable_info")
# dplyr::filter(!is.na(GENETYPE)) %>%
# dplyr::filter(GENETYPE == "protein-coding")

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
#     stype = "caucasian"
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
#     stype = "caucasian"
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
#     stype = "caucasian"
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
#     stype = "caucasian"
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
#     stype = "caucasian"
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
#     stype = "caucasian"
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
  dplyr::select(subject_id, age_diff, time_point, AgeAccelResid)

t2_result <-
  t2_result %>%
  dplyr::mutate(age_diff = RNAAge - ChronAge,
                time_point = "T2") %>%
  dplyr::select(subject_id, age_diff, time_point, AgeAccelResid)

t3_result <-
  t3_result %>%
  dplyr::mutate(age_diff = RNAAge - ChronAge,
                time_point = "T3") %>%
  dplyr::select(subject_id, age_diff, time_point, AgeAccelResid)

t4_result <-
  t4_result %>%
  dplyr::mutate(age_diff = RNAAge - ChronAge,
                time_point = "T4") %>%
  dplyr::select(subject_id, age_diff, time_point) %>%
  dplyr::mutate(AgeAccelResid = NA)

t5_result <-
  t5_result %>%
  dplyr::mutate(age_diff = RNAAge - ChronAge,
                time_point = "T5") %>%
  dplyr::select(subject_id, age_diff, time_point) %>%
  dplyr::mutate(AgeAccelResid = NA)

t6_result <-
  t6_result %>%
  dplyr::mutate(age_diff = RNAAge - ChronAge,
                time_point = "T6") %>%
  dplyr::select(subject_id, age_diff, time_point) %>%
  dplyr::mutate(AgeAccelResid = NA)

temp_data <-
  rbind(t1_result,
        t2_result,
        t3_result,
        t4_result,
        t5_result,
        t6_result)

library(plyr)

###remove the subjects without T1
temp_data <-
  temp_data %>%
  dplyr::filter(subject_id %in% t1_result$subject_id)

# temp_data <-
#   temp_data %>%
#   plyr::dlply(.variables = .(subject_id)) %>%
#   purrr::map(function(x) {
#     x <-
#       x %>%
#       dplyr::arrange(time_point)
#     x$age_diff <-
#       x$age_diff - x$age_diff[x$time_point == "T1"]
#     x
#   }) %>%
#   do.call(rbind, .) %>%
#   as.data.frame()

temp_data <-
  temp_data %>%
  dplyr::filter(subject_id %in% unique(temp_data$subject_id[duplicated(temp_data$subject_id)]))

depressed_subject_id <-
  sample_info %>%
  dplyr::filter(Time == "T1") %>%
  dplyr::filter(depressed == "Depressed") %>%
  pull(subject_id)

# temp_data %>%
#   ggplot(aes(time_point, age_diff)) +
#   geom_boxplot() +
#   geom_jitter()
#
# temp_data %>%
#   dplyr::filter(subject_id %in% depressed_subject_id) %>%
#   ggplot(aes(time_point, age_diff)) +
#   geom_boxplot() +
#   geom_jitter()
#
# temp_data %>%
#   dplyr::filter(subject_id %in% depressed_subject_id) %>%
#   ggplot(aes(time_point, age_diff)) +
#   geom_line(aes(group = subject_id,
#                 color = subject_id),
#             show.legend = FALSE) +
#   theme_bw()

# plot <-
#   temp_data %>%
#   dplyr::filter(subject_id %in% depressed_subject_id) %>%
#   dplyr::mutate(subject_id = factor(subject_id, levels = stringr::str_sort(unique(
#     temp_data$subject_id
#   ), numeric = TRUE))) %>%
#   ggplot(aes(time_point, age_diff)) +
#   geom_line(aes(group = subject_id,
#                 color = subject_id),
#             show.legend = FALSE) +
#   geom_point(aes(group = subject_id,
#                  color = subject_id),
#              show.legend = FALSE) +
#   theme_bw() +
#   labs(x = "",
#        y = "Age difference (Biological age - real age, years)") +
#   facet_wrap(facets = vars(subject_id))
#
# plot
# ggsave(plot,
#        filename = "age_difference_each_person.pdf",
#        width = 7,
#        height = 6)
# ggsave(plot,
#        filename = "age_difference_each_person.png",
#        width = 7,
#        height = 6)
###only remain the depressed subjects

# ###T1 vs T2
# temp1 <-
#   t1_result %>%
#   dplyr::filter(subject_id %in% depressed_subject_id) %>%
#   dplyr::filter(subject_id %in% intersect(t1_result$subject_id,
#                                           t2_result$subject_id)) %>%
#   dplyr::arrange(subject_id)
#
# temp2 <-
#   t2_result %>%
#   dplyr::filter(subject_id %in% depressed_subject_id) %>%
#   dplyr::filter(subject_id %in% intersect(t1_result$subject_id,
#                                           t2_result$subject_id)) %>%
#   dplyr::arrange(subject_id)
#
# t.test(temp1$age_diff,
#        temp2$age_diff,
#        paired = TRUE)
#
# ###T1 vs T3
# temp1 <-
#   t1_result %>%
#   dplyr::filter(subject_id %in% depressed_subject_id) %>%
#   dplyr::filter(subject_id %in% intersect(t1_result$subject_id,
#                                           t3_result$subject_id)) %>%
#   dplyr::arrange(subject_id)
#
# temp3 <-
#   t3_result %>%
#   dplyr::filter(subject_id %in% depressed_subject_id) %>%
#   dplyr::filter(subject_id %in% intersect(t1_result$subject_id,
#                                           t3_result$subject_id)) %>%
#   dplyr::arrange(subject_id)
#
# t.test(temp1$age_diff,
#        temp3$age_diff,
#        paired = TRUE)
#
# ###T1 vs T4
# temp1 <-
#   t1_result %>%
#   dplyr::filter(subject_id %in% depressed_subject_id) %>%
#   dplyr::filter(subject_id %in% intersect(t1_result$subject_id,
#                                           t4_result$subject_id)) %>%
#   dplyr::arrange(subject_id)
#
# temp4 <-
#   t4_result %>%
#   dplyr::filter(subject_id %in% depressed_subject_id) %>%
#   dplyr::filter(subject_id %in% intersect(t1_result$subject_id,
#                                           t4_result$subject_id)) %>%
#   dplyr::arrange(subject_id)
#
# t.test(temp1$age_diff,
#        temp4$age_diff,
#        paired = TRUE)
#
#
# ###T1 vs T5
# temp1 <-
#   t1_result %>%
#   dplyr::filter(subject_id %in% depressed_subject_id) %>%
#   dplyr::filter(subject_id %in% intersect(t1_result$subject_id,
#                                           t5_result$subject_id)) %>%
#   dplyr::arrange(subject_id)
#
# temp5 <-
#   t5_result %>%
#   dplyr::filter(subject_id %in% depressed_subject_id) %>%
#   dplyr::filter(subject_id %in% intersect(t1_result$subject_id,
#                                           t5_result$subject_id)) %>%
#   dplyr::arrange(subject_id)
#
# t.test(temp1$age_diff,
#        temp5$age_diff,
#        paired = TRUE)
#
# ###T1 vs T6
# temp1 <-
#   t1_result %>%
#   dplyr::filter(subject_id %in% depressed_subject_id) %>%
#   dplyr::filter(subject_id %in% intersect(t1_result$subject_id,
#                                           t6_result$subject_id)) %>%
#   dplyr::arrange(subject_id)
#
# temp6 <-
#   t6_result %>%
#   dplyr::filter(subject_id %in% depressed_subject_id) %>%
#   dplyr::filter(subject_id %in% intersect(t1_result$subject_id,
#                                           t6_result$subject_id)) %>%
#   dplyr::arrange(subject_id)
#
# t.test(temp1$age_diff,
#        temp6$age_diff,
#        paired = TRUE)
#
# ###T2 vs T3
# temp2 <-
#   t2_result %>%
#   dplyr::filter(subject_id %in% depressed_subject_id) %>%
#   dplyr::filter(subject_id %in% intersect(t2_result$subject_id,
#                                           t3_result$subject_id)) %>%
#   dplyr::arrange(subject_id)
#
# temp3 <-
#   t3_result %>%
#   dplyr::filter(subject_id %in% depressed_subject_id) %>%
#   dplyr::filter(subject_id %in% intersect(t2_result$subject_id,
#                                           t3_result$subject_id)) %>%
#   dplyr::arrange(subject_id)
#
# t.test(temp2$age_diff,
#        temp3$age_diff,
#        paired = TRUE)
#
# ###T2 vs T4
# temp2 <-
#   t2_result %>%
#   dplyr::filter(subject_id %in% depressed_subject_id) %>%
#   dplyr::filter(subject_id %in% intersect(t2_result$subject_id,
#                                           t4_result$subject_id)) %>%
#   dplyr::arrange(subject_id)
#
# temp4 <-
#   t4_result %>%
#   dplyr::filter(subject_id %in% depressed_subject_id) %>%
#   dplyr::filter(subject_id %in% intersect(t2_result$subject_id,
#                                           t4_result$subject_id)) %>%
#   dplyr::arrange(subject_id)
#
# t.test(temp2$age_diff,
#        temp4$age_diff,
#        paired = TRUE)


####AgeAccelResid
temp_data2 <-
  rbind(t1_result,
        t2_result,
        t3_result,
        t4_result,
        t5_result,
        t6_result) %>%
  dplyr::filter(!is.na(AgeAccelResid))

###only remain the subject2 with at least two samples
temp_data2 <-
  temp_data2 %>%
  dplyr::filter(subject_id %in% unique(temp_data2$subject_id[duplicated(temp_data2$subject_id)]))

##only remain the subjects with acceleration > 0 at time point 1
remain_subject_id <-
  temp_data2 %>%
  dplyr::filter(time_point == "T1") %>%
  dplyr::filter(AgeAccelResid > 0) %>%
  pull(subject_id)

temp_data2 <-
  temp_data2 %>%
  dplyr::filter(subject_id %in% remain_subject_id)

# temp_data2 %>%
#   dplyr::filter(subject_id %in% depressed_subject_id) %>%
#   ggplot(aes(time_point, AgeAccelResid)) +
#   geom_boxplot() +
#   geom_jitter()
#
# temp_data2 %>%
#   dplyr::filter(subject_id %in% depressed_subject_id) %>%
#   ggplot(aes(time_point, AgeAccelResid)) +
#   geom_line(aes(group = subject_id,
#                 color = subject_id),
#             show.legend = FALSE) +
#   theme_bw()

plot <-
  temp_data2 %>%
  # dplyr::filter(subject_id %in% depressed_subject_id) %>%
  dplyr::mutate(subject_id = factor(subject_id, levels = stringr::str_sort(unique(
    temp_data2$subject_id
  ), numeric = TRUE))) %>%
  dplyr::filter(!is.na(AgeAccelResid)) %>%
  ggplot(aes(time_point, AgeAccelResid)) +
  geom_line(aes(group = subject_id,
                color = subject_id),
            show.legend = FALSE) +
  geom_point(aes(group = subject_id,
                 color = subject_id),
             show.legend = FALSE) +
  theme_bw() +
  facet_wrap(facets = vars(subject_id))

plot

ggsave(plot,
       filename = "age_AccelResid_each_person.pdf",
       width = 7,
       height = 6)

library(ggsignif)

plot <-
  temp_data2 %>%
  # dplyr::filter(subject_id %in% depressed_subject_id) %>%
  dplyr::mutate(subject_id = factor(subject_id, levels = stringr::str_sort(unique(
    temp_data2$subject_id
  ), numeric = TRUE))) %>%
  dplyr::filter(!is.na(AgeAccelResid)) %>%
  ggplot(aes(time_point, AgeAccelResid)) +
  geom_boxplot(aes(group = time_point)) +
  geom_line(aes(group = subject_id,
                color = subject_id),
            show.legend = FALSE) +
  geom_point(aes(group = subject_id,
                 color = subject_id),
             show.legend = FALSE) +
  geom_signif(
    comparisons = list(c("T1", "T2")),
    map_signif_level = TRUE,
    y_position = 10.3,
    test = "wilcox.test"
  ) +
  geom_signif(
    comparisons = list(c("T1", "T3")),
    map_signif_level = TRUE,
    y_position = 12,
    test = "wilcox.test"
  ) +
  geom_signif(
    comparisons = list(c("T2", "T3")),
    map_signif_level = TRUE,
    y_position = 5,
    test = "wilcox.test"
  ) +
  theme_bw() +
  labs(x = "")

plot

ggsave(plot,
       filename = "age_AccelResid.pdf",
       width = 7,
       height = 6)


plot <-
  temp_data %>%
  # dplyr::filter(subject_id %in% depressed_subject_id) %>%
  dplyr::mutate(subject_id = factor(subject_id, levels = stringr::str_sort(unique(
    temp_data$subject_id
  ), numeric = TRUE))) %>%
  ggplot(aes(time_point, AgeAccelResid)) +
  geom_line(aes(group = subject_id,
                color = subject_id),
            show.legend = FALSE) +
  geom_point(aes(group = subject_id,
                 color = subject_id),
             show.legend = FALSE) +
  theme_bw() +
  labs(x = "") +
  theme(panel.grid.minor = element_blank()) +
  facet_wrap(facets = vars(subject_id))

ggsave(plot,
       filename = "age_AccelResid_for_each_participant.pdf",
       width = 7,
       height = 6)
