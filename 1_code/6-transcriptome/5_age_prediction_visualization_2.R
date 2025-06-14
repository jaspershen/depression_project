#####all the depressed participants with at least 2 samples
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

dir.create(
  "3_data_analysis/transcriptomics/age_prediction/all_depressed_participants",
  recursive = TRUE
)
setwd("3_data_analysis/transcriptomics/age_prediction/all_depressed_participants")

load("../result")

###summary
result <-
  result %>%
  tibble::rownames_to_column(var = "sample_id") %>%
  dplyr::mutate(subject_id = stringr::str_extract(sample_id, "[0-9]{1,2}")) %>%
  dplyr::mutate(time_point = stringr::str_extract(sample_id, "T[0-9]{1,2}")) %>%
  dplyr::arrange(subject_id, time_point)


result <-
  result %>%
  dplyr::mutate(age_diff = RNAAge - ChronAge) %>%
  dplyr::select(subject_id, age_diff, time_point, AgeAccelResid)

####AgeAccelResid
temp_data <-
  result %>%
  dplyr::filter(!is.na(AgeAccelResid))

###only remain the subject with at least two samples
temp_data <-
  temp_data %>%
  dplyr::filter(subject_id %in% unique(temp_data$subject_id[duplicated(temp_data$subject_id)]))

###only remain the depressed participants at T1
depressed_subject_id <-
  sample_info %>%
  dplyr::filter(Time == "T1") %>%
  dplyr::filter(bdi_total >= 14) %>%
  dplyr::pull(subject_id)

temp_data <-
  temp_data %>%
  dplyr::filter(subject_id %in% depressed_subject_id)

length(unique(temp_data$subject_id))

##21 subjects
plot <-
  temp_data %>%
  dplyr::mutate(subject_id = factor(subject_id, levels = stringr::str_sort(unique(
    temp_data$subject_id
  ), numeric = TRUE))) %>%
  dplyr::filter(!is.na(AgeAccelResid)) %>%
  ggplot(aes(time_point, AgeAccelResid)) +
  geom_line(aes(group = subject_id, color = subject_id), show.legend = FALSE) +
  geom_point(aes(group = subject_id, color = subject_id), show.legend = FALSE) +
  theme_bw() +
  facet_wrap(facets = vars(subject_id))

plot

ggsave(plot,
       filename = "age_AccelResid_each_person.pdf",
       width = 7,
       height = 6)

library(ggsignif)

plot <-
  temp_data %>%
  dplyr::filter(time_point %in% c("T1", "T2", "T3")) %>%
  dplyr::mutate(subject_id = factor(subject_id, levels = stringr::str_sort(unique(
    temp_data$subject_id
  ), numeric = TRUE))) %>%
  dplyr::filter(!is.na(AgeAccelResid)) %>%
  ggplot(aes(time_point, AgeAccelResid)) +
  geom_boxplot(aes(group = time_point)) +
  geom_line(aes(group = subject_id, color = subject_id), show.legend = FALSE) +
  geom_point(aes(group = subject_id, color = subject_id), show.legend = FALSE) +
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
       filename = "age_AccelResid_box_plot.pdf",
       width = 7,
       height = 6)

library(tidyr)

####wilxon test
t1_data <-
  temp_data %>%
  dplyr::filter(time_point == "T1") %>%
  dplyr::arrange(subject_id)

t2_data <-
  temp_data %>%
  dplyr::filter(time_point == "T2") %>%
  dplyr::arrange(subject_id)

t3_data <-
  temp_data %>%
  dplyr::filter(time_point == "T3") %>%
  dplyr::arrange(subject_id)

wilcox.test(
  t1_data %>%
    dplyr::filter(
      subject_id %in% intersect(t1_data$subject_id, t2_data$subject_id)
    ) %>%
    pull(AgeAccelResid),
  t2_data %>%
    dplyr::filter(
      subject_id %in% intersect(t1_data$subject_id, t2_data$subject_id)
    ) %>%
    pull(AgeAccelResid),
  paired = TRUE
)


wilcox.test(
  t1_data %>%
    dplyr::filter(
      subject_id %in% intersect(t1_data$subject_id, t3_data$subject_id)
    ) %>%
    pull(AgeAccelResid),
  t3_data %>%
    dplyr::filter(
      subject_id %in% intersect(t1_data$subject_id, t3_data$subject_id)
    ) %>%
    pull(AgeAccelResid),
  paired = TRUE
)


wilcox.test(
  t2_data %>%
    dplyr::filter(
      subject_id %in% intersect(t2_data$subject_id, t3_data$subject_id)
    ) %>%
    pull(AgeAccelResid),
  t3_data %>%
    dplyr::filter(
      subject_id %in% intersect(t2_data$subject_id, t3_data$subject_id)
    ) %>%
    pull(AgeAccelResid),
  paired = TRUE
)

