##
no_function()

setwd(r4projects::get_project_wd())
library(tidyverse)
rm(list = ls())
source("1-code/tools.R")

###load data
###cardopanel
load("3_data_analysis/Cardiovascular_Risk_Panel/data_preparation/expression_data")
load("3_data_analysis/Cardiovascular_Risk_Panel/data_preparation/sample_info")
load("3_data_analysis/Cardiovascular_Risk_Panel/data_preparation/variable_info")

load("3_data_analysis/Cardiovascular_Risk_Panel/marker/depression_association_pos")
load("3_data_analysis/Cardiovascular_Risk_Panel/marker/depression_association_neg")

expression_data[1,] =
  scale(as.numeric(expression_data[1,])) %>% 
  as.numeric()

range(expression_data)  

depression_association_pos$Variables
depression_association_neg$Variables

marker_name = c(depression_association_pos$Variables,
                depression_association_neg$Variables)

expression_data_cardopanel <-
  expression_data[marker_name, ]

variable_info_cardopanel <-
  variable_info[match(marker_name, variable_info$variable_id), , drop = FALSE] %>% 
  dplyr::mutate(class = "Cardiopanel")


expression_data_cardopanel %>%
  tibble::rownames_to_column(var = "variable_id") %>%
  tidyr::pivot_longer(cols = -variable_id,
                      names_to = "sample_id",
                      values_to = "value") %>%
  dplyr::left_join(sample_info, by = c("sample_id")) %>%
  dplyr::filter(!is.na(bdi_total)) %>%
  dplyr::mutate(subject_id = as.character(subject_id)) %>% 
  ggplot(aes(bdi_total, value)) +
  geom_smooth(aes(color = subject_id), se = FALSE, 
              show.legend = FALSE, method = "lm") +
  geom_point(aes(color = subject_id), 
             size = 2,
             shape = 16, show.legend = FALSE,
             alpha = 0.5) +
  facet_wrap(facets = ~ variable_id, scales = "free") +
  base_theme
