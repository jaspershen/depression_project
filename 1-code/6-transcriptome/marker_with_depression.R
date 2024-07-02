###
no_source()

setwd(r4projects::get_project_wd())
rm(list = ls())

load("3-data_analysis/transcriptomics/data_preparation/transcriptomics_data")

dir.create("3-data_analysis/transcriptomics/marker_depression")
setwd("3-data_analysis/transcriptomics/marker_depression")

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

transcriptomics_data <-
  transcriptomics_data %>%
  activate_mass_dataset(what = "sample_info") %>%
  dplyr::filter(!is.na(Time))

####remove some genes that have more than 50% zero values
transcriptomics_data <-
  transcriptomics_data %>%
  mutate_variable_zero_freq()

transcriptomics_data <-
  transcriptomics_data %>%
  activate_mass_dataset(what = "variable_info") %>%
  dplyr::filter(zero_freq < 0.5)

transcriptomics_data@sample_info$Time

expression_data <-
  transcriptomics_data@expression_data

sample_info <-
  transcriptomics_data@sample_info

variable_info <-
  transcriptomics_data@variable_info

new_expression_data <-
  sample_info$Time %>%
  unique() %>%
  sort() %>%
  purrr::map(function(x) {
    idx <-
      which(sample_info$Time == x)
    expression_data[, idx] %>%
      apply(1, mean)
  }) %>%
  do.call(cbind, .) %>%
  as.data.frame()

colnames(new_expression_data) <-
  sample_info$Time %>%
  unique() %>%
  sort()

new_sample_info <-
  sample_info %>%
  dplyr::distinct(Time, keep_all = TRUE) %>%
  dplyr::arrange(Time)

####remove some genes that have more than 50% zero values
zero_ratio <-
  new_expression_data %>%
  apply(1, function(x) {
    sum(x == 0) / ncol(new_expression_data)
  })

zero_ratio %>% plot()

lm_result <- 
seq_len(nrow(new_expression_data)) %>%
  purrr::map(function(idx) {
    temp_data <-
      data.frame(y = as.numeric(new_expression_data[idx,]),
                 x = 1:6)
    lm_reg <-
      lm(y ~ x, data = temp_data)
    summary_fit <- summary(lm_reg)
    coeffs <- summary_fit$coefficients[,1]
    p_values <- summary_fit$coefficients[,4]
    
    data.frame(intercept = unname(coeffs[1]),
               x = unname(coeffs[2]),
               p = unname(p_values[2]))
    
  }) %>% 
  do.call(rbind, .) %>% 
  as.data.frame()

lm_result$variable_id <-
  variable_info$variable_id

lm_result$p_adjust <- p.adjust(lm_result$p, method = "fdr")

save(lm_result, file = "lm_result")
  
