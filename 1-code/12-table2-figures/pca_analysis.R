no_function

setwd(r4projects::get_project_wd())

rm(list = ls())

data <-
  readxl::read_xlsx("3-data_analysis/some_figures/table2/data1.xlsx")

###metabolomics
load("3-data_analysis/metabolomics_data/data_preparation/metabolites/expression_data")
load("3-data_analysis/metabolomics_data/data_preparation/metabolites/sample_info")
load("3-data_analysis/metabolomics_data/data_preparation/metabolites/variable_info")

load("3-data_analysis/metabolomics_data/marker/depression_association_pos")
load("3-data_analysis/metabolomics_data/marker/depression_association_neg")

load("3-data_analysis/transcriptomics/data_preparation/transcriptomics_data")

transcriptomics_sample_info <-
  extract_sample_info(transcriptomics_data)

sample_info <-
  sample_info %>%
  dplyr::left_join(transcriptomics_sample_info[, c("subject_id", "depressed")] %>%
                     dplyr::distinct(subject_id, .keep_all = TRUE),
                   by = "subject_id") %>%
  dplyr::filter(!is.na(depressed))

dir.create("3-data_analysis/some_figures/table2/pca_analysis")
setwd("3-data_analysis/some_figures/table2/pca_analysis")

variable_info <-
  variable_info %>%
  dplyr::filter(Metabolite %in% data$Metabolite)

expression_data <-
  expression_data[variable_info$variable_id, sample_info$sample_id]

###pca_analysis
###only for T1
temp_sample_info <-
  sample_info %>%
  dplyr::filter(Time == "T1") %>%
  dplyr::mutate(class = "Subject")

temp_expression_data <-
  expression_data[, temp_sample_info$sample_id]

library(tidymass)

temp_data <-
  create_mass_dataset(
    expression_data = temp_expression_data,
    sample_info = temp_sample_info,
    variable_info = variable_info
  )

plot1 <-
massqc::massqc_pca(object = temp_data, color_by = "depressed")



###only for T2
temp_sample_info <-
  sample_info %>%
  dplyr::filter(Time == "T2") %>%
  dplyr::mutate(class = "Subject")

temp_expression_data <-
  expression_data[, temp_sample_info$sample_id]

library(tidymass)

temp_data <-
  create_mass_dataset(
    expression_data = temp_expression_data,
    sample_info = temp_sample_info,
    variable_info = variable_info
  )

plot2 <-
  massqc::massqc_pca(object = temp_data, color_by = "depressed")


###only for T3
temp_sample_info <-
  sample_info %>%
  dplyr::filter(Time == "T3") %>%
  dplyr::mutate(class = "Subject")

temp_expression_data <-
  expression_data[, temp_sample_info$sample_id]

library(tidymass)

temp_data <-
  create_mass_dataset(
    expression_data = temp_expression_data,
    sample_info = temp_sample_info,
    variable_info = variable_info
  )

plot3 <-
  massqc::massqc_pca(object = temp_data, color_by = "depressed")



###only for T4
temp_sample_info <-
  sample_info %>%
  dplyr::filter(Time == "T4") %>%
  dplyr::mutate(class = "Subject")

temp_expression_data <-
  expression_data[, temp_sample_info$sample_id]

library(tidymass)

temp_data <-
  create_mass_dataset(
    expression_data = temp_expression_data,
    sample_info = temp_sample_info,
    variable_info = variable_info
  )

plot4 <-
  massqc::massqc_pca(object = temp_data, color_by = "depressed")


###only for T5
temp_sample_info <-
  sample_info %>%
  dplyr::filter(Time == "T5") %>%
  dplyr::mutate(class = "Subject")

temp_expression_data <-
  expression_data[, temp_sample_info$sample_id]

library(tidymass)

temp_data <-
  create_mass_dataset(
    expression_data = temp_expression_data,
    sample_info = temp_sample_info,
    variable_info = variable_info
  )

plot5 <-
  massqc::massqc_pca(object = temp_data, color_by = "depressed")

ggsave(plot1, filename = "pca_t1.pdf", width = 9, height = 7)
ggsave(plot2, filename = "pca_t2.pdf", width = 9, height = 7)
ggsave(plot3, filename = "pca_t3.pdf", width = 9, height = 7)
ggsave(plot4, filename = "pca_t4.pdf", width = 9, height = 7)
ggsave(plot5, filename = "pca_t5.pdf", width = 9, height = 7)

ggsave(plot1, filename = "pca_t1.png", width = 9, height = 7)
ggsave(plot2, filename = "pca_t2.png", width = 9, height = 7)
ggsave(plot3, filename = "pca_t3.png", width = 9, height = 7)
ggsave(plot4, filename = "pca_t4.png", width = 9, height = 7)
ggsave(plot5, filename = "pca_t5.png", width = 9, height = 7)
