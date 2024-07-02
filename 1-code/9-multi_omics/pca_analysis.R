no_function

setwd(r4projects::get_project_wd())

rm(list = ls())

source("1-code/tools.R")

###load data

load("3-data_analysis/transcriptomics/data_preparation/transcriptomics_data")

transcriptomics_sample_info <-
  extract_sample_info(transcriptomics_data)

###cardopanel
load("3-data_analysis/Cardiovascular_Risk_Panel/data_preparation/expression_data")
load("3-data_analysis/Cardiovascular_Risk_Panel/data_preparation/sample_info")
load("3-data_analysis/Cardiovascular_Risk_Panel/data_preparation/variable_info")

load("3-data_analysis/Cardiovascular_Risk_Panel/marker/depression_association_pos")
load("3-data_analysis/Cardiovascular_Risk_Panel/marker/depression_association_neg")

expression_data[1, ] =
  scale(as.numeric(expression_data[1, ])) %>%
  as.numeric()

range(expression_data)

depression_association_pos$Variables
depression_association_neg$Variables

marker_name = c(depression_association_pos$Variables,
                depression_association_neg$Variables)

expression_data_cardopanel <-
  expression_data[marker_name,]

variable_info_cardopanel <-
  variable_info[match(marker_name, variable_info$variable_id), , drop = FALSE] %>%
  dplyr::mutate(class = "Cardiopanel")


###cytokine
load("3-data_analysis/Cytokines/data_preparation/expression_data")
load("3-data_analysis/Cytokines/data_preparation/sample_info")
load("3-data_analysis/Cytokines/data_preparation/variable_info")

load("3-data_analysis/Cytokines/marker/depression_association_pos")
load("3-data_analysis/Cytokines/marker/depression_association_neg")

depression_association_pos$Variables
depression_association_neg$Variables

marker_name = c(depression_association_pos$Variables,
                depression_association_neg$Variables)

marker_name

expression_data_cytokine <-
  expression_data[marker_name,]

variable_info_cytokine <-
  variable_info[match(marker_name, variable_info$variable_id), , drop = FALSE] %>%
  dplyr::mutate(class = "Cytokine")

range(expression_data_cytokine)

###lipid
load("3-data_analysis/Lipids/data_preparation/expression_data")
load("3-data_analysis/Lipids/data_preparation/sample_info")
load("3-data_analysis/Lipids/data_preparation/variable_info")

load("3-data_analysis/Lipids/marker/depression_association_pos")
load("3-data_analysis/Lipids/marker/depression_association_neg")

depression_association_pos$Variables
depression_association_neg$Variables

marker_name = c(depression_association_pos$Variables,
                depression_association_neg$Variables)

marker_name

expression_data_lipid <-
  expression_data[marker_name,]

variable_info_lipid <-
  variable_info[match(marker_name, variable_info$variable_id), , drop = FALSE] %>%
  dplyr::mutate(class = "Lipid")

range(expression_data_lipid)

###metabolic_panel
load("3-data_analysis/Metabolic_Panel/data_preparation/expression_data")
load("3-data_analysis/Metabolic_Panel/data_preparation/sample_info")
load("3-data_analysis/Metabolic_Panel/data_preparation/variable_info")

load("3-data_analysis/Metabolic_Panel/marker/depression_association_pos")
load("3-data_analysis/Metabolic_Panel/marker/depression_association_neg")

depression_association_pos$Variables
depression_association_neg$Variables

marker_name = c(depression_association_pos$Variables,
                depression_association_neg$Variables)

marker_name

expression_data_metabolic_panel <-
  expression_data[marker_name,]

variable_info_metabolic_panel <-
  variable_info[match(marker_name, variable_info$variable_id), , drop = FALSE] %>%
  dplyr::mutate(class = "Metabolic_panel")

range(expression_data_metabolic_panel)

###metabolomics
load("3-data_analysis/metabolomics_data/data_preparation/metabolites/expression_data")
load("3-data_analysis/metabolomics_data/data_preparation/metabolites/sample_info")
load("3-data_analysis/metabolomics_data/data_preparation/metabolites/variable_info")

load("3-data_analysis/metabolomics_data/marker/depression_association_pos")
load("3-data_analysis/metabolomics_data/marker/depression_association_neg")

depression_association_pos$Variables
depression_association_neg$Variables

marker_name = c(depression_association_pos$Variables,
                depression_association_neg$Variables)

marker_name

expression_data_metabolomics <-
  expression_data[marker_name,] %>%
  dplyr::filter(!is.na(`1_T1`))

variable_info_metabolomics <-
  variable_info[match(marker_name, variable_info$variable_id), , drop = FALSE] %>%
  dplyr::filter(!is.na(variable_id)) %>%
  dplyr::mutate(class = "Metabolomics")

rownames(expression_data_metabolomics) == variable_info_metabolomics$variable_id

range(expression_data_metabolomics)

setwd(r4projects::get_project_wd())
dir.create("3-data_analysis/multi_omics/pca_analysis/")
setwd("3-data_analysis/multi_omics/pca_analysis/")

intersect_name <-
  Reduce(f = intersect,
         x = list(
           colnames(expression_data_cardopanel),
           colnames(expression_data_cytokine),
           colnames(expression_data_lipid),
           colnames(expression_data_metabolic_panel),
           colnames(expression_data_metabolomics)
         ))
expression_data <-
  rbind(
    expression_data_cardopanel[, intersect_name],
    expression_data_cytokine[, intersect_name],
    expression_data_lipid[, intersect_name],
    expression_data_metabolic_panel[, intersect_name],
    expression_data_metabolomics[, intersect_name]
  )

variable_info <-
  variable_info_cardopanel %>%
  dplyr::full_join(variable_info_cytokine,
                   by = intersect(
                     colnames(variable_info_cardopanel),
                     colnames(variable_info_cytokine)
                   )) %>%
  dplyr::full_join(variable_info_lipid,
                   by = intersect(colnames(.),
                                  colnames(variable_info_lipid))) %>%
  dplyr::full_join(variable_info_metabolic_panel,
                   by = intersect(colnames(.),
                                  colnames(variable_info_metabolic_panel))) %>%
  dplyr::full_join(variable_info_metabolomics,
                   by = intersect(colnames(.),
                                  colnames(variable_info_metabolomics)))

variable_info$mol_name = variable_info$Metabolite
variable_info$mol_name[is.na(variable_info$mol_name)] =
  variable_info$variable_id[is.na(variable_info$mol_name)]

sample_info <-
  sample_info[match(colnames(expression_data), sample_info$sample_id), ]


dim(expression_data)
dim(sample_info)
dim(variable_info)

library(tidymass)

sample_info <-
  sample_info %>%  
  dplyr::mutate(class = "Subject") %>% 
  dplyr::left_join(transcriptomics_sample_info[,c("sample_id", "depressed")]) %>% 
  dplyr::filter(!is.na(depressed))

expression_data <-
  expression_data[,sample_info$sample_id]

temp_data <-
  create_mass_dataset(
    expression_data = expression_data,
    sample_info = sample_info,
    variable_info = variable_info
  )

plot <-
  massqc::massqc_pca(object = temp_data, color_by = "depressed")

plot


###only for T1
plot1 <-
  temp_data %>% 
  activate_mass_dataset(what = "sample_info") %>% 
  filter(Time == "T1") %>% 
  massqc::massqc_pca(color_by = "depressed")

plot1

ggsave(plot1,
       filename = "pca_t1.pdf",
       width = 9,
       height = 7)
ggsave(plot1,
       filename = "pca_t1.png",
       width = 9,
       height = 7)



###only for T2
plot2 <-
  temp_data %>% 
  activate_mass_dataset(what = "sample_info") %>% 
  filter(Time == "T2") %>% 
  massqc::massqc_pca(color_by = "depressed")

plot2

ggsave(plot2,
       filename = "pca_t2.pdf",
       width = 9,
       height = 7)
ggsave(plot2,
       filename = "pca_t2.png",
       width = 9,
       height = 7)



###only for T3
plot3 <-
  temp_data %>% 
  activate_mass_dataset(what = "sample_info") %>% 
  filter(Time == "T3") %>% 
  massqc::massqc_pca(color_by = "depressed")

plot3

ggsave(plot3,
       filename = "pca_t3.pdf",
       width = 9,
       height = 7)
ggsave(plot3,
       filename = "pca_t3.png",
       width = 9,
       height = 7)





###only for T4
plot4 <-
  temp_data %>% 
  activate_mass_dataset(what = "sample_info") %>% 
  filter(Time == "T4") %>% 
  massqc::massqc_pca(color_by = "depressed")

plot4

ggsave(plot4,
       filename = "pca_t4.pdf",
       width = 9,
       height = 7)
ggsave(plot4,
       filename = "pca_t4.png",
       width = 9,
       height = 7)



###only for T5
plot5 <-
  temp_data %>% 
  activate_mass_dataset(what = "sample_info") %>% 
  filter(Time == "T5") %>% 
  massqc::massqc_pca(color_by = "depressed")

plot5

ggsave(plot5,
       filename = "pca_t5.pdf",
       width = 9,
       height = 7)
ggsave(plot5,
       filename = "pca_t5.png",
       width = 9,
       height = 7)


