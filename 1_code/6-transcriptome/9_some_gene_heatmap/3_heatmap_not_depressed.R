###
no_source()

setwd(r4projects::get_project_wd())
rm(list = ls())

library(tidyverse)

load("3_data_analysis/transcriptomics/data_preparation/transcriptomics_data")

dir.create(
  "3_data_analysis/transcriptomics/heatmap_for_some_genes/non_depressed",
  recursive = TRUE
)
setwd("3_data_analysis/transcriptomics/heatmap_for_some_genes/non_depressed")

# data <- readr::read_csv("Gene.csv")

###19 pro-inflammatory indicator genes
gene_list1 <-
  c(
    "CXCL8",
    "FOS",
    "FOSB",
    "FOSL1",
    "FOSL2",
    "IL1A",
    "IL1B",
    "IL6",
    "JUN",
    "JUNB",
    "JUND",
    "NFKB1",
    "NFKB2",
    "PTGS1",
    "PTGS2",
    "REL",
    "RELA",
    "RELB",
    "TNF"
  )


###28 Type I interferon innate antiviral response genes as
gene_list2 <-
  c(
    "GBP1",
    "IFI16",
    "IFI27",
    "IFI27L1",
    "IFI27L2",
    "IFI30",
    "IFI35",
    "IFI44",
    "IFI44L",
    "IFI6",
    "IFIH1",
    "IFIT1",
    "IFIT2",
    "IFIT3",
    "IFIT5",
    "IFITM1",
    "IFITM2",
    "IFITM3",
    "IRF2",
    "IRF7",
    "IRF8",
    "JCHAIN",
    "MX1",
    "MX2",
    "OAS1",
    "OAS2",
    "OAS3",
    "OASL"
  )


##remove T6

transcriptomics_data <-
  transcriptomics_data %>%
  activate_mass_dataset(what = "sample_info") %>%
  dplyr::filter(!is.na(Time)) %>%
  dplyr::filter(!Time %in% "T6") %>%
  dplyr::filter(!is.na(depressed)) %>%
  dplyr::filter(depressed == "Not Depressed")

expression_data <-
  transcriptomics_data %>%
  extract_expression_data()

sample_info <-
  transcriptomics_data %>%
  extract_sample_info()

variable_info <-
  transcriptomics_data %>%
  extract_variable_info()

#####
new_expression_data <-
  unique(sample_info$Time) %>%
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
  unique(sample_info$Time) %>%
  sort()

new_expression_data <-
  new_expression_data %>%
  apply(1, function(x) {
    x <-
      (x - mean(x)) / sd(x)
    x <- x - x[1]
    x
  }) %>%
  t() %>%
  as.data.frame()

temp_data1 <-
  new_expression_data[gene_list1, ]

library(ComplexHeatmap)

col_fun <-
  circlize::colorRamp2(
    breaks = seq(-4, 4, length.out = 11),
    colors = rev(RColorBrewer::brewer.pal(n = 11, name = "RdBu"))
  )

plot1 <-
  Heatmap(
    temp_data1,
    cluster_columns = FALSE,
    col = col_fun,
    border = TRUE,
    name = "Z-score"
  )

library(ggplotify)
plot1 <-
  as.ggplot(plot1)
plot1

ggsave(plot1,
       filename = "gene_list1_heatmap.pdf",
       width = 4,
       height = 7)

ggsave(plot1,
       filename = "gene_list1_heatmap.png",
       width = 4,
       height = 7)

temp_data2 <-
  new_expression_data[gene_list2, ]

library(ComplexHeatmap)

col_fun <-
  circlize::colorRamp2(
    breaks = seq(-4, 4, length.out = 11),
    colors = rev(RColorBrewer::brewer.pal(n = 11, name = "RdBu"))
  )

plot2 <-
  Heatmap(
    temp_data2,
    cluster_columns = FALSE,
    col = col_fun,
    border = TRUE,
    name = "Z-score"
  )

library(ggplotify)
plot2 <-
  as.ggplot(plot2)
plot2

ggsave(plot2,
       filename = "gene_list2_heatmap.pdf",
       width = 4,
       height = 7)

ggsave(plot2,
       filename = "gene_list2_heatmap.png",
       width = 4,
       height = 7)

####Heatmap for all the participants
dim(expression_data)
temp_data1 <-
  expression_data[gene_list1,]

temp_data2 <-
  expression_data[gene_list2,]

sample_info <-
  data.frame(sample_id = colnames(temp_data1)) %>%
  dplyr::mutate(
    subject_id = as.numeric(stringr::str_extract(colnames(temp_data1), "[0-9]{1,3}")),
    time_point = stringr::str_extract(colnames(temp_data1), "T[0-9]{1,3}")
  ) %>%
  dplyr::arrange(time_point, subject_id) %>%
  dplyr::mutate(subject_id = as.character(subject_id))


##remove the subjects without T1 or only with one time point
library(plyr)
remove_subject_id <-
  sample_info %>%
  plyr::dlply(.variables = .(subject_id)) %>%
  purrr::map(function(x) {
    if (nrow(x) == 1) {
      return(x$subject_id[1])
    }
    if (all(x$time_point != "T1")) {
      return(x$subject_id[1])
    }
    return(NULL)
  }) %>%
  unlist()


if (length(remove_subject_id) > 0) {
  sample_info <-
    sample_info %>%
    dplyr::filter(!subject_id %in% remove_subject_id)
}

sample_info <-
  sample_info %>%
  dplyr::arrange(subject_id, time_point) %>%
  dplyr::mutate(subject_id = factor(subject_id, unique(subject_id)))


temp_data1 <-
  temp_data1[, sample_info$sample_id]

temp_data2 <-
  temp_data2[, sample_info$sample_id]

temp_data1 <-
  temp_data1 %>%
  apply(1, function(x) {
    (x - mean(x)) / sd(x)
  }) %>%
  t() %>%
  as.data.frame()

temp_data2 <-
  temp_data2 %>%
  apply(1, function(x) {
    (x - mean(x)) / sd(x)
  }) %>%
  t() %>%
  as.data.frame()

###normalize T1 for each participant
temp_data1 <-
  unique(sample_info$subject_id) %>%
  as.character() %>%
  purrr::map(function(id) {
    temp_sample_info <-
      sample_info %>%
      dplyr::filter(subject_id == id) %>%
      dplyr::arrange(time_point)
    temp_data1[, temp_sample_info$sample_id, drop = FALSE] %>%
      apply(1, function(x) {
        x <- x - x[1]
        x
      }) %>%
      t() %>%
      as.data.frame()
  }) %>%
  do.call(cbind, .) %>%
  as.data.frame()


temp_data1 <-
  temp_data1[, sample_info$sample_id]

###normalize T1 for each participant
temp_data2 <-
  unique(sample_info$subject_id) %>%
  as.character() %>%
  purrr::map(function(id) {
    temp_sample_info <-
      sample_info %>%
      dplyr::filter(subject_id == id) %>%
      dplyr::arrange(time_point)
    temp_data2[, temp_sample_info$sample_id, drop = FALSE] %>%
      apply(1, function(x) {
        x <- x - x[1]
        x
      }) %>%
      t() %>%
      as.data.frame()
  }) %>%
  do.call(cbind, .) %>%
  as.data.frame()

temp_data2 <-
  temp_data2[, sample_info$sample_id]

library(ComplexHeatmap)

column_split <-
  sample_info$time_point

subject_color <-
  colorRampPalette(colors = ggsci::pal_lancet()(n = 9))(length(levels(sample_info$subject_id)))

names(subject_color) <- levels(sample_info$subject_id)

column_anno <-
  HeatmapAnnotation(subject_id = sample_info$subject_id,
                    col = list(subject_id = subject_color))

library(circlize)

col_fun <-
  circlize::colorRamp2(
    breaks = seq(-4, 4, length.out = 11),
    colors = rev(RColorBrewer::brewer.pal(n = 11, name = "RdBu"))
  )

plot1 <-
  Heatmap(
    temp_data1,
    col = col_fun,
    cluster_columns = FALSE,
    cluster_rows = TRUE,
    show_column_names = FALSE,
    border = TRUE,
    column_split = column_split,
    name = "Z-score",
    top_annotation = column_anno
  )

plot1

plot1 <-
  ggplotify::as.ggplot(plot1)

ggsave(plot1,
       filename = "gene_list1_all_heatmap.pdf",
       width = 7,
       height = 7)

ggsave(plot1,
       filename = "gene_list1_all_heatmap.png",
       width = 7,
       height = 7)


plot2 <-
  Heatmap(
    temp_data2,
    col = col_fun,
    cluster_columns = FALSE,
    cluster_rows = TRUE,
    show_column_names = FALSE,
    border = TRUE,
    column_split = column_split,
    name = "Z-score",
    top_annotation = column_anno
  )

plot2

plot2 <-
  ggplotify::as.ggplot(plot2)
plot2
ggsave(plot2,
       filename = "gene_list2_all_heatmap.pdf",
       width = 7,
       height = 7)

ggsave(plot2,
       filename = "gene_list2_all_heatmap.png",
       width = 7,
       height = 7)

