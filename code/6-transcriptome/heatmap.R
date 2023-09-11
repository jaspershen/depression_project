###
no_source()

setwd(masstools::get_project_wd())
rm(list = ls())

library(tidyverse)

load("data_analysis/transcriptomics/data_preparation/transcriptomics_data")

dir.create("data_analysis/transcriptomics/heatmap_for_some_genes/",
           recursive = TRUE)
setwd("data_analysis/transcriptomics/heatmap_for_some_genes/")

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




#
#
#
#
#
#
# sample_info <-
#   data.frame(sample_id = colnames(temp_data1)) %>%
#   dplyr::mutate(
#     subject_id = as.numeric(stringr::str_extract(colnames(temp_data1), "[0-9]{1,3}")),
#     time_point = stringr::str_extract(colnames(temp_data1), "T[0-9]{1,3}")
#   ) %>%
#   dplyr::arrange(time_point, subject_id) %>%
#   dplyr::mutate(subject_id = as.character(subject_id)) %>%
#   dplyr::mutate(subject_id = factor(subject_id, levels = unique(subject_id)))
#
# temp_data1 <-
#   temp_data1[, sample_info$sample_id]
#
# temp_data1 <-
#   temp_data1 %>%
#   apply(1, function(x) {
#     (x - mean(x)) / sd(x)
#   }) %>%
#   t() %>%
#   as.data.frame()
#
# library(ComplexHeatmap)
#
# column_split <-
#   sample_info$time_point
#
# subject_color <-
#   colorRampPalette(colors = ggsci::pal_lancet()(n = 9))(length(levels(sample_info$subject_id)))
# names(subject_color) <- levels(sample_info$subject_id)
#
# column_anno <-
#   HeatmapAnnotation(subject_id = sample_info$subject_id,
#                     col = list(subject_id = subject_color))
#
# library(circlize)
#
# col_fun <-
#   circlize::colorRamp2(
#     breaks = seq(-4, 4, length.out = 11),
#     colors = rev(RColorBrewer::brewer.pal(n = 11, name = "RdBu"))
#   )
#
# Heatmap(
#   temp_data1,
#   col = col_fun,
#   cluster_columns = FALSE,
#   cluster_rows = FALSE,
#   show_column_names = FALSE,
#   border = TRUE,
#   column_split = column_split,
#   name = "Z-score",
#   top_annotation = column_anno
# )
#
# temp_data <-
#   temp_data1 %>%
#   tibble::rownames_to_column(var = "variable_id") %>%
#   tidyr::pivot_longer(cols = -variable_id,
#                       names_to = "sample_id",
#                       values_to = "value") %>%
#   tidyr::separate(
#     col = "sample_id",
#     into = c("subject_id", "time_point"),
#     sep = "-"
#   )
#
# ##remove some subjects
# remain_subject_id <-
#   temp_data %>%
#   dplyr::filter(time_point == "T1") %>%
#   pull(subject_id) %>%
#   unique()
#
# remove_subject_id <-
#   temp_data %>%
#   dplyr::distinct(subject_id, time_point, .keep_all = FALSE) %>%
#   dplyr::count(subject_id) %>%
#   dplyr::filter(n == 1) %>%
#   pull(subject_id)
#
# temp_data <-
#   temp_data %>%
#   dplyr::filter(subject_id %in% c(remain_subject_id)) %>%
#   dplyr::filter(!subject_id %in% c(remove_subject_id)) %>%
#   dplyr::mutate(subject_id = factor(subject_id, labels = stringr::str_sort(unique(subject_id), numeric = TRUE)))
#
# all_plot <-
#   purrr::map(stringr::str_sort(unique(temp_data$subject_id), numeric = TRUE), function(x) {
#     temp <-
#       temp_data %>%
#       dplyr::filter(subject_id == x)
#     left_time_point <-
#       setdiff(paste0("T", 1:6),
#               unique(temp$time_point))
#     if (length(left_time_point) > 0) {
#       left_data <-
#         purrr::map(left_time_point, function(y) {
#           data.frame(
#             variable_id = unique(temp$variable_id),
#             subject_id = unique(temp$subject_id),
#             time_point = y,
#             value = NA
#           )
#         }) %>%
#         do.call(rbind, .) %>%
#         as.data.frame()
#       temp <-
#         rbind(temp, left_data)
#     }
#
#     temp$value[which(temp$value > 3)] <- 3
#     temp$value[which(temp$value < -3)] <- -3
#
#     library(plyr)
#     temp <-
#       temp %>%
#       plyr::dlply(.variables = .(subject_id, variable_id)) %>%
#       purrr::map(function(z) {
#         z$value <-
#           z$value - z$value[z$time_point == "T1"]
#         z
#       }) %>%
#       do.call(rbind, .) %>%
#       as.data.frame()
#
#     temp %>%
#       ggplot(aes(time_point, variable_id)) +
#       geom_tile(aes(fill = value),
#                 color = "black",
#                 show.legend = TRUE) +
#       facet_wrap(facets = vars(subject_id)) +
#       scale_fill_gradient2(
#         low = "blue",
#         mid = "white",
#         high = "red",
#         midpoint = 0,
#         na.value = "grey"
#       ) +
#       theme_bw() +
#       theme(panel.grid = element_blank(), plot.margin = margin(0, 0, 0, 0)) +
#       scale_x_discrete(expand = expansion(mult = c(0, 0))) +
#       scale_y_discrete(expand = expansion(mult = c(0, 0))) +
#       labs(x = "", y = "")
#   })
#
# library(patchwork)
# names(all_plot) <-
#   stringr::str_sort(unique(temp_data$subject_id), numeric = TRUE)
#
# length(all_plot)
#
# dir.create("for each participant")
# 1:length(all_plot) %>%
#   purrr::map(function(idx) {
#     plot <-
#       all_plot[[idx]]
#     ggsave(
#       plot,
#       filename = file.path("for each participant", paste0(names(all_plot)[idx], ".pdf")),
#       width = 7,
#       height = 5
#     )
#   })
