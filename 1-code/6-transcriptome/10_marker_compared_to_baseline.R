###
no_source()

setwd(r4projects::get_project_wd())
rm(list = ls())

load("3-data_analysis/transcriptomics/data_preparation/transcriptomics_data")

dir.create("3-data_analysis/transcriptomics/marker_compared_to_baseline")
setwd("3-data_analysis/transcriptomics/marker_compared_to_baseline")

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


###only remain the depressed people
transcriptomics_data <-
  transcriptomics_data %>%
  activate_mass_dataset(what = "sample_info") %>%
  dplyr::filter(depressed == "Depressed")

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


#find all the genes in different time points
subject_data <-
  sample_info$Time %>% 
  sort() %>% 
  unique() %>% 
  purrr::map(function(x){
    expression_data[,which(sample_info$Time == x)]
  })

names(subject_data) <-
  sample_info$Time %>% 
  sort() %>% 
  unique()

table(sample_info$Time)

# fc_p_value <-
#   pbapply::pblapply(subject_data[-1], function(x) {
#     base_data <-
#       subject_data[[1]]
#     
#     compare_data <-
#       x
#     
#     sample_info1 <-
#     sample_info %>% 
#       dplyr::filter(sample_id %in% colnames(base_data)) %>% 
#       dplyr::arrange(subject_id)
#     
#     sample_info2 <-
#       sample_info %>% 
#       dplyr::filter(sample_id %in% colnames(x)) %>% 
#       dplyr::arrange(subject_id)
#     
#     intersect_subject_id <- 
#     intersect(sample_info1$subject_id,
#               sample_info2$subject_id)
#     
#     sample_info1 <-
#       sample_info1 %>% 
#       dplyr::filter(subject_id %in% intersect_subject_id)
#     
#     sample_info2 <-
#       sample_info2 %>% 
#       dplyr::filter(subject_id %in% intersect_subject_id)
#     
#     
#     base_data <-
#       base_data[,sample_info1$sample_id]
#     
#     x <-
#       x[,sample_info2$sample_id]
#     
#     p_value <- lapply(1:nrow(x), function(idx) {
#       wilcox.test(as.numeric(x[idx, ]), as.numeric(base_data[idx, ]))$p.value
#     }) %>%
#       unlist()
# 
#     fdr <- p.adjust(p_value, method = "fdr")
# 
#     fc <- lapply(1:nrow(x), function(idx) {
#       mean(as.numeric(x[idx, ])) / mean(as.numeric(base_data[idx, ]))
#     }) %>%
#       unlist()
# 
#     fc[is.infinite(fc)] <- max(fc[!is.infinite(fc)])
# 
#     data.frame(
#       variable_id = variable_info$variable_id,
#       p_value,
#       fdr,
#       fc,
#       stringsAsFactors = FALSE
#     )
#   })
# 
# save(fc_p_value, file = "fc_p_value")

load("fc_p_value")

# ##find marker for each time points
# marker_each_point <-
#   lapply(fc_p_value, function(x) {
#     idx1 <- which(x$p_value < 0.05 & x$fc > 1)
#     idx2 <- which(x$p_value < 0.05 & x$fc < 1)
# 
#     gene1 <-
#       try(data.frame(x[idx1, ],
#                      class = "increase",
#                      stringsAsFactors = FALSE),
#           silent = TRUE)
# 
#     if (class(gene1) == "try-error") {
#       gene1 <- NULL
#     }
# 
#     gene2 <-
#       try(data.frame(x[idx2, ],
#                      class = "decrease",
#                      stringsAsFactors = FALSE),
#           silent = TRUE)
# 
#     if (class(gene2) == "try-error") {
#       gene2 <- NULL
#     }
# 
#     rbind(gene1, gene2)
#   })
# 
# 
# save(marker_each_point, file = "marker_each_point")
load("marker_each_point")

#####a sankey
marker_each_point %>%
  lapply(
    FUN = function(x) {
      if (is.null(x)) {
        return(0)
      } else{
        return(nrow(x))
      }
    }
  ) %>%
  unlist()

all_marker_name <-
  lapply(marker_each_point, function(x) {
    x$variable_id
  }) %>%
  unlist() %>%
  unique()

length(all_marker_name)

library(ggalluvial)

temp_data <-
  lapply(marker_each_point, function(x) {
    if (is.null(x)) {
      x <-
        data.frame(
          variable_id = all_marker_name,
          class = "no",
          freq = 1,
          stringsAsFactors = FALSE
        )
    }
    x <-
      data.frame(variable_id = all_marker_name,
                 stringsAsFactors = FALSE) %>%
      left_join(x, by = "variable_id") %>%
      dplyr::select(variable_id, class)
    
    x$class[is.na(x$class)] <- "no"
    x$freq <- 1
    x
  })

names(temp_data) <-
  names(subject_data)[-1]
  
temp_data <-
  purrr::map2(
    .x = temp_data,
    .y = names(temp_data),
    .f = function(x, y) {
      if (is.null(x)) {
        return(NULL)
      }
      data.frame(x, point = y, stringsAsFactors = FALSE)
    }
  )

temp_data <-
  do.call(rbind, temp_data)

temp_data$point <-
  factor(temp_data$point, levels = unique(temp_data$point))

RColorBrewer::display.brewer.all()

plot1 <-
  ggplot(
    temp_data,
    aes(
      x = point,
      y = freq,
      stratum = class,
      alluvium = variable_id,
      fill = class,
      label = class
    )
  ) +
  scale_x_discrete(expand = c(0, 0)) +
  ggalluvial::geom_flow() +
  labs(x = "", y = "") +
  scale_fill_manual(
    values = c(
      "increase" = ggsci::pal_aaas()(10)[2],
      "decrease" = ggsci::pal_aaas()(10)[1],
      "no" = "azure2"
    )
  ) +
  ggalluvial::geom_stratum(alpha = 1, color = NA) +
  # geom_text(stat = "stratum", size = 3) +
  theme_bw() +
  theme(
    legend.position = "top",
    # panel.border = element_blank(),
    panel.grid = element_blank(),
    axis.ticks.x = element_blank(),
    # axis.text.y = element_blank(),
    axis.text.x = element_text(
      angle = 45,
      hjust = 1,
      vjust = 2
    ),
    panel.background = element_rect(fill = "transparent", color = NA),
    plot.background = element_rect(fill = "transparent", color = NA),
    legend.background = element_rect(fill = "transparent", color = NA)
  ) +
  scale_y_continuous(
    position = "right",
    expand = expansion(mult = c(0, 0)),
    labels = scales::scientific
  )

plot1

# ggsave(
#   plot1,
#   file = file.path("gene_sankey_light.pdf"),
#   width = 7,
#   height = 3,
#   bg = "transparent"
# )
