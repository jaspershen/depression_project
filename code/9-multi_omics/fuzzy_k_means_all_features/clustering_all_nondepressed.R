##
no_function()

setwd(r4projects::get_project_wd())
library(tidyverse)
rm(list = ls())
source("code/tools.R")

###load data
###cardopanel
load("data_analysis/Cardiovascular_Risk_Panel/data_preparation/expression_data")
load("data_analysis/Cardiovascular_Risk_Panel/data_preparation/sample_info")
load("data_analysis/Cardiovascular_Risk_Panel/data_preparation/variable_info")

variable_info <-
  variable_info %>%
  dplyr::mutate(mol_name = variable_id) %>%
  dplyr::mutate(variable_id = paste("cardopanel", 1:nrow(variable_info), sep = "_"))

rownames(expression_data) <-
  variable_info$variable_id

# load("data_analysis/Cardiovascular_Risk_Panel/marker/depression_association_pos")
# load("data_analysis/Cardiovascular_Risk_Panel/marker/depression_association_neg")

expression_data[1,] =
  scale(as.numeric(expression_data[1,])) %>%
  as.numeric()

range(expression_data)

# depression_association_pos$Variables
# depression_association_neg$Variables
#
# marker_name = c(depression_association_pos$Variables,
#                 depression_association_neg$Variables)

expression_data_cardopanel <-
  expression_data
# expression_data[marker_name, ]

variable_info_cardopanel <-
  variable_info %>%
  # variable_info[match(marker_name, variable_info$variable_id), , drop = FALSE] %>%
  dplyr::mutate(class = "Cardiopanel")

sample_info_cardopanel <-
  sample_info

###cytokine
load("data_analysis/Cytokines/data_preparation/expression_data")
load("data_analysis/Cytokines/data_preparation/sample_info")
load("data_analysis/Cytokines/data_preparation/variable_info")

load("data_analysis/Cytokines/marker/depression_association_pos")
load("data_analysis/Cytokines/marker/depression_association_neg")

variable_info <-
  variable_info %>%
  dplyr::mutate(mol_name = variable_id) %>%
  dplyr::mutate(variable_id = paste("cytokine", 1:nrow(variable_info), sep = "_"))

rownames(expression_data) <-
  variable_info$variable_id

depression_association_pos$Variables
depression_association_neg$Variables

marker_name = c(depression_association_pos$Variables,
                depression_association_neg$Variables)

marker_name

expression_data_cytokine <-
  expression_data
# expression_data[marker_name, ]

variable_info_cytokine <-
  variable_info %>%
  # variable_info[match(marker_name, variable_info$variable_id), ,drop = FALSE] %>%
  dplyr::mutate(class = "Cytokine")

range(expression_data_cytokine)

sample_info_cytokine <-
  sample_info

###lipid
load("data_analysis/Lipids/data_preparation/expression_data")
load("data_analysis/Lipids/data_preparation/sample_info")
load("data_analysis/Lipids/data_preparation/variable_info")

load("data_analysis/Lipids/marker/depression_association_pos")
load("data_analysis/Lipids/marker/depression_association_neg")

variable_info <-
  variable_info %>%
  dplyr::mutate(mol_name = variable_id) %>%
  dplyr::mutate(variable_id = paste("lipid", 1:nrow(variable_info), sep = "_"))

rownames(expression_data) <-
  variable_info$variable_id

depression_association_pos$Variables
depression_association_neg$Variables

marker_name = c(depression_association_pos$Variables,
                depression_association_neg$Variables)

marker_name

expression_data_lipid <-
  expression_data
# expression_data[marker_name, ]

variable_info_lipid <-
  variable_info %>%
  # variable_info[match(marker_name, variable_info$variable_id), ,drop = FALSE] %>%
  dplyr::mutate(class = "Lipid")

range(expression_data_lipid)

sample_info_lipid <-
  sample_info

###metabolic_panel
load("data_analysis/Metabolic_Panel/data_preparation/expression_data")
load("data_analysis/Metabolic_Panel/data_preparation/sample_info")
load("data_analysis/Metabolic_Panel/data_preparation/variable_info")

load("data_analysis/Metabolic_Panel/marker/depression_association_pos")
load("data_analysis/Metabolic_Panel/marker/depression_association_neg")

variable_info <-
  variable_info %>%
  dplyr::mutate(mol_name = variable_id) %>%
  dplyr::mutate(variable_id = paste("metabolic_panel", 1:nrow(variable_info), sep = "_"))

rownames(expression_data) <-
  variable_info$variable_id

depression_association_pos$Variables
depression_association_neg$Variables

marker_name = c(depression_association_pos$Variables,
                depression_association_neg$Variables)

marker_name

expression_data_metabolic_panel <-
  expression_data
# expression_data[marker_name, ]

variable_info_metabolic_panel <-
  variable_info %>%
  # variable_info[match(marker_name, variable_info$variable_id), ,drop = FALSE] %>%
  dplyr::mutate(class = "Metabolic_panel")

range(expression_data_metabolic_panel)

sample_info_metabolic_panel <-
  sample_info

###metabolomics
load("data_analysis/metabolomics_data/data_preparation/metabolites/expression_data")
load("data_analysis/metabolomics_data/data_preparation/metabolites/sample_info")
load("data_analysis/metabolomics_data/data_preparation/metabolites/variable_info")

load("data_analysis/metabolomics_data/marker/depression_association_pos")
load("data_analysis/metabolomics_data/marker/depression_association_neg")

depression_association_pos$Variables
depression_association_neg$Variables

marker_name = c(depression_association_pos$Variables,
                depression_association_neg$Variables)

marker_name

expression_data_metabolomics <-
  expression_data %>%
  # expression_data[marker_name, ] %>%
  dplyr::filter(!is.na(`1_T1`))

variable_info_metabolomics <-
  variable_info %>%
  # variable_info[match(marker_name, variable_info$variable_id), ,drop = FALSE] %>%
  dplyr::filter(!is.na(variable_id)) %>%
  dplyr::mutate(class = "Metabolomics")

rownames(expression_data_metabolomics) == variable_info_metabolomics$variable_id

range(expression_data_metabolomics)

sample_info_metabolomics <-
  sample_info

###transcriptomics
load("data_analysis/transcriptomics/data_preparation/transcriptomics_data")

dim(transcriptomics_data)

###only remain protein-coding
library(tidymass)

transcriptomics_data <-
  transcriptomics_data %>%
  activate_mass_dataset(what = "variable_info") %>%
  dplyr::filter(!is.na(GENETYPE)) %>%
  dplyr::filter(GENETYPE == "protein-coding")

colnames(transcriptomics_data)

###remove the genes < 50% samples
transcriptomics_data <-
  transcriptomics_data %>%
  mutate_variable_zero_freq(according_to_samples = "all") %>%
  activate_mass_dataset(what = "variable_info") %>%
  dplyr::filter(zero_freq < 0.5)

##remove samples without Time information
transcriptomics_data <-
  transcriptomics_data %>%
  activate_mass_dataset(what = "sample_info") %>%
  dplyr::filter(!is.na(Time))

expression_data_transcriptomics <-
  transcriptomics_data@expression_data

sample_info_transcriptomics <-
  extract_sample_info(transcriptomics_data)

variable_info_transcriptomics <-
  extract_variable_info(transcriptomics_data)

variable_info_transcriptomics <-
  variable_info_transcriptomics %>%
  dplyr::mutate(mol_name = variable_id) %>%
  dplyr::mutate(variable_id = paste(
    "transcriptomics",
    1:nrow(variable_info_transcriptomics),
    sep = "_"
  ))

rownames(expression_data_transcriptomics) <-
  variable_info_transcriptomics$variable_id


setwd(r4projects::get_project_wd())
dir.create("data_analysis/multi_omics_all_features/non-depressed/",
           recursive = TRUE)
setwd("data_analysis/multi_omics_all_features/non-depressed/")

sample_info_cardopanel <-
  sample_info_cardopanel %>%
  dplyr::left_join(sample_info_transcriptomics[, c("subject_id", "depressed")],
                   by = "subject_id") %>%
  dplyr::filter(!is.na(depressed)) %>%
  dplyr::filter(depressed == "Not Depressed")

expression_data_cardopanel <-
  expression_data_cardopanel[, sample_info_cardopanel$sample_id]


expression_data_cardopanel2 <-
  sort(unique(sample_info_cardopanel$Time)) %>%
  purrr::map(function(x) {
    idx <-
      which(sample_info_cardopanel$Time == x)
    expression_data_cardopanel[, idx] %>%
      apply(1, mean)
  }) %>%
  do.call(cbind, .) %>%
  as.data.frame()

colnames(expression_data_cardopanel2) <-
  sort(unique(sample_info_cardopanel$Time))


sample_info_cytokine <-
  sample_info_cytokine %>%
  dplyr::left_join(sample_info_transcriptomics[, c("subject_id", "depressed")],
                   by = "subject_id") %>%
  dplyr::filter(!is.na(depressed)) %>%
  dplyr::filter(depressed == "Not Depressed")

expression_data_cytokine <-
  expression_data_cytokine[, sample_info_cytokine$sample_id]

expression_data_cytokine2 <-
  sort(unique(sample_info_cytokine$Time)) %>%
  purrr::map(function(x) {
    idx <-
      which(sample_info_cytokine$Time == x)
    expression_data_cytokine[, idx] %>%
      apply(1, mean)
  }) %>%
  do.call(cbind, .) %>%
  as.data.frame()

colnames(expression_data_cytokine2) <-
  sort(unique(sample_info_cytokine$Time))


sample_info_lipid <-
  sample_info_lipid %>%
  dplyr::left_join(sample_info_transcriptomics[, c("subject_id", "depressed")],
                   by = "subject_id") %>%
  dplyr::filter(!is.na(depressed)) %>%
  dplyr::filter(depressed == "Not Depressed")

expression_data_lipid <-
  expression_data_lipid[, sample_info_lipid$sample_id]

expression_data_lipid2 <-
  sort(unique(sample_info_lipid$Time)) %>%
  purrr::map(function(x) {
    idx <-
      which(sample_info_lipid$Time == x)
    expression_data_lipid[, idx] %>%
      apply(1, mean)
  }) %>%
  do.call(cbind, .) %>%
  as.data.frame()

colnames(expression_data_lipid2) <-
  sort(unique(sample_info_lipid$Time))

sample_info_metabolic_panel <-
  sample_info_metabolic_panel %>%
  dplyr::filter(Time != "T4.5")

sample_info_metabolic_panel <-
  sample_info_metabolic_panel %>%
  dplyr::left_join(sample_info_transcriptomics[, c("subject_id", "depressed")],
                   by = "subject_id") %>%
  dplyr::filter(!is.na(depressed)) %>%
  dplyr::filter(depressed == "Not Depressed")

expression_data_metabolic_panel <-
  expression_data_metabolic_panel[, sample_info_metabolic_panel$sample_id]

expression_data_metabolic_panel2 <-
  sort(unique(sample_info_metabolic_panel$Time)) %>%
  purrr::map(function(x) {
    idx <-
      which(sample_info_metabolic_panel$Time == x)
    expression_data_metabolic_panel[, idx] %>%
      apply(1, mean)
  }) %>%
  do.call(cbind, .) %>%
  as.data.frame()

colnames(expression_data_metabolic_panel2) <-
  sort(unique(sample_info_metabolic_panel$Time))

sample_info_metabolomics <-
  sample_info_metabolomics %>%
  dplyr::left_join(sample_info_transcriptomics[, c("subject_id", "depressed")],
                   by = "subject_id") %>%
  dplyr::filter(!is.na(depressed)) %>%
  dplyr::filter(depressed == "Not Depressed")

expression_data_metabolomics <-
  expression_data_metabolomics[, sample_info_metabolomics$sample_id]

expression_data_metabolomics2 <-
  sort(unique(sample_info_metabolomics$Time)) %>%
  purrr::map(function(x) {
    idx <-
      which(sample_info_metabolomics$Time == x)
    expression_data_metabolomics[, idx] %>%
      apply(1, mean)
  }) %>%
  do.call(cbind, .) %>%
  as.data.frame()

colnames(expression_data_metabolomics2) <-
  sort(unique(sample_info_metabolomics$Time))


expression_data_transcriptomics <-
  expression_data_transcriptomics %>%
  apply(1, function(x) {
    (x - mean(x)) / sd(x)
  }) %>%
  t() %>%
  as.data.frame()

sample_info_transcriptomics <-
  sample_info_transcriptomics %>%
  dplyr::filter(Time != "T6")

sample_info_transcriptomics <-
  sample_info_transcriptomics %>%
  dplyr::filter(!is.na(depressed)) %>%
  dplyr::filter(depressed == "Not Depressed")

expression_data_transcriptomics <-
  expression_data_transcriptomics[, sample_info_transcriptomics$sample_id]

expression_data_transcriptomics2 <-
  sort(unique(sample_info_transcriptomics$Time)) %>%
  purrr::map(function(x) {
    idx <-
      which(sample_info_transcriptomics$Time == x)
    expression_data_transcriptomics[, idx] %>%
      apply(1, mean)
  }) %>%
  do.call(cbind, .) %>%
  as.data.frame()

colnames(expression_data_transcriptomics2) <-
  sort(unique(sample_info_transcriptomics$Time))

expression_data <-
  rbind(
    expression_data_cardopanel2,
    expression_data_cytokine2,
    expression_data_lipid2,
    expression_data_metabolic_panel2,
    expression_data_metabolomics2,
    expression_data_transcriptomics2
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
                                  colnames(variable_info_metabolomics))) %>%
  dplyr::full_join(variable_info_transcriptomics,
                   by = intersect(colnames(.),
                                  colnames(variable_info_transcriptomics)))

variable_info$variable_id == rownames(expression_data)
sum(variable_info$variable_id == rownames(expression_data))

variable_info$mol_name = variable_info$Metabolite
variable_info$mol_name[is.na(variable_info$mol_name)] =
  variable_info$variable_id[is.na(variable_info$mol_name)]

sample_info <-
  data.frame(sample_id = colnames(expression_data))

dim(expression_data)
dim(sample_info)
dim(variable_info)

library(Mfuzz)
library(e1071)

sample_info$sample_id == colnames(expression_data)
sample_info$Time

###clustering
library(Mfuzz)

temp_data <-
  expression_data

time <- colnames(temp_data)

temp_data <- rbind(time, temp_data)

row.names(temp_data)[1] <- "time"
rownames(temp_data)

# write.table(
#   temp_data,
#   file = "temp_data.txt",
#   sep = '\t',
#   quote = FALSE,
#   col.names = NA
# )

#read it back in as an expression set
data <- table2eset(filename = "temp_data.txt")
data.s <-
  # standardise(data)
  data
m1 <- mestimate(data.s)
m1

plot <-
  Dmin(
    data.s,
    m = m1,
    crange = seq(2, 40, 2),
    repeats = 3,
    visu = TRUE
  )

plot <-
  plot %>%
  data.frame(distance = plot,
             k = seq(2, 40, 2)) %>%
  ggplot(aes(k, distance)) +
  geom_point(shape = 21,
             size = 4,
             fill = "black") +
  # geom_smooth() +
  geom_segment(aes(
    x = k,
    y = 0,
    xend = k,
    yend = distance
  )) +
  theme_bw() +
  theme(
    # legend.position = c(0, 1),
    # legend.justification = c(0, 1),
    panel.grid = element_blank(),
    axis.title = element_text(size = 13),
    axis.text = element_text(size = 12),
    panel.background = element_rect(fill = "transparent", color = NA),
    plot.background = element_rect(fill = "transparent", color = NA),
    legend.background = element_rect(fill = "transparent", color = NA)
  ) +
  labs(x = "Cluster number",
       y = "Min. centroid distance") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1)))

plot

# ggsave(plot,
#        filename = "distance_k_number.pdf",
#        width = 7,
#        height = 7)

cluster_number <- 8

c <- mfuzz(data.s, c = cluster_number, m = m1)

save(c, file = "c")
load("c")
# ####any two clusters with correlation > 0.8 should be considered as one
library(corrplot)
layout(1)
# center <- c$centers

membership_cutoff <- 0.5

center <-
  get_mfuzz_center(data = data.s,
                   c = c,
                   membership_cutoff = 0.5)
rownames(center) <- paste("Cluster", rownames(center), sep = ' ')

corrplot::corrplot(
  corr = cor(t(center)),
  type = "full",
  diag = TRUE,
  order = "hclust",
  hclust.method = "ward.D",
  # addrect = 5,
  col = colorRampPalette(colors = rev(
    RColorBrewer::brewer.pal(n = 11, name = "Spectral")
  ))(n = 100),
  number.cex = .7,
  addCoef.col = "black"
)

mfuzz.plot(
  eset = data.s,
  min.mem = 0.5,
  cl = c,
  mfrow = c(3, 3),
  time.labels = time,
  new.window = FALSE
)

library(ComplexHeatmap)

cluster_info <-
  data.frame(
    variable_id = names(c$cluster),
    c$membership,
    cluster = c$cluster,
    stringsAsFactors = FALSE
  ) %>%
  arrange(cluster)

####plot for each cluster
idx <- 1

temp_data <-
  data.s %>% as.data.frame() %>%
  t() %>% as.data.frame()

for (idx in 1:cluster_number) {
  cat(idx, " ")
  
  cluster_data <-
    cluster_info %>%
    # dplyr::filter(cluster == idx) %>%
    dplyr::select(1, 1 + idx, cluster)
  
  colnames(cluster_data)[2] <- c("membership")
  
  cluster_data <-
    cluster_data %>%
    dplyr::filter(membership > membership_cutoff)
  
  path <- paste("cluster", idx, sep = "_")
  dir.create(path)
  
  openxlsx::write.xlsx(
    cluster_data,
    file = file.path(path, paste("cluster", idx, ".xlsx", sep = "")),
    asTable = TRUE,
    overwrite = TRUE
  )
  
  temp_center <-
    center[idx, , drop = TRUE] %>%
    unlist() %>%
    data.frame(time = names(.),
               value = .,
               stringsAsFactors = FALSE)
  
  temp <-
    temp_data[cluster_data$variable_id, ] %>%
    data.frame(
      membership = cluster_data$membership,
      .,
      stringsAsFactors = FALSE,
      check.names = FALSE
    ) %>%
    tibble::rownames_to_column(var = "variable_id") %>%
    tidyr::pivot_longer(
      cols = -c(variable_id, membership),
      names_to = "time",
      values_to = "value"
    )
  
  plot <-
    temp %>%
    dplyr::arrange(membership, variable_id) %>%
    dplyr::mutate(variable_id = factor(variable_id, levels = unique(variable_id))) %>%
    ggplot(aes(time, value, group = variable_id)) +
    geom_line(aes(color = membership), alpha = 0.7) +
    scale_x_discrete(expand = expansion(mult = c(0.05, 0.05))) +
    theme_bw() +
    theme(
      legend.position = "bottom",
      legend.justification = c(0, 1),
      panel.grid = element_blank(),
      axis.title = element_text(size = 13),
      axis.text = element_text(size = 12),
      axis.text.x = element_text(size = 12),
      panel.background = element_rect(fill = "transparent", color = NA),
      plot.background = element_rect(fill = "transparent", color = NA),
      legend.background = element_rect(fill = "transparent", color = NA)
    ) +
    labs(
      x = "",
      y = "Z-score",
      title = paste("Cluster ",
                    idx,
                    " (",
                    nrow(cluster_data),
                    " molecules)",
                    sep = "")
    ) +
    geom_line(
      mapping = aes(time, value, group = 1),
      data = temp_center,
      size = 2
    ) +
    geom_hline(yintercept = 0) +
    scale_color_continuous(type = "viridis")
  # scale_color_manual(values = omics_color)
  # viridis::scale_color_viridis()
  
  plot
  
  ggsave(
    plot,
    filename = file.path(path, paste("cluster", idx, ".pdf", sep = "")),
    width = 8,
    height = 7
  )
  ggsave(
    plot,
    filename = file.path(path, paste("cluster", idx, ".png", sep = "")),
    width = 8,
    height = 7
  )
}

table(cluster_info$cluster)

cluster_info <-
  unique(cluster_info$cluster) %>%
  purrr::map(function(x) {
    temp <-
      cluster_info %>%
      # dplyr::filter(cluster == x) %>%
      dplyr::select(variable_id, paste0("X", x), cluster)
    colnames(temp)[2] <- "membership"
    temp <-
      temp %>%
      dplyr::filter(membership >= membership_cutoff)
    temp <-
      temp %>%
      dplyr::mutate(cluster_raw = cluster) %>%
      dplyr::mutate(cluster = x)
    temp
  }) %>%
  dplyr::bind_rows() %>%
  as.data.frame()

cluster_info %>%
  dplyr::count(cluster)

cluster_info %>%
  dplyr::filter(membership > 0.5) %>%
  dplyr::count(cluster)

final_cluster_info <-
  cluster_info

save(final_cluster_info, file = "final_cluster_info")

openxlsx::write.xlsx(
  final_cluster_info,
  file = "final_cluster_info.xlsx",
  asTable = TRUE,
  overwrite = TRUE
)

dim(final_cluster_info)

dim(expression_data)

###heatmap
library(ComplexHeatmap)
library(circlize)
col_fun = colorRamp2(c(-1, 0, 1),
                     c(
                       viridis::viridis(n = 3)[1],
                       viridis::viridis(n = 3)[2],
                       viridis::viridis(n = 3)[3]
                     ))

temp_cluster_info <-
  final_cluster_info
# dplyr::filter(cluster %in% c("1", "4"))

library(ComplexHeatmap)

plot <-
  temp_data[temp_cluster_info$variable_id,] %>%
  Heatmap(
    show_row_names = FALSE,
    cluster_columns = FALSE,
    col = col_fun,
    border = TRUE,
    name = "Z-score",
    row_split = temp_cluster_info$cluster
  )

plot

plot <-
  ggplotify::as.ggplot(plot)
plot

ggsave(plot,
       file = "cluster_heatmap.pdf",
       width = 3,
       height = 7)

ggsave(plot,
       file = "cluster_heatmap.png",
       width = 3,
       height = 7)

