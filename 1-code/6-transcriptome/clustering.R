###
no_source()

setwd(r4projects::get_project_wd())
rm(list = ls())

setwd(r4projects::get_project_wd())

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

expression_data <-
  transcriptomics_data@expression_data

sample_info <-
  extract_sample_info(transcriptomics_data)

variable_info <-
  extract_variable_info(transcriptomics_data)

sample_info$Time

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

###clustering
library(Mfuzz)
temp_data <-
  new_expression_data

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
data.s <- standardise(data)
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

# save(c, file = "c")
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
    temp_data[cluster_data$variable_id,] %>%
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
col_fun = colorRamp2(c(-2, 0, 2),
                     c(
                       viridis::viridis(n = 3)[1],
                       viridis::viridis(n = 3)[2],
                       viridis::viridis(n = 3)[3]
                     ))

temp_cluster_info <-
  final_cluster_info %>%
  dplyr::filter(cluster %in% c("1", "4"))

library(ComplexHeatmap)

plot <-
  temp_data[temp_cluster_info$variable_id, ] %>%
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
