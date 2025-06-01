###no_function
setwd(r4projects::get_project_wd())
rm(list = ls())

source("1_code/tools.R")

setwd("3_data_analysis/Figure_1/")

library(tidyverse)

ibsr_data <- read.csv('Psychometric Data.csv')
# remove 66, 67, and 70-72
ibsr_data <- filter(ibsr_data, !(id %in% c(66, 67, 70, 71, 72)))

# Make a depression plot.
# Color by initial depression.

# Make Time into a factor
ibsr_data$Time <- as.factor(ibsr_data$Time)
# add a column based on the BDI_total
ibsr_data$Status <-
  as.factor(ifelse(ibsr_data$bdi_total < 14, "Not Depressed", "Depressed"))

# Also add a column based on INITIAL depression status
initial_status <- function(ID) {
  initial_dep <-
    filter(ibsr_data, id == ID, Time == 'T1')$bdi_total[1]
  return(ifelse(
    initial_dep < 14,
    "Not Initially Depressed",
    "Initially Depressed"
  ))
}

ibsr_data$`Initial Depression` <-
  factor(sapply(ibsr_data$id, initial_status))

# filter out the data for people with unknown initial depression status
ibsr_data <- filter(ibsr_data, !(is.na(`Initial Depression`)))

# added outlier.shape = NA because otherwise it was giving two points per outlier

plot <-
  ggplot(ibsr_data, aes(x = Time, y = bdi_total)) +
  geom_jitter(
    aes(fill = `Initial Depression`),
    size = 4,
    shape = 21,
    na.rm = TRUE
  ) +
  geom_boxplot(na.rm = TRUE,
               fill = "transparent",
               outlier.shape = NA) +
  labs(y = "Depression Score (BDI-II)", x = "") +
  theme_bw() +
  scale_fill_manual(values = initial_depressed_color) +
  base_theme +
  geom_hline(yintercept = 14,
             colour = "black",
             linetype = "dotted") +
  scale_x_discrete(
    labels = c(
      "T1" = "Baseline",
      "T2" = "Mid-Retreat",
      "T3" = "End of Retreat",
      "T4" = "One Month",
      "T5" = "Three Months",
      "T6" = "Six Months",
      "T7" = "One Year"
    )
  ) +
  theme(
    legend.position = c(1, 1),
    legend.justification = c(1, 1),
    axis.text = element_text(size = 15),
    axis.title = element_text(size = 18)
  )
plot
ggsave('Depression Boxplot.pdf', height = 6, width = 10)

ggsave('Depression Boxplot.png', height = 4.5, width = 10)


plot <-
  ggplot(ibsr_data, aes(x = Time, y = bdi_total)) +
  geom_jitter(
    aes(fill = `Initial Depression`),
    size = 4,
    shape = 21,
    na.rm = TRUE
  ) +
  geom_boxplot(na.rm = TRUE,
               fill = "transparent",
               outlier.shape = NA) +
  labs(y = "Depression Score (BDI-II)", x = "") +
  theme_bw() +
  scale_fill_manual(values = initial_depressed_color) +
  base_theme +
  geom_hline(yintercept = 14,
             colour = "black",
             linetype = "dotted") +
  scale_x_discrete(
    labels = c(
      "T1" = "Baseline",
      "T2" = "Mid-Retreat",
      "T3" = "End of Retreat",
      "T4" = "One Month",
      "T5" = "Three Months",
      "T6" = "Six Months",
      "T7" = "One Year"
    )
  ) +
  theme(
    legend.position = c(1, 1),
    legend.justification = c(1, 1),
    axis.text = element_text(size = 15),
    axis.title = element_text(size = 18)
  ) +
  facet_wrap(facets = vars(`Initial Depression`), nrow = 2)
plot

ggsave('Depression Boxplot2.pdf',
       height = 9,
       width = 10)

ggsave('Depression Boxplot2.png',
       height = 9,
       width = 10)

# # Next, make primals plots. Also color by initial depression.
#
# primals_plotting_func <- function(primal) {
#   theme_set(theme_bw())
#   ggplot(ibsr_data, aes_string(x = "Time", y = primal)) +
#     geom_boxplot(na.rm =
#                    TRUE, outlier.shape = NA) + geom_point(
#                      aes(colour = `Initial Depression`),
#                      size = 2.5,
#                      position = position_jitter(0.2, seed = 2),
#                      na.rm = TRUE
#                    ) + labs(y = paste(str_to_title(primal), "(PI-18)")) + theme(
#                      text = element_text(size = 18, family = "Arial"),
#                      legend.title.align = 0.5,
#                      panel.grid.major = element_blank(),
#                      panel.grid.minor = element_blank(),
#                      panel.background = element_rect(colour =
#                                                        "black"),
#                      axis.text.x = element_text(colour = "black"),
#                      axis.text.y = element_text(colour = "black"),
#                      axis.title.x = element_blank(),
#                      axis.title.y = element_text(face = "bold"),
#                      legend.position = "none"
#                    ) + scale_colour_discrete(type = c("#43CC47", "#1982FC"),
#                                              na.translate = F) + scale_x_discrete(
#                                                labels = c(
#                                                  "T1" = "Baseline",
#                                                  "T2" = "Mid-Retreat",
#                                                  "T3" = "End of Retreat",
#                                                  "T4" = "One Month",
#                                                  "T5" = "Three Months",
#                                                  "T6" = "Six Months",
#                                                  "T7" = "One Year"
#                                                )
#                                              )
#
#   ggsave(
#     paste(str_to_title(primal), 'Boxplot.png'),
#     height = 6,
#     width = 10,
#     units = "in"
#   )
#
#
#   sapply(c("alive", "safe", "good", "enticing"),
#          primals_plotting_func)
#
#
#   # Finally, make plots for everything else (no T2/T3 data). Also color by initial depression.
#
#
#   # Filter data frame
#   ibsr_data_no23 <-
#     filter(ibsr_data, !(Time == "T2" | Time == "T3"))
#
#   # Define measures of interest
#   measures <-
#     read.csv('Measures of interest.csv', header = FALSE)$V2
#   measures <- measures[6:length(measures)]
#
#   # Nice names
#   nice_measures <-
#     read.csv('Measures of interest.csv', header = FALSE)$V1
#   nice_measures <- nice_measures[6:length(nice_measures)]
#
#   no23_plotting_func <- function(measure, nice_measure) {
#     theme_set(theme_bw())
#     ggplot(ibsr_data_no23, aes_string(x = "Time", y = measure)) +
#       geom_boxplot(na.rm =
#                      TRUE, outlier.shape = NA) + geom_point(
#                        aes(colour = `Initial Depression`),
#                        size = 2.5,
#                        position = position_jitter(0.2, seed = 2),
#                        na.rm = TRUE
#                      ) + labs(y = nice_measure) + theme(
#                        text = element_text(size = 18, family = "Arial"),
#                        legend.title.align = 0.5,
#                        panel.grid.major = element_blank(),
#                        panel.grid.minor = element_blank(),
#                        panel.background = element_rect(colour =
#                                                          "black"),
#                        axis.text.x = element_text(colour = "black"),
#                        axis.text.y = element_text(colour = "black"),
#                        axis.title.x = element_blank(),
#                        axis.title.y = element_text(face = "bold"),
#                        legend.position = "none"
#                      ) + scale_colour_discrete(type = c("#43CC47", "#1982FC"),
#                                                na.translate = F) + scale_x_discrete(
#                                                  labels = c(
#                                                    "T1" = "Baseline",
#                                                    "T4" = "One Month",
#                                                    "T5" = "Three Months",
#                                                    "T6" = "Six Months",
#                                                    "T7" = "One Year"
#                                                  )
#                                                )
#
#     ggsave(
#       paste(nice_measure, 'Boxplot.png'),
#       height = 6,
#       width = 10,
#       units = "in"
#     )
#
#
#     mapply(no23_plotting_func, measures, nice_measures)
#
