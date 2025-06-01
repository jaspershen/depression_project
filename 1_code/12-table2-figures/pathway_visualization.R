no_function

setwd(r4projects::get_project_wd())

rm(list = ls())

library(tidyverse)

data <-
  readxl::read_xlsx("3_data_analysis/some_figures/table2/data1.xlsx")

dir.create("3_data_analysis/some_figures/table2/pathway_visualization")
setwd("3_data_analysis/some_figures/table2/pathway_visualization")

data$qvalue[data$qvalue == "<0.001"] <- 0.001

data$qvalue <-
  as.numeric(data$qvalue)

###network to show the correlations between disease and molecules
node_data <-
  rbind(
    data %>%
      dplyr::select(node = Metabolite,
                    correlation,
                    qvalue) %>%
      dplyr::mutate(node_class = "metabolite"),
    data %>%
      dplyr::select(node = Class) %>%
      dplyr::mutate(
        correlation = NA,
        qvalue = 0.01,
        node_class = "metabolite_class"
      ),
    data %>%
      dplyr::select(node = Pathway) %>%
      dplyr::mutate(
        correlation = NA,
        qvalue = 0.01,
        node_class = "metabolite_pathway"
      )
  ) %>%
  dplyr::distinct(node, .keep_all = TRUE)

edge_data <-
  rbind(
    data %>%
      dplyr::select(from = Metabolite,
                    to = Class) %>%
      dplyr::mutate(
        from_class = "metabolite",
        to_class = "metabolite_class",
        edge_class = "class_metabolite"
      ),
    data %>%
      dplyr::select(from = Metabolite,
                    to = Pathway) %>%
      dplyr::mutate(
        from_class = "metabolite",
        to_class = "metabolite_pathway",
        edge_class = "pathway_metabolite"
      )
  )

library(tidygraph)
library(ggraph)
library(igraph)

total_graph <-
  tidygraph::tbl_graph(nodes = node_data,
                       edges = edge_data,
                       directed = FALSE) %>%
  dplyr::mutate(Degree = centrality_degree(mode = 'all'))

g <- total_graph

V(g)$type <- bipartite_mapping(g)$type

coords <-
  create_layout(g, layout = "bipartite") %>%
  dplyr::select(node, node_class, correlation, qvalue, x, y)

coords$index = 1:nrow(coords)

coords$y[coords$y == 0] <- 0.5

coords <-
  coords %>%
  dplyr::select(x, y) %>%
  dplyr::mutate(
    theta = x / (max(x) + 1) * 2 * pi,
    r = y + 1,
    x = r * cos(theta),
    y = r * sin(theta)
  )

my_graph <-
  create_layout(
    graph = g,
    layout = "manual",
    x = coords$x,
    y = coords$y
    # node.position = coords
  )

plot =
  ggraph(my_graph,
         layout = 'bipartite') +
  geom_edge_diagonal(
    aes(color = edge_class),
    width = 0.5,
    alpha = 0.5,
    show.legend = TRUE
  ) +
  geom_node_point(
    aes(fill = correlation,
        size = -log(qvalue, 10)),
    shape = 21,
    alpha = 1,
    show.legend = TRUE
  ) +
  geom_node_text(
    aes(
      x = x * 1.03,
      y = y * 1.03,
      label = node,
      hjust = ifelse(node_class == "metabolite", 'outward', "inward"),
      angle = -((-node_angle(x, y) + 90) %% 180) + 90
    ),
    repel = FALSE,
    size = 3,
    alpha = 1,
    show.legend = FALSE
  ) +
  scale_size_continuous(range = c(1, 5)) +
  scale_fill_gradientn(colours = c(alpha("#3B4992FF", 0.7),
                                   "white",
                                   alpha("#EE0000FF", 0.7)),
                       na.value = "black") +
  ggraph::theme_graph() +
  theme(
    plot.background = element_rect(fill = "transparent", color = NA),
    panel.background = element_rect(fill = "transparent", color = NA),
    legend.position = "right",
    legend.background = element_rect(fill = "transparent", color = NA)
  )

plot

extrafont::loadfonts()

ggsave(plot,
       filename = "pathway_network.pdf",
       width = 8.6,
       height = 7)

ggsave(plot,
       filename = "pathway_network.png",
       width = 8.6,
       height = 7)



###horizon plot
total_graph <-
  tidygraph::tbl_graph(nodes = node_data,
                       edges = edge_data,
                       directed = FALSE) %>%
  dplyr::mutate(Degree = centrality_degree(mode = 'all'))

g <- total_graph

V(g)$type <- bipartite_mapping(g)$type

coords <-
  create_layout(g, layout = "bipartite") %>%
  dplyr::select(node, node_class, correlation, qvalue, x, y)

coords$index = 1:nrow(coords)

coords$y[coords$y == 0] <- 0.5

# coords <-
#   coords %>%
#   dplyr::select(x, y) %>%
#   dplyr::mutate(
#     theta = x / (max(x) + 1) * 2 * pi,
#     r = y + 1,
#     x = r * cos(theta),
#     y = r * sin(theta)
#   )

my_graph <-
  create_layout(
    graph = g,
    layout = "manual",
    x = coords$x,
    y = coords$y
    # node.position = coords
  )

plot =
  ggraph(my_graph,
         layout = 'bipartite') +
  geom_edge_diagonal(
    aes(color = edge_class),
    width = 0.5,
    alpha = 0.5,
    show.legend = TRUE
  ) +
  geom_node_point(
    aes(fill = correlation,
        size = -log(qvalue, 10)),
    shape = 21,
    alpha = 1,
    show.legend = TRUE
  ) +
  geom_node_text(
    aes(
      x = x * 1.03,
      y = y * 1.03,
      label = node,
      hjust = ifelse(node_class == "metabolite", 'outward', "inward"),
      angle = 90
    ),
    repel = FALSE,
    size = 3,
    alpha = 1,
    show.legend = FALSE
  ) +
  scale_size_continuous(range = c(1, 5)) +
  scale_fill_gradientn(colours = c(alpha("#3B4992FF", 0.7),
                                   "white",
                                   alpha("#EE0000FF", 0.7)),
                       na.value = "black") +
  ggraph::theme_graph() +
  theme(
    plot.background = element_rect(fill = "transparent", color = NA),
    panel.background = element_rect(fill = "transparent", color = NA),
    legend.position = "right",
    legend.background = element_rect(fill = "transparent", color = NA)
  )

plot

extrafont::loadfonts()

ggsave(plot,
       filename = "pathway_network2.pdf",
       width = 8.6,
       height = 7)

ggsave(plot,
       filename = "pathway_network2.png",
       width = 8.6,
       height = 7)


# coords <-
#   create_layout(g, layout = "bipartite") %>%
#   dplyr::select(node, node_class, correlation, qvalue, x, y)
# 
# coords <-
#   coords %>%
#   dplyr::mutate(y1 = x) %>%
#   dplyr::mutate(x1 = case_when(y == 1 ~ 0,
#                                y == 0 ~ 1)) %>%
#   dplyr::select(-c(x, y)) %>%
#   dplyr::select(x = x1, y = y1, everything())
# 
# coords$x[coords$node_class == "metabolite"] = 0
# coords$x[coords$node_class == "metabolite_class"] = -2
# coords$x[coords$node_class == "metabolite_pathway"] = 2
# 
# # x = seq(from = 0,  to = 2 * pi,  by = 0.1)
# # y = sin(x)
# # plot(x, y, type = "l")
# #
# # coords$y[coords$x == 0]
# 
# coords <-
#   coords %>%
#   dplyr::mutate(y2 = y)
# 
# coords$y2 <-
#   ((coords$y2 - min(coords$y2)) / (max(coords$y2) -  min(coords$y2))) * (2 * pi - 0) + 0
# 
# coords$x[coords$x == 0] <-
#   sin(coords$y2[coords$x == 0])
# 
# my_graph <-
#   create_layout(
#     graph = g,
#     layout = "manual",
#     x = coords$x,
#     y = coords$y
#     # node.position = coords
#   )
# 
# plot =
#   ggraph(total_graph,
#          layout = "manual",
#          x = coords$x,
#          y = coords$y) +
#   geom_edge_diagonal(
#     aes(color = edge_class),
#     width = 0.5,
#     alpha = 0.5,
#     show.legend = TRUE
#   ) +
#   geom_node_point(
#     aes(fill = correlation,
#         size = -log(qvalue, 10)),
#     shape = 21,
#     alpha = 1,
#     show.legend = TRUE
#   ) +
#   geom_text(aes(
#     x = x,
#     y = y,
#     label = ifelse(node_class != "metabolite", node, NA)
#   ),
#   color = "black") +
#   scale_size_continuous(range = c(1, 3)) +
#   scale_fill_gradientn(colours = c(alpha("#3B4992FF", 0.7),
#                                    "white",
#                                    alpha("#EE0000FF", 0.7)),
#                        na.value = "black") +
#   ggraph::theme_graph() +
#   theme(
#     plot.background = element_rect(fill = "transparent", color = NA),
#     panel.background = element_rect(fill = "transparent", color = NA),
#     legend.position = "right",
#     legend.background = element_rect(fill = "transparent", color = NA)
#   )
# 
# plot
# 
# extrafont::loadfonts()
# 
# ggsave(plot,
#        filename = "pathway_network.pdf",
#        width = 7,
#        height = 7)
# 
# ggsave(plot,
#        filename = "pathway_network.png",
#        width = 7,
#        height = 7)
