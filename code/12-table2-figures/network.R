no_function

setwd(r4projects::get_project_wd())

rm(list = ls())

data <-
  readxl::read_xlsx("data_analysis/some_figures/table2/data2.xlsx")

dir.create("data_analysis/some_figures/table2/network")
setwd("data_analysis/some_figures/table2/network")

###network to show the correlations between disease and molecules

edge_data <-
  data %>%
  dplyr::select(Disease, Metabolite, correlation, qvalue, Class) %>%
  dplyr::rename(node1 = Disease,
                node2 = Metabolite,
                node2_class = Class) %>%
  dplyr::mutate(node1_class = "disease")

node_data <-
  data.frame(
    node = c(edge_data$node1[1],
             edge_data$node2),
    class = c(edge_data$node1_class[1],
              edge_data$node2_class)
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
  dplyr::select(node, class, x, y)

coords =
  coords %>%
  dplyr::left_join(
    edge_data %>% dplyr::select(node2, correlation, qvalue, node2_class),
    by = c("node" = "node2")
  )


coords <-
  coords %>%
  dplyr::mutate(y1 = x) %>%
  dplyr::mutate(x1 = case_when(y == 1 ~ 0,
                               y == 0 ~ 1)) %>%
  dplyr::select(-c(x, y)) %>%
  dplyr::select(x = x1, y = y1, everything())

coords$x[coords$correlation < 0] = -1

my_graph <-
  create_layout(
    graph = g,
    layout = "manual",
    x = coords$x,
    y = coords$y
    # node.position = coords
  )

plot =
  ggraph(total_graph,
         layout = "manual",
         x = coords$x,
         y = coords$y) +
  geom_edge_diagonal(aes(color = correlation,
                         width = -log(qvalue, 10)),
                     alpha = 0.5,
                     show.legend = TRUE) +
  geom_node_point(
    aes(fill = class,
        size = Degree),
    shape = 21,
    alpha = 1,
    show.legend = TRUE
  ) +
  shadowtext::geom_shadowtext(
    aes(
      x = x,
      y = y,
      label = node,
      color = class
    ),
    bg.color = "white",
    size = 3,
    show.legend = FALSE
  ) +
  scale_size_continuous(range = c(5, 10)) +
  # scale_fill_manual(values = c(class_color, wearable_color)) +
  # scale_color_manual(values = c(class_color, wearable_color)) +
  guides(
    linetype = "none",
    color = guide_colorbar(title = "Correlation",
                           override.aes = list(linetype = "none")),
    size = guide_legend(
      title = "Degree",
      override.aes = list(
        linetype = NA,
        fill = "transparent",
        shape = 21,
        color = "black"
      )
    ),
    fill = guide_legend(
      title = "Class",
      override.aes = list(
        shape = 21,
        size = 3,
        alpha = 1
      )
    )
  ) +
  scale_edge_width_continuous(range = c(1, 4)) +
  ggraph::scale_edge_color_gradientn(colours = c(alpha("#3B4992FF", 0.7),
                                                 "white",
                                                 alpha("#EE0000FF", 0.7))) +
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
       filename = "network.pdf",
       width = 7,
       height = 7)

ggsave(plot,
       filename = "network.png",
       width = 7,
       height = 7)
