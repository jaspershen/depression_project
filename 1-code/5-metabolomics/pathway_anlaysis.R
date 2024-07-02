###
no_source()

setwd(r4projects::get_project_wd())
rm(list = ls())
library(tidyverse)
library(data.table)

load("3-data_analysis/metabolomics_data/data_preparation/expression_data")
load("3-data_analysis/metabolomics_data/data_preparation/sample_info")
load("3-data_analysis/metabolomics_data/data_preparation/variable_info")

setwd("3-data_analysis/metabolomics_data/Pathway Analysis/")

combined_dep_pathways =
  readr::read_csv("combined_dep_pathways.csv")

temp_data =
  head(combined_dep_pathways, 10)

temp_data = 
purrr::map2(
  .x = temp_data$Cpd.Hits,
  .y = temp_data$X1,
  .f = function(x, y) {
    x = stringr::str_split(x, ";")[[1]]
    idx = lapply(x, function(y) {
      grep(y, variable_info$KEGG_ID)[1]
    }) %>% 
      unlist()
    idx = idx[!is.na(idx)]
    if(length(idx) == 0){
      return(NULL)
    }
    data.frame(
      pathway = y,
      variable_id = variable_info$Name[idx],
      expression_data[idx, ],
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
  }
) %>% 
  do.call(rbind, .) %>% 
  as.data.frame()

###important metabolite heatmap
library(ComplexHeatmap)

plot =
  proteome_gene_sample_level_heatmap(
    result_go_cluster = result_go_cluster,
    result_kegg_cluster = result_kegg_cluster,
    result_reactome_cluster = result_reactome_cluster,
    subject_data = subject_data,
    protein_marker = protein_marker,
    sample_info = sample_info,
    path_level = path_level
  )

plot

# ggsave(plot,
#        filename = "top10_pathway_gene_level_heatmap.pdf",
#        width = 16, height = 8)


####graph and heatmap
    library(plyr)

    edge_data <- temp_data[,c(1:2)]
    node_data <- temp_data[,c(1:2)]
    
    node_data1 <-
      node_data %>%
      dplyr::select(pathway) %>%
      dplyr::rename(name = pathway) %>%
      dplyr::distinct(name, .keep_all = TRUE) %>%
      dplyr::mutate(class = "pathway") %>%
      dplyr::mutate(path_name = name)
  
    node_data2 <-
      node_data %>%
      dplyr::select(variable_id) %>%
      dplyr::rename(
        name = variable_id,
      ) %>%
      dplyr::distinct(name, .keep_all = TRUE) %>%
      dplyr::mutate(class = "metabolite") %>%
      dplyr::mutate(path_name = NA)

    node_data <- rbind(node_data1, node_data2)
    node_data$path_name[is.na(node_data$path_name)] = "metabolite"
    
    edge_data <-
      edge_data %>%
      dplyr::select(from = pathway,
                    to = variable_id) %>%
      dplyr::mutate(path = from)
    
    library(tidygraph)
    library(ggraph)
    
    total_graph <-
      tidygraph::tbl_graph(nodes = node_data,
                           edges = edge_data,
                           directed = FALSE) %>%
      dplyr::mutate(Degree = centrality_degree(mode = 'all'))
    
    library(igraph)
    
    g <- total_graph
    
    V(g)$type <- bipartite_mapping(g)$type
    
    coords <-
      create_layout(g, layout = "bipartite") %>%
      dplyr::select(name, class, x, y)
    
    coords$index = 1:nrow(coords)
    
    coords =
      coords %>%
      plyr::dlply(.variables = .(class)) %>%
      purrr::map(
        .f = function(x) {
          x =
            x %>%
            dplyr::arrange(x)
          x$x =
            seq(from = 1,
                to = 100,
                length.out = nrow(x))
          
          x %>%
            dplyr::arrange(index)
          
        }
      ) %>%
      do.call(rbind, .) %>%
      as.data.frame() %>%
      dplyr::arrange(index)
    
    coords =
      coords %>%
      dplyr::mutate(x1 = y, y1 = x) %>%
      dplyr::select(-c(x, y)) %>%
      dplyr::rename(x = x1, y = y1) 
    
    coords$x[coords$x == 0] <- 2
    
    my_graph <-
      create_layout(
        graph = g,
        layout = "manual",
        x = coords$x,
        y = coords$y
        # node.position = coords
      )
    
    library(ggraph)
    
    # RColorBrewer::display.brewer.all()
    path_col <-
      colorRampPalette(colors =
                         RColorBrewer::brewer.pal(n = 11,
                                                  name = "Spectral"))(length(unique(edge_data$path)))
    
    names(path_col) <- unique(edge_data$path)
    path_col = c(path_col, "metabolite" = "grey")
    plot1 <-
      ggraph(my_graph,
             layout = 'bipartite') +
      geom_edge_diagonal(
        strength = 1,
        aes(color = path),
        edge_width = 0.5,
        alpha = 1,
        show.legend = FALSE
      ) +
      geom_node_point(
        aes(fill = path_name,
            size = Degree),
        shape = 21,
        alpha = 1,
        show.legend = FALSE
      ) +
      geom_node_text(aes(
        x = x,
        y = y,
        hjust = 1,
        size = 3,
        label = ifelse(class == "pathway", name, NA)
      ),
      show.legend = FALSE) +
      ggraph::scale_edge_color_manual(values = path_col) +
      scale_size_continuous(range = c(4, 8)) +
      scale_fill_manual(values = path_col) +
      ggraph::theme_graph() +
      theme(
        plot.background = element_rect(fill = "transparent", color = NA),
        panel.background = element_rect(fill = "transparent", color = NA),
        legend.position = "right",
        legend.background = element_rect(fill = "transparent", color = NA)
      )
    
    plot1
    
    ###important genes heatmap
    library(ComplexHeatmap)
    
    col_fun = colorRamp2(c(-3, 0, 3), c("#4292C6", "white", "red"))
    ###quantitative pathways
    
    library(plyr)
    rownames(temp_data) = NULL
    
    temp_data2 = 
      temp_data %>% 
      dplyr::distinct(variable_id, .keep_all = TRUE) %>% 
      tibble::column_to_rownames(var = "variable_id") %>% 
      dplyr::select(one_of(sample_info$sample_id)) %>% 
      t() %>%
      as.data.frame() %>%
      data.frame(time = sample_info$Time, ., check.names = FALSE) %>% 
    plyr::dlply(.variables = .(time)) %>% 
      purrr::map(function(x){
        temp = x[,-1] %>%
          apply(1, function(y){
            y = as.numeric(y)
            (y )/sd(y)
          }) %>% 
          t() %>% 
          colMeans()
      names(temp) = colnames(x)[-1]    
      temp
      }) %>% 
      do.call(rbind, .) %>% 
      t() %>% 
      as.data.frame()
    
    temp_data2 = 
    temp_data2[coords %>%
                 dplyr::filter(class == "metabolite") %>%
                 dplyr::arrange(desc(y)) %>%
                 pull(name), ] %>% 
      apply(1, function(x){
        (x - mean(x))/sd(x)
      }) %>% 
      t() %>% 
      as.data.frame()
    
    if (abs(range(temp_data2)[1]) > range(temp_data2)[2]) {
      temp_data2[temp_data2 < -range(temp_data2)[2]] <-
        -range(temp_data2)[2]
    } else{
      temp_data2[temp_data2 > -range(temp_data2)[1]] <-
        -range(temp_data2)[1]
    }
    
    library(circlize)
    
    plot2 <-
      Heatmap(
        as.matrix(temp_data2),
        cluster_columns = FALSE,
        cluster_rows = FALSE,
        show_row_names = TRUE,
        show_column_names = TRUE,
        border = TRUE,
        col = col_fun,
        name = "z-score",
        clustering_method_rows = "ward.D",
        rect_gp = gpar(col = "white")
      )
    
    library(ggplotify)
    
    plot2 <- as.ggplot(plot2)
    
    library(patchwork)
    
    plot =
      plot1 + plot2 + patchwork::plot_layout(widths = c(1, 4))
    
  plot


ggsave(plot,
       filename = "top10_pathway_gene_level_ga_heatmap.pdf",
       width = 10,
       height = 7)
