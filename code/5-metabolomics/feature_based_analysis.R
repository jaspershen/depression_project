###
no_source()

setwd(masstools::get_project_wd())
rm(list = ls())
library(tidyverse)
library(data.table)

###load data
load("data_analysis/metabolomics_data/data_preparation/peaks/variable_info")

##load marker
load("data_analysis/metabolomics_data/marker/depression_association_pos")
load("data_analysis/metabolomics_data/marker/depression_association_neg")

dim(depression_association_pos)
dim(depression_association_neg)

setwd("data_analysis/metabolomics_data/PIUMet")

piumet_file = rbind(
  depression_association_pos %>%
    data.frame(polarity = "positive", .),
  depression_association_neg %>%
    data.frame(polarity = "negative", .)
) %>%
  dplyr::select(name = Variables,
                polarity, correlation) %>%
  dplyr::mutate(correlation = abs(correlation))

piumet_file = 
piumet_file %>%
  dplyr::left_join(variable_info[, c("variable_id", "mz")], 
                   by = c("name" = "variable_id")) %>% 
  dplyr::select(mz, polarity, correlation)


colnames(piumet_file) <- NULL

write.table(
  piumet_file,
  "piumet_file.txt",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)








load("PIUMet/piumet_output_up/Result/node_data")

hmdb_id <-
  node_data %>%
  dplyr::filter(!is.null(HMDB_ID)) %>%
  dplyr::filter(HMDB_ID != " " & HMDB_ID != "NULL") %>%
  dplyr::pull(HMDB_ID) %>%
  unique()

kegg_id <-
  lapply(hmdb_id, function(x){
    metflow2::transID(query = x,
                      from = "Human Metabolome Database",
                      to = "KEGG", top = 1)
  }) %>%
  do.call(rbind, .)

kegg_id

load("PIUMet/hmdbMS1Database0.0.1")

hmdb_data <- hmdbMS1Database0.0.1@spectra.info

kegg_id2 <-
  hmdb_data$KEGG.ID[match(kegg_id$`Human Metabolome Database`, hmdb_data$HMDB.ID)]

kegg_id <- data.frame(kegg_id, kegg_id2, stringsAsFactors = FALSE)

colnames(kegg_id) <- c("HMDB_ID", "KEGG_ID1", "KEGG_ID2")

KEGG_ID <-
  apply(kegg_id, 1, function(x){
    x <- as.character(x)
    if(is.na(x[2]) & is.na(x[3])){
      return(NA)
    }

    if(!is.na(x[2])){
      return(x[2])
    }

    if(!is.na(x[3])){
      return(x[3])
    }


  })

kegg_id <-
  kegg_id %>%
  dplyr::mutate(KEGG_ID = KEGG_ID) %>%
  dplyr::select(-c(KEGG_ID1, KEGG_ID2))

kegg_id

node_data <-
  node_data %>%
  dplyr::left_join(kegg_id, by = c("HMDB_ID"))

annotation_result <- node_data

load("PIUMet/piumet_output_up/Result/edge_data")

edge_data <-
edge_data %>%
  dplyr::filter(stringr::str_detect(from, "POS") |
                  stringr::str_detect(from, "NEG") |
                  stringr::str_detect(to, "POS") |
                  stringr::str_detect(to, "NEG")
                  ) %>%
  dplyr::select(from, to) %>%
  t() %>%
  as.data.frame() %>%
  purrr::map(.f = function(x){
    if(stringr::str_detect(x[1], "POS") | stringr::str_detect(x[1], "NEG")){
      return(x)
    }else{
      return(rev(x))
    }
  }) %>%
  do.call(rbind, .) %>%
  as.data.frame() %>%
  dplyr::arrange(V1) %>%
  dplyr::distinct()

colnames(edge_data) <- c('peak', "metabolite")

annotation_result <-
annotation_result %>%
  dplyr::filter(node_class == "Metabolite") %>%
  dplyr::select(node, HMDB_ID, KEGG_ID, super.class)

annotation_result_up <-
edge_data %>%
  dplyr::left_join(annotation_result, by = c("metabolite" = "node"))

save(annotation_result_up,
     file = "PIUMet/piumet_output_up/Result/annotation_result_up")
