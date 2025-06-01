###
no_source()

setwd(r4projects::get_project_wd())
rm(list = ls())

load("3_data_analysis/sample_info/sample_info")

total_sample_info <-
  sample_info

setwd(r4projects::get_project_wd())
setwd("3_data_analysis/transcriptomics/")

library(tidyverse)
library(data.table)

data = readr::read_delim("Ganz - Gene CPM Log2.txt")

colnames(data)

data <-
  data %>%
  dplyr::distinct(Gene, .keep_all = TRUE)

variable_info = data[, c(1)] %>%
  dplyr::rename(variable_id = Gene)

##remove duplicated gene
expression_data = data[-c(1, 140)] %>%
  as.data.frame()

sample_info =
  data.frame(sample_id = colnames(expression_data))

sample_info <-
  sample_info %>% 
  dplyr::mutate(sample_id = stringr::str_replace(sample_id, "-", "_"))

colnames(expression_data) <-
  sample_info$sample_id

match(sample_info$sample_id,
      total_sample_info$sample_id)

rownames(expression_data) <-
  variable_info$variable_id

###add more information
library(clusterProfiler)
library(org.Hs.eg.db)

variable_info <-
  variable_info %>%
  dplyr::mutate(SYMBOL = variable_id)

ENTREZID <-
  clusterProfiler::bitr(
    geneID = variable_info$SYMBOL,
    fromType = "SYMBOL",
    toType = "ENTREZID",
    OrgDb = org.Hs.eg.db
  ) %>%
  dplyr::distinct(SYMBOL, .keep_all = TRUE)

variable_info <-
  variable_info %>%
  dplyr::left_join(ENTREZID, by = "SYMBOL")

GENENAME <-
  clusterProfiler::bitr(
    geneID = variable_info$SYMBOL,
    fromType = "SYMBOL",
    toType = "GENENAME",
    OrgDb = org.Hs.eg.db
  ) %>%
  dplyr::distinct(SYMBOL, .keep_all = TRUE)

variable_info <-
  variable_info %>%
  dplyr::left_join(GENENAME, by = "SYMBOL")

ALIAS <-
  clusterProfiler::bitr(
    geneID = variable_info$SYMBOL,
    fromType = "SYMBOL",
    toType = "ALIAS",
    OrgDb = org.Hs.eg.db
  ) %>%
  dplyr::distinct(SYMBOL, .keep_all = TRUE)

variable_info <-
  variable_info %>%
  dplyr::left_join(ALIAS, by = "SYMBOL")

GENETYPE <-
  clusterProfiler::bitr(
    geneID = variable_info$SYMBOL,
    fromType = "SYMBOL",
    toType = "GENETYPE",
    OrgDb = org.Hs.eg.db
  ) %>%
  dplyr::distinct(SYMBOL, .keep_all = TRUE)

variable_info <-
  variable_info %>%
  dplyr::left_join(GENETYPE, by = "SYMBOL")

ACCNUM <-
  clusterProfiler::bitr(
    geneID = variable_info$SYMBOL,
    fromType = "SYMBOL",
    toType = "ACCNUM",
    OrgDb = org.Hs.eg.db
  ) %>%
  dplyr::distinct(SYMBOL, .keep_all = TRUE)

variable_info <-
  variable_info %>%
  dplyr::left_join(ACCNUM, by = "SYMBOL")

UCSCKG <-
  clusterProfiler::bitr(
    geneID = variable_info$SYMBOL,
    fromType = "SYMBOL",
    toType = "UCSCKG",
    OrgDb = org.Hs.eg.db
  ) %>%
  dplyr::distinct(SYMBOL, .keep_all = TRUE)

variable_info <-
  variable_info %>%
  dplyr::left_join(UCSCKG, by = "SYMBOL")

ENSEMBL <-
  clusterProfiler::bitr(
    geneID = variable_info$SYMBOL,
    fromType = "SYMBOL",
    toType = "ENSEMBL",
    OrgDb = org.Hs.eg.db
  ) %>%
  dplyr::distinct(SYMBOL, .keep_all = TRUE)

variable_info <-
  variable_info %>%
  dplyr::left_join(ENSEMBL, by = "SYMBOL")

UNIPROT <-
  clusterProfiler::bitr(
    geneID = variable_info$SYMBOL,
    fromType = "SYMBOL",
    toType = "UNIPROT",
    OrgDb = org.Hs.eg.db
  ) %>%
  dplyr::distinct(SYMBOL, .keep_all = TRUE)

variable_info <-
  variable_info %>%
  dplyr::left_join(UNIPROT, by = "SYMBOL")

dir.create("data_preparation")
setwd("data_preparation")

sample_info <-
  sample_info %>%
  dplyr::mutate(class = "Subject")

sample_info <-
sample_info %>% 
  dplyr::left_join(total_sample_info, by = "sample_id")

save(expression_data, file = "expression_data")
save(sample_info, file = "sample_info")
save(variable_info, file = "variable_info")

transcriptomics_data <-
  massdataset::create_mass_dataset(expression_data = expression_data, 
                                   sample_info = sample_info, 
                                   variable_info = variable_info)

save(transcriptomics_data, file = "transcriptomics_data")
