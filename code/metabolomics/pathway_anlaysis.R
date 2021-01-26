###
no_source()

sxtTools::setwd_project()
setwd("data_analysis/metabolomics_data/Pathway Analysis/")

neg_dep_pathway =
  readr::read_csv("neg_dep_pathways.csv")
