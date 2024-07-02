###
###test sync
###test snyc in windows
###test again
no_source()


###RPLC 
###positive mode
# setwd(r4projects::get_project_wd())
# setwd("3-data_analysis/metabolomics_data/MSMS Annotation/RPLC pos/")
# library(tidyverse)
# library(data.table)
# library(metID)
# 
# # peak_table = readr::read_csv("Ariel_Study_RPLC_Pos.csv")
# # 
# # peak_table = 
# # peak_table %>% 
# #   dplyr::select(name = Compound, mz = `m/z`, rt = `Retention time (min)`) %>% 
# #   dplyr::mutate(rt = rt * 60)
# # 
# # write.csv(peak_table, "peak_table.csv", row.names = FALSE)
# # 
# # ###level 1
# # param1 <-
# #   identify_metabolites_params(
# #     ms1.match.ppm = 25,
# #     rt.match.tol = 30,
# #     polarity = "positive",
# #     ce = "all",
# #     column = "rp",
# #     total.score.tol = 0.5,
# #     candidate.num = 3,
# #     threads = 3,
# #     database = "msDatabase_rplc0.0.2"
# #   )
# # 
# # ###level 2
# # param2 <-
# #   identify_metabolites_params(
# #     ms1.match.ppm = 25,
# #     rt.match.tol = 1000000,
# #     polarity = "positive",
# #     ce = "all",
# #     column = "rp",
# #     total.score.tol = 0.5,
# #     candidate.num = 3,
# #     threads = 3,
# #     database = "hmdbDatabase0.0.2"
# #   )
# # 
# # param3 <-
# #   identify_metabolites_params(
# #     ms1.match.ppm = 25,
# #     rt.match.tol = 1000000,
# #     polarity = "positive",
# #     ce = "all",
# #     column = "rp",
# #     total.score.tol = 0.5,
# #     candidate.num = 3,
# #     threads = 3,
# #     database = "massbankDatabase0.0.2"
# #   )
# # 
# # param4 <-
# #   identify_metabolites_params(
# #     ms1.match.ppm = 25,
# #     rt.match.tol = 1000000,
# #     polarity = "positive",
# #     ce = "all",
# #     column = "rp",
# #     total.score.tol = 0.5,
# #     candidate.num = 3,
# #     threads = 3,
# #     database = "metlinDatabase0.0.2"
# #   )
# # 
# # param5 <-
# #   identify_metabolites_params(
# #     ms1.match.ppm = 25,
# #     rt.match.tol = 1000000,
# #     polarity = "positive",
# #     ce = "all",
# #     column = "rp",
# #     total.score.tol = 0.5,
# #     candidate.num = 3,
# #     threads = 3,
# #     database = "monaDatabase0.0.2"
# #   )
# # 
# # param6 <-
# #   identify_metabolites_params(
# #     ms1.match.ppm = 25,
# #     rt.match.tol = 1000000,
# #     polarity = "positive",
# #     ce = "all",
# #     column = "rp",
# #     total.score.tol = 0.5,
# #     candidate.num = 3,
# #     threads = 3,
# #     database = "nistDatabase0.0.2"
# #   )
# # 
# # param7 <-
# #   identify_metabolites_params(
# #     ms1.match.ppm = 25,
# #     rt.match.tol = 1000000,
# #     polarity = "positive",
# #     ce = "all",
# #     column = "rp",
# #     total.score.tol = 0.5,
# #     candidate.num = 3,
# #     threads = 3,
# #     database = "orbitrapDatabase0.0.1"
# #   )
# # 
# # ##level 3
# # param8 <-
# #   identify_metabolites_params(
# #     ms1.match.ppm = 25,
# #     rt.match.tol = 1000000,
# #     polarity = "positive",
# #     ce = "all",
# #     column = "rp",
# #     total.score.tol = 0.5,
# #     candidate.num = 3,
# #     threads = 3,
# #     database = "hmdbMS1Database0.0.1"
# #   )
# # 
# # 
# # 
# # result_rp_pos25 <- identify_metabolite_all(
# #   ms1.data = "peak_table.csv",
# #   ms2.data = c("QC1_MSMS_NCE25.mzXML"),
# #   parameter.list = c(param1, param2, param3, param4, param5, param6, param7, param8),
# #   path = "."
# # )
# # 
# # save(result_rp_pos25, file = "result_rp_pos25")
# # 
# # result_rp_pos50 <- identify_metabolite_all(
# #   ms1.data = "peak_table.csv",
# #   ms2.data = c("QC2_MSMS_NCE50.mzXML"),
# #   parameter.list = c(param1, param2, param3, param4, param5, param6, param7),
# #   path = "."
# # )
# # 
# # save(result_rp_pos50, file = "result_rp_pos50")
# 
# load("result_rp_pos25")
# load("result_rp_pos50")
# 
# annotation_table1 <-
#   get_identification_table(result_rp_pos25[[1]], 
#                            result_rp_pos50[[1]], 
#                            type = "new", 
#                            candidate.num = 1)
# 
# annotation_table1 = 
# annotation_table1 %>% 
#   dplyr::filter(!is.na(Compound.name))
# 
# 
# annotation_table2 <-
#   get_identification_table(result_rp_pos25[[2]], 
#                            result_rp_pos25[[3]], 
#                            result_rp_pos25[[4]], 
#                            result_rp_pos25[[5]], 
#                            result_rp_pos25[[6]], 
#                            result_rp_pos25[[7]], 
#                            result_rp_pos50[[2]], 
#                            result_rp_pos50[[3]], 
#                            result_rp_pos50[[4]], 
#                            result_rp_pos50[[5]], 
#                            result_rp_pos50[[6]], 
#                            result_rp_pos50[[7]], 
#                            type = "new", 
#                            candidate.num = 1)
# 
# annotation_table2 = 
#   annotation_table2 %>% 
#   dplyr::filter(!is.na(Compound.name)) %>% 
#   dplyr::filter(!name %in% annotation_table1$name)
# 
# 
# annotation_table3 <-
#   get_identification_table(result_rp_pos25[[8]], 
#                            type = "new", 
#                            candidate.num = 1)
# 
# annotation_table3 = 
#   annotation_table3 %>% 
#   dplyr::filter(!is.na(Compound.name)) %>% 
#   dplyr::filter(!name %in% annotation_table1$name) %>% 
#   dplyr::filter(!name %in% annotation_table2$name)
# 
# 
# annotation_table1 <-
#   data.frame(annotation_table1,
#              Level = 1,
#              stringsAsFactors = FALSE)
# 
# annotation_table2 <-
#   data.frame(annotation_table2,
#              Level = 2,
#              stringsAsFactors = FALSE)
# 
# annotation_table3 <-
#   data.frame(annotation_table3,
#              Level = 3,
#              stringsAsFactors = FALSE)
# 
# annotation_table = 
# rbind(annotation_table1,
#       annotation_table2) %>%
#   dplyr::full_join(annotation_table3, by = intersect(colnames(annotation_table1),
#                                                      colnames(annotation_table3)))
# 
# write.csv(annotation_table, file = "annotation_table.csv", row.names = FALSE)
# 
# 
# 
# 
# ###RPLC 
# ###negative mode
# setwd(r4projects::get_project_wd())
# setwd("3-data_analysis/metabolomics_data/MSMS Annotation/RPLC neg/")
# library(tidyverse)
# library(data.table)
# library(metID)
# 
# # peak_table = readr::read_csv("Ariel_Study_RPLC_neg.csv")
# # 
# # peak_table = 
# #   peak_table %>% 
# #   dplyr::select(name = Compound, mz = `m/z`, rt = `Retention time (min)`) %>% 
# #   dplyr::mutate(rt = rt * 60)
# # 
# # write.csv(peak_table, "peak_table.csv", row.names = FALSE)
# # 
# # ###level 1
# # param1 <-
# #   identify_metabolites_params(
# #     ms1.match.ppm = 25,
# #     rt.match.tol = 30,
# #     polarity = "negative",
# #     ce = "all",
# #     column = "rp",
# #     total.score.tol = 0.5,
# #     candidate.num = 3,
# #     threads = 3,
# #     database = "msDatabase_rplc0.0.2"
# #   )
# # 
# # ###level 2
# # param2 <-
# #   identify_metabolites_params(
# #     ms1.match.ppm = 25,
# #     rt.match.tol = 1000000,
# #     polarity = "negative",
# #     ce = "all",
# #     column = "rp",
# #     total.score.tol = 0.5,
# #     candidate.num = 3,
# #     threads = 3,
# #     database = "hmdbDatabase0.0.2"
# #   )
# # 
# # param3 <-
# #   identify_metabolites_params(
# #     ms1.match.ppm = 25,
# #     rt.match.tol = 1000000,
# #     polarity = "negative",
# #     ce = "all",
# #     column = "rp",
# #     total.score.tol = 0.5,
# #     candidate.num = 3,
# #     threads = 3,
# #     database = "massbankDatabase0.0.2"
# #   )
# # 
# # param4 <-
# #   identify_metabolites_params(
# #     ms1.match.ppm = 25,
# #     rt.match.tol = 1000000,
# #     polarity = "negative",
# #     ce = "all",
# #     column = "rp",
# #     total.score.tol = 0.5,
# #     candidate.num = 3,
# #     threads = 3,
# #     database = "metlinDatabase0.0.2"
# #   )
# # 
# # param5 <-
# #   identify_metabolites_params(
# #     ms1.match.ppm = 25,
# #     rt.match.tol = 1000000,
# #     polarity = "negative",
# #     ce = "all",
# #     column = "rp",
# #     total.score.tol = 0.5,
# #     candidate.num = 3,
# #     threads = 3,
# #     database = "monaDatabase0.0.2"
# #   )
# # 
# # param6 <-
# #   identify_metabolites_params(
# #     ms1.match.ppm = 25,
# #     rt.match.tol = 1000000,
# #     polarity = "negative",
# #     ce = "all",
# #     column = "rp",
# #     total.score.tol = 0.5,
# #     candidate.num = 3,
# #     threads = 3,
# #     database = "nistDatabase0.0.2"
# #   )
# # 
# # param7 <-
# #   identify_metabolites_params(
# #     ms1.match.ppm = 25,
# #     rt.match.tol = 1000000,
# #     polarity = "negative",
# #     ce = "all",
# #     column = "rp",
# #     total.score.tol = 0.5,
# #     candidate.num = 3,
# #     threads = 3,
# #     database = "orbitrapDatabase0.0.1"
# #   )
# # 
# # ##level 3
# # param8 <-
# #   identify_metabolites_params(
# #     ms1.match.ppm = 25,
# #     rt.match.tol = 1000000,
# #     polarity = "negative",
# #     ce = "all",
# #     column = "rp",
# #     total.score.tol = 0.5,
# #     candidate.num = 3,
# #     threads = 3,
# #     database = "hmdbMS1Database0.0.1"
# #   )
# # 
# # 
# # 
# # result_rp_neg25 <- identify_metabolite_all(
# #   ms1.data = "peak_table.csv",
# #   ms2.data = c("QC1_MSMS_NCE25.mzXML"),
# #   parameter.list = c(param1, param2, param3, param4, param5, param6, param7, param8),
# #   path = "."
# # )
# # 
# # save(result_rp_neg25, file = "result_rp_neg25")
# 
# 
# result_rp_neg50 <- identify_metabolite_all(
#   ms1.data = "peak_table.csv",
#   ms2.data = c("QC2_MSMS_NCE50.mzXML"),
#   parameter.list = c(param1, param2, param3, param4, param5, param6, param7),
#   path = "."
# )
# 
# save(result_rp_neg50, file = "result_rp_neg50")
# 
# 
# annotation_table1 <-
#   get_identification_table(result_rp_neg25[[1]], 
#                            result_rp_neg50[[1]], 
#                            type = "new", 
#                            candidate.num = 1)
# 
# annotation_table1 = 
#   annotation_table1 %>% 
#   dplyr::filter(!is.na(Compound.name))
# 
# 
# annotation_table2 <-
#   get_identification_table(result_rp_neg25[[2]], 
#                            result_rp_neg25[[3]], 
#                            result_rp_neg25[[4]], 
#                            result_rp_neg25[[5]], 
#                            result_rp_neg25[[6]], 
#                            result_rp_neg25[[7]], 
#                            result_rp_neg50[[2]], 
#                            result_rp_neg50[[3]], 
#                            result_rp_neg50[[4]], 
#                            result_rp_neg50[[5]], 
#                            result_rp_neg50[[6]], 
#                            result_rp_neg50[[7]], 
#                            type = "new", 
#                            candidate.num = 1)
# 
# annotation_table2 = 
#   annotation_table2 %>% 
#   dplyr::filter(!is.na(Compound.name)) %>% 
#   dplyr::filter(!name %in% annotation_table1$name)
# 
# 
# annotation_table3 <-
#   get_identification_table(result_rp_neg25[[8]], 
#                            type = "new", 
#                            candidate.num = 1)
# 
# annotation_table3 = 
#   annotation_table3 %>% 
#   dplyr::filter(!is.na(Compound.name)) %>% 
#   dplyr::filter(!name %in% annotation_table1$name) %>% 
#   dplyr::filter(!name %in% annotation_table2$name)
# 
# annotation_table1 <-
#   data.frame(annotation_table1,
#              Level = 1,
#              stringsAsFactors = FALSE)
# 
# annotation_table2 <-
#   data.frame(annotation_table2,
#              Level = 2,
#              stringsAsFactors = FALSE)
# 
# annotation_table3 <-
#   data.frame(annotation_table3,
#              Level = 3,
#              stringsAsFactors = FALSE)
# 
# annotation_table = 
#   rbind(annotation_table1,
#         annotation_table2) %>%
#   dplyr::full_join(annotation_table3, by = intersect(colnames(annotation_table1),
#                                                      colnames(annotation_table3)))
# 
# write.csv(annotation_table, file = "annotation_table.csv", row.names = FALSE)
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# ###HILIC 
# ###positive mode
# setwd(r4projects::get_project_wd())
# setwd("3-data_analysis/metabolomics_data/MSMS Annotation/HILIC pos/")
# library(tidyverse)
# library(data.table)
# library(metID)
# 
# peak_table = readr::read_csv("Ariel_Study_HILIC_Pos.csv")
# 
# peak_table =
# peak_table %>%
#   dplyr::select(name = Compound, mz = `m/z`, rt = `Retention time (min)`) %>%
#   dplyr::mutate(rt = rt * 60)
# 
# write.csv(peak_table, "peak_table.csv", row.names = FALSE)
# 
# ###level 1
# param1 <-
#   identify_metabolites_params(
#     ms1.match.ppm = 25,
#     rt.match.tol = 30,
#     polarity = "positive",
#     ce = "all",
#     column = "hilic",
#     total.score.tol = 0.5,
#     candidate.num = 3,
#     threads = 8,
#     database = "msDatabase_hilic0.0.2"
#   )
# 
# ###level 2
# param2 <-
#   identify_metabolites_params(
#     ms1.match.ppm = 25,
#     rt.match.tol = 1000000,
#     polarity = "positive",
#     ce = "all",
#     column = "hilic",
#     total.score.tol = 0.5,
#     candidate.num = 3,
#     threads = 8,
#     database = "hmdbDatabase0.0.2"
#   )
# 
# param3 <-
#   identify_metabolites_params(
#     ms1.match.ppm = 25,
#     rt.match.tol = 1000000,
#     polarity = "positive",
#     ce = "all",
#     column = "hilic",
#     total.score.tol = 0.5,
#     candidate.num = 3,
#     threads = 8,
#     database = "massbankDatabase0.0.2"
#   )
# 
# param4 <-
#   identify_metabolites_params(
#     ms1.match.ppm = 25,
#     rt.match.tol = 1000000,
#     polarity = "positive",
#     ce = "all",
#     column = "hilic",
#     total.score.tol = 0.5,
#     candidate.num = 3,
#     threads = 8,
#     database = "metlinDatabase0.0.2"
#   )
# 
# param5 <-
#   identify_metabolites_params(
#     ms1.match.ppm = 25,
#     rt.match.tol = 1000000,
#     polarity = "positive",
#     ce = "all",
#     column = "hilic",
#     total.score.tol = 0.5,
#     candidate.num = 3,
#     threads = 8,
#     database = "monaDatabase0.0.2"
#   )
# 
# param6 <-
#   identify_metabolites_params(
#     ms1.match.ppm = 25,
#     rt.match.tol = 1000000,
#     polarity = "positive",
#     ce = "all",
#     column = "hilic",
#     total.score.tol = 0.5,
#     candidate.num = 3,
#     threads = 8,
#     database = "nistDatabase0.0.2"
#   )
# 
# param7 <-
#   identify_metabolites_params(
#     ms1.match.ppm = 25,
#     rt.match.tol = 1000000,
#     polarity = "positive",
#     ce = "all",
#     column = "hilic",
#     total.score.tol = 0.5,
#     candidate.num = 3,
#     threads = 8,
#     database = "orbitrapDatabase0.0.1"
#   )
# 
# ##level 3
# param8 <-
#   identify_metabolites_params(
#     ms1.match.ppm = 25,
#     rt.match.tol = 1000000,
#     polarity = "positive",
#     ce = "all",
#     column = "hilic",
#     total.score.tol = 0.5,
#     candidate.num = 3,
#     threads = 8,
#     database = "hmdbMS1Database0.0.1"
#   )
# 
# 
# 
# result_hilic_pos25 <- identify_metabolite_all(
#   ms1.data = "peak_table.csv",
#   ms2.data = c("QC1_MSMS_NCE25.mzXML"),
#   parameter.list = c(param1, param2, param3, param4, param5, param6, param7, param8),
#   path = "."
# )
# 
# save(result_hilic_pos25, file = "result_hilic_pos25")
# 
# result_hilic_pos50 <- identify_metabolite_all(
#   ms1.data = "peak_table.csv",
#   ms2.data = c("QC2_MSMS_NCE35.mzXML"),
#   parameter.list = c(param1, param2, param3, param4, param5, param6, param7),
#   path = "."
# )
# 
# save(result_hilic_pos50, file = "result_hilic_pos50")
# 
# load("result_hilic_pos25")
# load("result_hilic_pos50")
# 
# annotation_table1 <-
#   get_identification_table(result_hilic_pos25[[1]], 
#                            result_hilic_pos50[[1]], 
#                            type = "new", 
#                            candidate.num = 1)
# 
# annotation_table1 = 
#   annotation_table1 %>% 
#   dplyr::filter(!is.na(Compound.name))
# 
# 
# annotation_table2 <-
#   get_identification_table(result_hilic_pos25[[2]], 
#                            result_hilic_pos25[[3]], 
#                            result_hilic_pos25[[4]], 
#                            result_hilic_pos25[[5]], 
#                            result_hilic_pos25[[6]], 
#                            result_hilic_pos25[[7]], 
#                            result_hilic_pos50[[2]], 
#                            result_hilic_pos50[[3]], 
#                            result_hilic_pos50[[4]], 
#                            result_hilic_pos50[[5]], 
#                            result_hilic_pos50[[6]], 
#                            result_hilic_pos50[[7]], 
#                            type = "new", 
#                            candidate.num = 1)
# 
# annotation_table2 = 
#   annotation_table2 %>% 
#   dplyr::filter(!is.na(Compound.name)) %>% 
#   dplyr::filter(!name %in% annotation_table1$name)
# 
# 
# annotation_table3 <-
#   get_identification_table(result_hilic_pos25[[8]], 
#                            type = "new", 
#                            candidate.num = 1)
# 
# annotation_table3 = 
#   annotation_table3 %>% 
#   dplyr::filter(!is.na(Compound.name)) %>% 
#   dplyr::filter(!name %in% annotation_table1$name) %>% 
#   dplyr::filter(!name %in% annotation_table2$name)
# 
# 
# annotation_table1 <-
#   data.frame(annotation_table1,
#              Level = 1,
#              stringsAsFactors = FALSE)
# 
# annotation_table2 <-
#   data.frame(annotation_table2,
#              Level = 2,
#              stringsAsFactors = FALSE)
# 
# annotation_table3 <-
#   data.frame(annotation_table3,
#              Level = 3,
#              stringsAsFactors = FALSE)
# 
# annotation_table = 
#   rbind(annotation_table1,
#         annotation_table2) %>%
#   dplyr::full_join(annotation_table3, by = intersect(colnames(annotation_table1),
#                                                      colnames(annotation_table3)))
# 
# write.csv(annotation_table, file = "annotation_table.csv", row.names = FALSE)
# 
# 
# 
# 
# ###HILIC 
# ###negative mode
# setwd(r4projects::get_project_wd())
# setwd("3-data_analysis/metabolomics_data/MSMS Annotation/HILIC neg/")
# library(tidyverse)
# library(data.table)
# library(metID)
# 
# peak_table = readr::read_csv("Ariel_Study_HILIC_neg.csv")
# 
# peak_table =
#   peak_table %>%
#   dplyr::select(name = Compound, mz = `m/z`, rt = `Retention time (min)`) %>%
#   dplyr::mutate(rt = rt * 60)
# 
# write.csv(peak_table, "peak_table.csv", row.names = FALSE)
# 
# ###level 1
# param1 <-
#   identify_metabolites_params(
#     ms1.match.ppm = 25,
#     rt.match.tol = 30,
#     polarity = "negative",
#     ce = "all",
#     column = "hilic",
#     total.score.tol = 0.5,
#     candidate.num = 3,
#     threads = 8,
#     database = "msDatabase_hilic0.0.2"
#   )
# 
# ###level 2
# param2 <-
#   identify_metabolites_params(
#     ms1.match.ppm = 25,
#     rt.match.tol = 1000000,
#     polarity = "negative",
#     ce = "all",
#     column = "hilic",
#     total.score.tol = 0.5,
#     candidate.num = 3,
#     threads = 8,
#     database = "hmdbDatabase0.0.2"
#   )
# 
# param3 <-
#   identify_metabolites_params(
#     ms1.match.ppm = 25,
#     rt.match.tol = 1000000,
#     polarity = "negative",
#     ce = "all",
#     column = "hilic",
#     total.score.tol = 0.5,
#     candidate.num = 3,
#     threads = 8,
#     database = "massbankDatabase0.0.2"
#   )
# 
# param4 <-
#   identify_metabolites_params(
#     ms1.match.ppm = 25,
#     rt.match.tol = 1000000,
#     polarity = "negative",
#     ce = "all",
#     column = "hilic",
#     total.score.tol = 0.5,
#     candidate.num = 3,
#     threads = 8,
#     database = "metlinDatabase0.0.2"
#   )
# 
# param5 <-
#   identify_metabolites_params(
#     ms1.match.ppm = 25,
#     rt.match.tol = 1000000,
#     polarity = "negative",
#     ce = "all",
#     column = "hilic",
#     total.score.tol = 0.5,
#     candidate.num = 3,
#     threads = 8,
#     database = "monaDatabase0.0.2"
#   )
# 
# param6 <-
#   identify_metabolites_params(
#     ms1.match.ppm = 25,
#     rt.match.tol = 1000000,
#     polarity = "negative",
#     ce = "all",
#     column = "hilic",
#     total.score.tol = 0.5,
#     candidate.num = 3,
#     threads = 8,
#     database = "nistDatabase0.0.2"
#   )
# 
# param7 <-
#   identify_metabolites_params(
#     ms1.match.ppm = 25,
#     rt.match.tol = 1000000,
#     polarity = "negative",
#     ce = "all",
#     column = "hilic",
#     total.score.tol = 0.5,
#     candidate.num = 3,
#     threads = 8,
#     database = "orbitrapDatabase0.0.1"
#   )
# 
# ##level 3
# param8 <-
#   identify_metabolites_params(
#     ms1.match.ppm = 25,
#     rt.match.tol = 1000000,
#     polarity = "negative",
#     ce = "all",
#     column = "hilic",
#     total.score.tol = 0.5,
#     candidate.num = 3,
#     threads = 8,
#     database = "hmdbMS1Database0.0.1"
#   )
# 
# 
# 
# result_hilic_neg25 <- identify_metabolite_all(
#   ms1.data = "peak_table.csv",
#   ms2.data = c("QC1_MSMS_NCE25.mzXML"),
#   parameter.list = c(param1, param2, param3, param4, param5, param6, param7, param8),
#   path = "."
# )
# 
# save(result_hilic_neg25, file = "result_hilic_neg25")
# 
# 
# result_hilic_neg50 <- identify_metabolite_all(
#   ms1.data = "peak_table.csv",
#   ms2.data = c("QC2_MSMS_NCE35.mzXML"),
#   parameter.list = c(param1, param2, param3, param4, param5, param6, param7),
#   path = "."
# )
# 
# save(result_hilic_neg50, file = "result_hilic_neg50")
# 
# 
# annotation_table1 <-
#   get_identification_table(result_hilic_neg25[[1]], 
#                            result_hilic_neg50[[1]], 
#                            type = "new", 
#                            candidate.num = 1)
# 
# annotation_table1 = 
#   annotation_table1 %>% 
#   dplyr::filter(!is.na(Compound.name))
# 
# 
# annotation_table2 <-
#   get_identification_table(result_hilic_neg25[[2]], 
#                            result_hilic_neg25[[3]], 
#                            result_hilic_neg25[[4]], 
#                            result_hilic_neg25[[5]], 
#                            result_hilic_neg25[[6]], 
#                            result_hilic_neg25[[7]], 
#                            result_hilic_neg50[[2]], 
#                            result_hilic_neg50[[3]], 
#                            result_hilic_neg50[[4]], 
#                            result_hilic_neg50[[5]], 
#                            result_hilic_neg50[[6]], 
#                            result_hilic_neg50[[7]], 
#                            type = "new", 
#                            candidate.num = 1)
# 
# annotation_table2 = 
#   annotation_table2 %>% 
#   dplyr::filter(!is.na(Compound.name)) %>% 
#   dplyr::filter(!name %in% annotation_table1$name)
# 
# 
# annotation_table3 <-
#   get_identification_table(result_hilic_neg25[[8]], 
#                            type = "new", 
#                            candidate.num = 1)
# 
# annotation_table3 = 
#   annotation_table3 %>% 
#   dplyr::filter(!is.na(Compound.name)) %>% 
#   dplyr::filter(!name %in% annotation_table1$name) %>% 
#   dplyr::filter(!name %in% annotation_table2$name)
# 
# annotation_table1 <-
#   data.frame(annotation_table1,
#              Level = 1,
#              stringsAsFactors = FALSE)
# 
# annotation_table2 <-
#   data.frame(annotation_table2,
#              Level = 2,
#              stringsAsFactors = FALSE)
# 
# annotation_table3 <-
#   data.frame(annotation_table3,
#              Level = 3,
#              stringsAsFactors = FALSE)
# 
# annotation_table = 
#   rbind(annotation_table1,
#         annotation_table2) %>%
#   dplyr::full_join(annotation_table3, by = intersect(colnames(annotation_table1),
#                                                      colnames(annotation_table3)))
# 
# write.csv(annotation_table, file = "annotation_table.csv", row.names = FALSE)



setwd(r4projects::get_project_wd())
setwd("3-data_analysis/metabolomics_data/MSMS Annotation/")
rplc_pos = readr::read_csv("RPLC pos/annotation_table.csv")
rplc_neg = readr::read_csv("RPLC neg/annotation_table.csv")
hilic_pos = readr::read_csv("HILIC pos/annotation_table.csv")
hilic_neg = readr::read_csv("HILIC neg/annotation_table.csv")


annotation_result =
  rbind(
    data.frame(rplc_pos, class = "rplc_pos") %>% 
      dplyr::mutate(name = paste("rplc_pos", name, sep = "_")),  
    
    data.frame(rplc_neg, class = "rplc_neg") %>% 
      dplyr::mutate(name = paste("rplc_neg", name, sep = "_")),
    
    data.frame(hilic_pos, class = "hilic_pos") %>% 
      dplyr::mutate(name = paste("hilic_pos", name, sep = "_")),
    
    data.frame(hilic_neg, class = "hilic_neg") %>% 
      dplyr::mutate(name = paste("hilic_neg", name, sep = "_"))
  )

# ###remove duplicated metabolites
# library(plyr)
# annotation_result1 = 
# annotation_result %>% 
# plyr::dlply(.variables = .(Compound.name)) %>% 
#   purrr::map(function(x){
#     if(nrow(x) == 1){
#       return(x)
#     }
#     
#     if(all(x$Level == 3)){
#       return(x %>% 
#         dplyr::filter(mz.error == min(mz.error)))
#     }
#     
#     x = 
#       x %>% 
#       dplyr::filter(Level != 3)
#     
#     if(any(unique(x$Level)) == 1){
#       x = 
#         x %>% 
#         dplyr::filter(Level == 1) %>% 
#         dplyr::filter(SS == max(SS)) %>% 
#         dplyr::filter(Total.score == max(Total.score)) %>% 
#         dplyr::filter(mz.match.score == max(mz.match.score))
#       return(x)
#     }
#     
#     x = 
#       x %>% 
#       dplyr::filter(SS == max(SS)) %>% 
#       dplyr::filter(Total.score == max(Total.score)) %>% 
#       dplyr::filter(mz.match.score == max(mz.match.score))
#     
#     return(x)
#   })


write.csv(annotation_result, file = "annotation_result.csv", row.names = FALSE)
  





  
