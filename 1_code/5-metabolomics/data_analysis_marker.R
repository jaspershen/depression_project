no_function()

setwd(r4projects::get_project_wd())
library(tidyverse)
library(plyr)
source("1-code/tools.R")

setwd("3_data_analysis/metabolomics_data/")

# ###load data
# load("data_preparation/expression_data")
# load("data_preparation/sample_info")
# load("data_preparation/variable_info")

##read marker
###depression
###x change is the time vs 
depression_association_pos = 
  try(
    readr::read_table(file = "results_0307_metabolites/depression/metabolite_depression_association_pos.txt")
    ,silent = TRUE
  )

if(class(depression_association_pos)[1] == "try-error"){
  depression_association_pos = NULL
}else{
  depression_association_pos = 
    depression_association_pos[,-1]
}


depression_association_neg = 
  try( readr::read_table(file = "results_0307_metabolites/depression/metabolite_depression_association_neg.txt"),
       silent = TRUE)

if(class(depression_association_neg)[1] == "try-error"){
  depression_association_neg = NULL
}else{
  depression_association_neg = 
    depression_association_neg[,-1]
}
 
depression_slope_up = 
  try(readr::read_table2(file = "results_0307_metabolites/depression/metabolite_depression_slope_up.txt"), 
      silent = TRUE)
  
if(class(depression_slope_up)[1] == "try-error"){
  depression_slope_up = NULL
}else{
  depression_slope_up = 
    depression_slope_up[,-1]
}

depression_slope_up = 
depression_slope_up %>% 
  purrr::map(function(x){
    x %>% 
    stringr::str_replace_all(pattern = '\"', "")
  }) %>% 
  do.call(cbind, .) %>% 
  as.data.frame()

depression_slope_down = 
  try(readr::read_table(file = "results_0307_metabolites/depression/metabolite_depression_slope_down.txt"), 
      silent = TRUE)

if(class(depression_slope_down)[1] == "try-error"){
  depression_slope_down = NULL
}else{
  depression_slope_down = 
    depression_slope_down[,-1]
}

depression_slope_down = 
  depression_slope_down %>% 
  purrr::map(function(x){
    x %>% 
      stringr::str_replace_all(pattern = '\"', "")
  }) %>% 
  do.call(cbind, .) %>% 
  as.data.frame()

###safety
safety_association_pos = 
  try(
    readr::read_table(file = "results_0307_metabolites/safety/metabolite_safety_association_pos.txt")
    ,silent = TRUE
  )

if(class(safety_association_pos)[1] == "try-error"){
  safety_association_pos = NULL
}else{
  safety_association_pos = 
    safety_association_pos[,-1]
}

safety_association_neg = 
  try( readr::read_table(file = "results_0307_metabolites/safety/metabolite_safety_association_neg.txt"),
       silent = TRUE)

if(class(safety_association_neg)[1] == "try-error"){
  safety_association_neg = NULL
}else{
  safety_association_neg = 
    safety_association_neg[,-1]
}

safety_slope_up = 
  try(readr::read_table(file = "results_0307_metabolites/safety/metabolite_safety_slope_up.txt"), 
      silent = TRUE)

if(class(safety_slope_up)[1] == "try-error"){
  safety_slope_up = NULL
}else{
  safety_slope_up = 
    safety_slope_up[,-1]
}

safety_slope_up = 
  safety_slope_up %>% 
  purrr::map(function(x){
    x %>% 
      stringr::str_replace_all(pattern = '\"', "")
  }) %>% 
  do.call(cbind, .) %>% 
  as.data.frame()


safety_slope_down = 
  try(readr::read_table(file = "results_0307_metabolites/safety/metabolite_safety_slope_down.txt"), 
      silent = TRUE)

if(class(safety_slope_down)[1] == "try-error"){
  safety_slope_down = NULL
}else{
  safety_slope_down = 
    safety_slope_down[,-1]
}

depression_association_pos
depression_association_neg
depression_slope_up
depression_slope_down

safety_association_pos
safety_association_neg
safety_slope_up
safety_slope_down


#####association is the makers changes with time
#####slope is the accelerated

dir.create("marker")

save(depression_association_pos, file = "marker/depression_association_pos")
save(depression_association_neg, file = "marker/depression_association_neg")
save(safety_association_pos, file = "marker/safety_association_pos")
save(safety_association_neg, file = "marker/safety_association_neg")

save(depression_slope_up, file = "marker/depression_slope_up")
save(depression_slope_down, file = "marker/depression_slope_down")
save(safety_slope_up, file = "marker/safety_slope_up")
save(safety_slope_down, file = "marker/safety_slope_down")
























