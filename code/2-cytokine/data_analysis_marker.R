no_function()

setwd(r4projects::get_project_wd())
library(tidyverse)
library(plyr)
source("code/tools.R")

setwd("data_analysis/Cytokines/")

###load data
load("data_preparation/expression_data")
load("data_preparation/sample_info")
load("data_preparation/variable_info")

##read marker
###depression

###x change is the time vs 

depression_association_pos = 
readr::read_table(file = "results_0328_cytokines/cytokine_depression_association_pos.txt")

depression_association_pos = 
  depression_association_pos[,-1]

depression_association_neg = NULL

depression_slope_up = 
  readr::read_table(file = "results_0328_cytokines/cytokine_depression_slope_up.txt")

depression_slope_up = 
  depression_slope_up[,-1]

depression_slope_down = NULL


###safety
safety_association_pos = 
  readr::read_table2(file = "results_0328_cytokines/cytokine_safety_association_pos.txt")

safety_association_pos = 
  safety_association_pos[,-1]

safety_association_neg = 
  readr::read_table2(file = "results_0328_cytokines/cytokine_safety_association_neg.txt")

safety_association_neg = 
  safety_association_neg[,-1]

safety_slope_up = 
  readr::read_table(file = "results_0328_cytokines/cytokine_safety_slope_up.txt")

safety_slope_up = 
  safety_slope_up[,-1]

safety_slope_down = NULL


dim(expression_data)

#####association is the makers changes with time
#####slope is the accelerated

depression_association_pos
depression_association_neg
safety_association_pos
safety_association_neg

dir.create("marker")

save(depression_association_pos, file = "marker/depression_association_pos")
save(depression_association_neg, file = "marker/depression_association_neg")
save(safety_association_pos, file = "marker/safety_association_pos")
save(safety_association_neg, file = "marker/safety_association_neg")

save(depression_slope_up, file = "marker/depression_slope_up")
save(depression_slope_down, file = "marker/depression_slope_down")
save(safety_slope_up, file = "marker/safety_slope_up")
save(safety_slope_down, file = "marker/safety_slope_down")
























