setwd("/Users/daisyding/Dropbox/daisy/depression_analysis")
library(zoo)
library(qvalue)
library(ggplot2)
library(tidyverse)
library(ggpubr)
library(rstatix)

cyt = read.csv("Cytokines1_new_psych_data_no_duplicate.csv") 

#remove T1 for patients 19 and 34
T1_19 = which(cyt[,1] == 19 & cyt[,2] == "T1")
T1_34 = which(cyt[,1] == 34 & cyt[,2] == "T1")
cyt_fil = cyt[-c(T1_19, T1_34),]

#calculate the correlation of rank of psychometrics scores vs. time, for each patient
cal_rank_cor <- function(df, y){
  ### Input ###
  # df: dataframe
  # y: psychometrics score vector
  
  y_change = data.frame()
  
  for (i in unique(df$id)){
    row_ind = which(df$id == i)
    y_i = y[row_ind]
    #remove NA scores
    y_i_ind = which(!is.na(y_i))
    y_i_fil = y_i[y_i_ind]
    
    #only calcualte the correlation if the psychometrics scores for the patient has more than one data point
    if (length(unique(y_i_fil)) != 1){
      #extract time
      Ts = df$Time[row_ind][y_i_ind]
      time = sapply(seq(length(Ts)), function(i) as.double(substr(Ts[i],2,2)))
      time_rank = seq(length(time))
      
      #rank psychometircs scores
      score_rank = order(y_i_fil) 
      #print(time_rank)
      #print(score_rank)
      score_cor <- cor(time_rank, score_rank, method="pearson")
      y_change = rbind(y_change, t(score_cor))
    }
  }
  return(y_change)
}

#depression
dep_rank_cor = cal_rank_cor(cyt_fil, cyt_fil$bdi_total) #want it to be negative
#two-sided t-test
t.test(dep_rank_cor, mu=0, alternative="two.sided")

#plot the correlation scores for all  patients with t-test
dep_rank_cor %>% get_summary_stats(V1, type = "mean_sd")
box_dep = ggboxplot(dep_rank_cor$V1, width = 0.5, add = c("mean", "jitter"), 
  ylab = "Correlation", xlab = FALSE, main="Depression: Rank of Depression vs Time")
test_dep <- dep_rank_cor %>% t_test(V1 ~ 1, mu = 0)
box_dep + 
  labs(subtitle = get_test_label(test_dep, detailed = TRUE))
  
#perceived safety
safety_rank_cor = cal_rank_cor(cyt_fil, cyt_fil$safe) #want it to be positive
t.test(safety_rank_cor, mu=0, alternative="two.sided")

safety_rank_cor %>% get_summary_stats(V1, type = "mean_sd")
box_safe = ggboxplot(safety_rank_cor$V1, width = 0.5, add = c("mean", "jitter"), 
                    ylab = "Correlation", xlab = FALSE, main="Perceived Safety: Rank of Safety vs Time")
test_safe <- safety_rank_cor %>% t_test(V1 ~ 1, mu = 0)
box_safe + 
  labs(subtitle = get_test_label(test_safe, detailed = TRUE))
