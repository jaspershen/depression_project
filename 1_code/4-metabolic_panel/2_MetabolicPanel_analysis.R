no_function()

setwd(r4projects::get_project_wd())

rm(list = ls())

source("1-code/tools.R")

setwd("3_data_analysis/Metabolic_Panel/")

library(samr)
library(zoo)
library(qvalue)
library(tidyverse)

#have manually removed those patients with both depression and safety missing, 
#and duplicate time point entries
metabolic = read.csv("MetabolicPanel_psych_data.csv") 

#adjust for the effect of another variable
get_adjsuted_y <- function(x, y){
  fit = lm(y ~ x)
  correct_coef = fit$coefficients[2]
  y_adjust = y - x * correct_coef
  return(y_adjust)
}

correct_by_CHEX4 <- function(df, start_ind, end_ind){
  var = df[,start_ind:end_ind]
  chex4 = df$CHEX.4
  var_adjusted = sapply(seq(dim(var)[2]), function(i) get_adjsuted_y(chex4, var[,i]))
  var_new = df
  var_new[,start_ind:end_ind] = var_adjusted
  return(var_new)
}

#correct metabolic panel measurements by CHEX4
metabolic[,8:dim(metabolic)[2]] = scale(metabolic[,8:dim(metabolic)[2]])
#plot(metabolic$CHEX.4, metabolic$GLP.1.total)
metabolic_adjusted = correct_by_CHEX4(metabolic, 12, dim(metabolic)[2])
#plot(metabolic_adjusted$CHEX.4, metabolic_adjusted$GLP.1.total)

save(metabolic_adjusted, metabolic, file="adjusted_metabolicPanel.RData")
load("adjusted_metabolicPanel.RData")

#remove T1 for patients 19 and 34, have mannually checked, no these two points

### PART I ###
#calculate slope of x and y vs. time for each patient
get_slopes <- function(df, start_ind, end_ind, y){
  y_change = data.frame()
  x_change = data.frame()
  y_x_change = data.frame()
  
  for (i in unique(df$id)){
    row_ind = which(df$id == i)
    y_i = y[row_ind]
    #remove the rows with NA in the y's
    y_i_ind = which(!is.na(y_i))
    y_i_fil = y_i[y_i_ind]
    
    #only calculate the slopes vs time if we have more than one y
    if (length(unique(y_i_fil)) != 1){
      #extract time
      Ts = df$Time[row_ind][y_i_ind]
      time = sapply(seq(length(Ts)), function(i) as.double(substr(Ts[i],2,2)))
      y_slope <- cor(time, y_i_fil, method="pearson")
      
      x = df[,start_ind:end_ind][row_ind,][y_i_ind,]
      x_slope = as.matrix(sapply(seq(dim(x)[2]), function(i) cor(x[,i], time, method="pearson")))
      y_x_slope = as.matrix(sapply(seq(dim(x)[2]), function(i) cor(x[,i], y_i_fil, method="pearson")))
      rownames(x_slope) = colnames(x)
      rownames(y_x_slope) = colnames(x)
      
      y_change = rbind(y_change, y_slope)
      x_change = rbind(x_change, t(x_slope))
      y_x_change = rbind(y_x_change, t(y_x_slope))
    }
  }
  return(list(y_change = y_change, x_change=x_change, y_x_change=y_x_change))
}

#how the changes of y and x are associated (permute across patients)
metabolic_depression_change_all = get_slopes(metabolic_adjusted, 12, dim(metabolic_adjusted)[2], metabolic_adjusted$bdi_total)
metabolic_safety_change_all = get_slopes(metabolic_adjusted, 12, dim(metabolic_adjusted)[2], metabolic_adjusted$safe)

samfit_dep = SAM(t(data.matrix(metabolic_depression_change_all$x_change)), 
                 y=data.matrix(metabolic_depression_change_all$y_change), 
              resp.type = "Quantitative",
              genenames = colnames(metabolic_depression_change_all$x_change),
              random.seed = 6)

sink('./results_0108_metabolic/depression/metabolic_depression_slope.txt', append=TRUE)
print(samfit_dep)
sink()

samfit_safe = SAM(t(data.matrix(metabolic_safety_change_all$x_change)), 
                  y=data.matrix(metabolic_safety_change_all$y_change), 
              resp.type ="Quantitative",
              genenames = colnames(metabolic_safety_change_all$x_change),
              random.seed=6)

sink('./results_0108_metabolic/safety/metabolic_safety_slope.txt', append=TRUE)
print(samfit_safe)
sink()

### PART II ###
#how y and x are associated
#permute within patients, permute y, and get new slopes, and see if the permuted slopes are more significant
#than the true slopes

insert.value <- function(vec, newval, pos) {
  if (pos == 1)
    return(c(newval, vec))
  lvec <- length(vec)
  if (pos > lvec)
    return(c(vec, newval))
  return(c(vec[1:pos - 1], newval, vec[pos:lvec]))
}

permute <- function(elem) {
  # generates all perms of the vector elem
  if (!missing(elem)) {
    if (length(elem) == 2)
      return(matrix(c(elem, elem[2], elem[1]), nrow = 2))
    last.matrix <- permute(elem[-1])
    dim.last <- dim(last.matrix)
    new.matrix <- matrix(0, nrow = dim.last[1] * (dim.last[2] +
                                                    1), ncol = dim.last[2] + 1)
    for (row in 1:(dim.last[1])) {
      for (col in 1:(dim.last[2] + 1)) 
        new.matrix[row + (col - 1) * dim.last[1], ] <- insert.value(last.matrix[row,], elem[1], col)
    }
    return(new.matrix)
  }
  else cat("Usage: permute(elem)\n\twhere elem is a vector\n")
}

#true average association
rho_true_depression = colMeans(metabolic_depression_change_all$y_x_change)
rho_true_safety = colMeans(metabolic_safety_change_all$y_x_change)

#permute within patient and calculate average association
calc_perm_corr <- function(df, i, y, start_ind, end_ind){
  #print(i)
  row_ind = which(df[,1] == i)
  y_i = y[row_ind]
  y_i_ind = which(!is.na(y_i))
  y_i_fil = y_i[y_i_ind]
  Ts = df[,2][row_ind][y_i_ind]
  time = sapply(seq(length(Ts)), function(i) as.double(substr(Ts[i],2,2)))
  x = t(df[,start_ind:end_ind][row_ind,][y_i_ind,])
  
  if (length(unique(y_i_fil)) != 1){
    all_perms = permute(y_i_fil)
    y_perms = all_perms[sample(nrow(all_perms), 1),]
    x_y_slope = as.matrix(sapply(seq(dim(x)[1]), function(i) cor(x[i,], y_perms, method="pearson")))
    rownames(x_y_slope) = rownames(x)
  } else {
    #make up an invalid number
    x_y_slope = rep(NA, dim(x)[1])
  }
  return(x_y_slope)
}

do_perm <- function(df, start_ind, end_ind, type="depression", nperms=1000){
  res = data.frame()
  for (j in seq(nperms)){
    print(j)
    each_perm_res = data.frame()
    
    if (type == "depression"){
      y = df$bdi_total
    } else if (type == "safety"){
      y = df$safe
    }
    
    perm_corr = t(sapply(unique(df[,1]), function(i) calc_perm_corr(df, i, y, start_ind, end_ind)))
    res = rbind(res, t(colMeans(perm_corr, na.rm=TRUE)))
  }
  return(res)
}

perm_dep = do_perm(metabolic_adjusted, 12, dim(metabolic_adjusted)[2], type="depression", nperms=1000)
perm_safety = do_perm(metabolic_adjusted, 12, dim(metabolic_adjusted)[2], type="safety", nperms=1000)
colnames(perm_dep) = colnames(metabolic_depression_change_all$y_x_change)
colnames(perm_safety) = colnames(metabolic_safety_change_all$y_x_change)
save(perm_dep, perm_safety, file="permutation_metabolic.RData")
load("permutation_metabolic.RData")

get_result <- function(true_rho, perm_rho, nperms, qvalue_threshold=0.3, pvalue_threshold=0.1){
  pos_ind = which(true_rho >= 0)
  neg_ind = which(true_rho < 0)
  p_value = sapply(seq(true_rho), function(i) sum(abs(true_rho[i]) <= abs(perm_rho[,i])) / nperms)
  q_value = qvalue(p = p_value, pi0=1)$qvalues
  
  q_value_pos = q_value[pos_ind]
  q_value_pos_gene = colnames(perm_rho)[pos_ind]
  q_value_neg = q_value[neg_ind]
  q_value_neg_gene = colnames(perm_rho)[neg_ind]
  
  q_value_ind_pos_ord = order(q_value_pos)
  q_value_pos_ord = q_value_pos[q_value_ind_pos_ord]
  gene_pos_ord = q_value_pos_gene[q_value_ind_pos_ord]
  p_value_pos_ord = p_value[pos_ind][q_value_ind_pos_ord]
  gene_up_qvalue_sig = data.frame(gene_pos_ord, q_value_pos_ord, p_value_pos_ord)[which(q_value_pos_ord<qvalue_threshold
                                                                                        & p_value_pos_ord<pvalue_threshold),] 
  colnames(gene_up_qvalue_sig) = c("Variables", "qvalue", "unadjusted pvalue")
  
  q_value_ind_neg_ord = order(q_value_neg)
  q_value_neg_ord = q_value_neg[q_value_ind_neg_ord]
  gene_neg_ord = q_value_neg_gene[q_value_ind_neg_ord]
  p_value_neg_ord = p_value[neg_ind][q_value_ind_neg_ord]
  gene_down_qvalue_sig = data.frame(gene_neg_ord, q_value_neg_ord, p_value_neg_ord)[which(q_value_neg_ord<qvalue_threshold
                                                                                          & p_value_neg_ord<pvalue_threshold),] 
  colnames(gene_down_qvalue_sig) = c("Variables", "qvalue", "unadjusted pvalue")
  return(list(up=gene_up_qvalue_sig, down=gene_down_qvalue_sig))
}

dep_res = get_result(rho_true_depression, perm_dep, 1000, qvalue_threshold=0.3, pvalue_threshold=0.1)
safety_res = get_result(rho_true_safety, perm_safety, 1000, qvalue_threshold=0.3, pvalue_threshold=0.1)

options(max.print=1000000)
sink('./results_0108_metabolic/depression/metabolic_depression_association.txt', append=TRUE)
print("How variables are associated with depression? (permute within patients)")
print("Positively associated:")
print(dep_res$up)
print("")
print("Negatively associated:")
print(dep_res$down)
sink()

options(max.print=1000000)
sink('./results_0108_metabolic/safety/metabolic_safety_association.txt', append=TRUE)
print("How variables are associated with perceived safety? (permute within patients)")
print("Positively associated:")
print(safety_res$up)
print("")
print("Negatively associated:")
print(safety_res$down)
sink()





