########################################################################################################################################################
############################          Code used in Restricted mean survival time in cluster randomized trials with          ############################  
############################    a small number of clusters: Improving variance estimation of the pseudo-values regression   ############################  
########################################################################################################################################################

#--------------------------------------------------------------------------------------------------------------------------------------------------------# 
#------------------------------------------ 3 - Estimation of the performance measures - Variance corrections -------------------------------------------#  
#--------------------------------------------------------------------------------------------------------------------------------------------------------# 

#library
library(dplyr)

# R functions 
##########################################################################################################################################################
setwd("~/your/path/to/Rfiles/") # set directory where you saved the R files with the necessary functions for the simulation
source("pm_estimation_variance.R")
source("performance_measures.R")
source("true_rmst_difference.R")
source("survival_function.R")

# Estimation of the performance measures 
##########################################################################################################################################################
# set the directory and the name of the file where all the estimations (difference in RMST, variance, 95% confidence interval) were saved in step 2 
# for the variance corrections
dataset <- read.table("~/your/path/to/analysis/name_of_your_file.txt", sep = " ", dec = ".",  header=T)

pm <- pm_estimation_variance(dataset)       
