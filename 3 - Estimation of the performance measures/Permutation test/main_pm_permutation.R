########################################################################################################################################################
############################          Code used in Restricted mean survival time in cluster randomized trials with          ############################  
############################    a small number of clusters: Improving variance estimation of the pseudo-values regression   ############################  
########################################################################################################################################################

#--------------------------------------------------------------------------------------------------------------------------------------------------------# 
#---------------------------------------------- 3 - Estimation of the performance measures - Permutation ------------------------------------------------#  
#--------------------------------------------------------------------------------------------------------------------------------------------------------# 

# R functions 
##########################################################################################################################################################
setwd("~/your/path/to/Rfiles/") # set directory where you saved the R files with the necessary functions for the simulation
source("pm_estimation_permutation.R")
source("performance_measures.R")
source("survival_function.R")
source("true_rmst_difference.R")


# Estimation of the coverage rate
##########################################################################################################################################################
# set the directory and the name of the file where all the estimations of the confidence intervals were saved in step 2 for the permutation test
dataset <- read.table("~/your/path/to/analysis/name_of_your_file.txt", sep = " ", dec = ".",  header=T)

pm <- pm_estimation_permutation(dataset)       
