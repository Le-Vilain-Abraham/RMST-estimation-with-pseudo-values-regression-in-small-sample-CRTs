########################################################################################################################################################
############################          Code used in Restricted mean survival time in cluster randomized trials with          ############################  
############################    a small number of clusters: Improving variance estimation of the pseudo-values regression   ############################  
########################################################################################################################################################

#--------------------------------------------------------------------------------------------------------------------------------------------------------# 
#------------------------------------------------------------ 1 - Generation of the datasets ------------------------------------------------------------#  
#--------------------------------------------------------------------------------------------------------------------------------------------------------# 

# Packages
##########################################################################################################################################################
library(doParallel)
library(doRNG)

# Functions needed 
##########################################################################################################################################################
setwd("~/your/path/to/Rfiles/") # set directory where you saved the R files with the necessary functions
source("sim_data.R")
source("generate_data.R")
source("generate_cluster.R")

# Parameters
##########################################################################################################################################################
#Parameter for the data generation
table_parameter <- data.frame("D" = 1000,          #number of simulated dataset
                              "K" = 10,            #total number of clusters (should be pair)
                              "m" = 80,            #mean cluster size
                              "cv" = 0.6,          #coefficient of variatin of the cluster sizes
                              "lambda" = 0.000016, #scale parameter of the Weibull distribution 
                              "rho" = 2,           #shape parameter of the Weibull distribution
                              "tau" = 0.1,         #Kendall's tau
                              "HR"= 0.5,           #intervention effect
                              "censoring" = 0.2,   #censoring rate (between 0 and 1)
                              "seed" = 1598)       #seed


# Generation of the data and analysis 
##########################################################################################################################################################
### Create a folder to save the datasets 
setwd("~/your/path/to/Simulated_datasets/") # set directory where the simulated datasets will be saved

name.file <- paste("HR=", table_parameter[ ,"HR"],
                   "-tau=", table_parameter[ ,"tau"], 
                   "-K=", table_parameter[ ,"K"], 
                   "-m=", table_parameter[ ,"m"], 
                   "-cv=", table_parameter[ , "cv"],
                   "-censoring=",table_parameter[ ,"censoring"], 
                   sep = "")    

dir.create(name.file)

###Parallelisation
registerDoParallel(cores = 8) #set the number of cores
set.seed(table_parameter[,"seed"])
res <- foreach(d = 1: table_parameter[ ,"D"],
               .packages = c("pseudo", "gee", "foreach", "geesmv", "matrixcalc")) %dorng% sim_data(K = table_parameter[ ,"K"], 
                                                                                                   m = table_parameter[ ,"m"],
                                                                                                   cv = table_parameter[ ,"cv"],
                                                                                                   lambda = table_parameter[ ,"lambda"], 
                                                                                                   rho = table_parameter[ ,"rho"], 
                                                                                                   tau = table_parameter[ ,"tau"], 
                                                                                                   HR = table_parameter[ ,"HR"], 
                                                                                                   censoring = table_parameter[ ,"censoring"], 
                                                                                                   d, 
                                                                                                   name.file) 
