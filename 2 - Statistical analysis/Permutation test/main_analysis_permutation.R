########################################################################################################################################################
############################          Code used in Restricted mean survival time in cluster randomized trials with          ############################  
############################    a small number of clusters: Improving variance estimation of the pseudo-values regression   ############################  
########################################################################################################################################################

#--------------------------------------------------------------------------------------------------------------------------------------------------------# 
#----------------------------------------------------- 2 - Analysis of the datasets - Permutation test --------------------------------------------------#  
#--------------------------------------------------------------------------------------------------------------------------------------------------------# 

# Packages
##########################################################################################################################################################
library(doParallel)
library(doRNG)

# Functions needed 
#######################################################################################
setwd("~/your/path/to/Rfiles/") # set directory where you saved the R files with the necessary functions
source("sim_analysis_permutation.R")
source("ci.R")
source("ci_permutation.R")
source("initialisation.R")
source("initialisation_ci.R")
source("allocation.R")
source("update_bound.R")


# Parameters
###########################################################################################################
#Parameter of the scenario
table_parameter <- data.frame("D" = 1000,             #number of simulated dataset
                              "K" = 10,            #total number of clusters (should be pair)
                              "m" = 80,            #mean cluster size
                              "cv" = 0.6,          #coefficient of variatin of the cluster sizes
                              "lambda" = 0.000016, #scale parameter of the Weibull distribution 
                              "rho" = 2,           #shape parameter of the Weibull distribution
                              "tau" = 0.1,         #Kendall's tau
                              "HR"= 0.5,    #intervention effect (beta)
                              "censoring" = 0.2,   #censoring rate (between 0 and 1)
                              "seed" = 1598)       #seed

# Horizon time (t*)
t_star <- 365 


# Analysis
###########################################################################################################
#set directory where the folder with the simulated datasets have been saved and where the estimations will be saved
setwd("~/your/path/to/Simulated_datasets/") 

###Create file to save data 
name.file <- paste("HR=", table_parameter[ ,"HR"],
                   "-tau=", table_parameter[ ,"tau"], 
                   "-K=", table_parameter[ ,"K"], 
                   "-m=", table_parameter[ ,"m"], 
                   "-cv=", table_parameter[ , "cv"],
                   "-censoring=",table_parameter[ ,"censoring"], 
                   sep = "")   


write.table(data.frame("d" = "d",
                       "K" = "K",
                       "m" = "m",
                       "cv" = "cv",
                       "HR" = "HR",
                       "tau" = "tau",
                       "censoring"="censoring",
                       "matrix" = "matrix",
                       "ci_low" = "ci_low",
                       "ci_up" = "ci_up"), 
            file= paste(name.file, ".txt", sep = ""), 
            sep = " ", dec = ".",
            col.names = FALSE, row.names = FALSE)

###Parallelisation
registerDoParallel(cores = 8) #set the number of cores
set.seed(table_parameter[ ,"seed"])
res <- foreach(d = 1: table_parameter[ ,"D"],
             .packages = c("pseudo", "gee", "survRM2", "foreach")) %dorng% sim_analysis_permutation(K = table_parameter[ ,"K"], 
                                                                                                    m = table_parameter[ ,"m"], 
                    									                                                              cv = table_parameter[ ,"cv"],
                                                                                                    lambda = table_parameter[ ,"lambda"], 
                    									                                                              rho = table_parameter[ ,"rho"], 
                                                                                                    tau  = table_parameter[ ,"tau"], 
                                                                                                    HR = table_parameter[ ,"HR"], 
                                                                                                    censoring = table_parameter[ ,"censoring"], 
                                                                                                    t_star = t_star, 
                                                                                                    d, 
                                                                                                    name.file) 
