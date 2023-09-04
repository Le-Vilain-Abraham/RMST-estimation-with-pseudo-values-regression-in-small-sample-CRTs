##########################################################################################################################################################
# Estimation of the coverage rate for the permutation-based confidence intervals
##########################################################################################################################################################

#####  Arguments  #####
# dataset: dataset of the 1000 estimations of difference in RMST, its variance and 95% confidence interval (data.frame)

####  Values ####
# Coverage rate and ype I eo/Powe (data.frame)


pm_estimation_permutation <- function(data) {
  
  results <- data.frame()

  results <- performance_measures(365, 
                                    lambda = 0.000016, 
                                    rho = 2, 
                                    beta = log(data$HR[1]), 
                                    data = data, 
                                    gamma = ifelse(data$tau[1]==0, 0, (1-data$tau[1])/(2*data$tau[1])))

  return(cbind(tstar = 365,
                   K = data$K[1],
                   m = data$m[1],
                   cv = data$cv[1],
                   HR = data$HR[1],
                   tau = data$tau[1],
                   censoring = data$censoring[1],
                   matrix = data$matrix[1],
                   results))
  
}
