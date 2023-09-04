##########################################################################################################################################################
# Estimation of the performance measures of the 45 aiances coecion, 2 disiuion, 2 coelaion maices
##########################################################################################################################################################

#####  Arguments  #####
# data: dataset of the 1000 estimations of difference in RMST, its variance and 95% confidence interval (data.frame)

####  Values ####
# relative bias (RE), relative error (RE), the coverage rate (coverage), type I error rate (for scenario where HR=1) (rejection.rate) and the number 
# of simulation iterations which converged (D) for the 5 methods (data.frame)


pm_estimation_variance  <- function(data) {

  #data pour lesquels ça a convergé
  data <- data[which((data$delta.rmst==data$var)==FALSE),]  
  
  data$meth <-  paste(data$method_correction, data$distribution, data$matrix, sep = "_")
  
  results <- data.frame()
  for(meth in unique(paste(data$method_correction, data$distribution,data$matrix, sep = "_"))) {
    results <-rbind(results, cbind(method_total = meth,
                                   method = data[which(data$meth==meth),]$method[1],
                                   distribution = data[which(data$meth==meth),]$distribution[1],
                                   matrix = data[which(data$meth==meth),]$matrix[1],
                                   performance_measures(365, 
                                                         lambda = 0.000016, 
                                                         rho = 2, 
                                                         beta = log(data$HR[1]), 
                                                         data = data[which(data$meth==meth),], 
                                                         gamma = ifelse(data$tau[1]==0, 0, (1-data$tau[1])/(2*data$tau[1])))))
  }

  
  results <- cbind(tstar = 365,
                    K = data$K[1],
                    m = data$m[1],
                    cv = data$cv[1],
                    HR = data$HR[1],
                    tau = data$tau[1],
                    censoring = data$censoring[1],
                    results)
  return(results)
}
