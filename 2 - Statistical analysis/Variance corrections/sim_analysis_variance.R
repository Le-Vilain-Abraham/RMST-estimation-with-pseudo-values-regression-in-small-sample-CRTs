#########################################################################################
# Simulate and analyse one dataset
########################################################################################
#
#input: K: number of clusters (int > 0)
#       m: cluster size or mean cluster size (int > 0)
#       lambda and rho: parameters of the Weibull distribution (float > 0)
#       tau: Kendall's tau (float > 0)
#       HR: hazard ratio (float)
#       censoring: censoring rate (float between 0 and 1)
#       t_star: horizon time (float > 0)
#
#ouput: estimations for all the methods for one simulated dataset (data.frame)


sim_analysis_variance <- function(K, m, cv, lambda, rho, tau, HR, censoring, t_star, d, name.file){
  
  ###### Load one dataset ######
  # set directory where the simulated datasets have been saved
  dataset <- read.table(file = paste("~/your/path/to/Simulated_datasets/",name.file, "/dataset", d, ".csv", sep = ""),
                        header = T, sep = ";", dec = ",")
  dataset <- dataset[order(dataset[,"arm"],decreasing=T), ]
  
  ###### Analyses using all the methods #######
  results <-RMST(dataset, t_star) 
  
  ###### Save results ######
  write.table(cbind(d = d,
 		                K = K,
                    m = m,
                    cv = cv,
                    HR = HR,
                    tau = tau,
                    censoring = censoring,
                    results), 
              file= paste(name.file, ".txt", sep = ""), 
              append = TRUE, 
              sep = " ", dec = ".", 
              col.names = FALSE, row.names = FALSE)
}

