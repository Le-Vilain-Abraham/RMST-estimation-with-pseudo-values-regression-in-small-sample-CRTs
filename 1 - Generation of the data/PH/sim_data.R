##########################################################################################################################################################
# Simulate one dataset
##########################################################################################################################################################

#####  Arguments  #####
# K: number of clusters (int > 0)
# m: mean cluster size (int > 0)
# cv: coefficient of variation of the cluster sizes (float between 0 and 1)
# lambda and rho: parameters of the Weibull distribution (floats > 0)
# tau: Kendall's tau (float > 0)
# HR: hazard ratio (float)
# censoring: censoring rate (float between 0 and 1)
# d: iteration number
# name.file: name of the folder and the file where the dataset and the estimation will be saved

####  Values ####
# one simulated dataset saved in a txt file


sim_data <- function(K, m, cv, lambda, rho, tau, HR, censoring, d, name.file){
  
  ###### Generate one dataset ######
  dataset <- generate_data(K, m, cv, lambda, rho, (1-tau)/(2*tau), log(HR), censoring)
  
  while(min(max(dataset[which(dataset$arm == 1), ]$time), 
            max(dataset[which(dataset$arm == 0), ]$time))<365){
    dataset <- generate_data(K, m, cv, lambda, rho, (1-tau)/(2*tau), log(HR), censoring)
  }

  dataset <- dataset[order(dataset[, "cluster"], decreasing = T), ]
  
  ###### Save dataset ######
  write.csv2(dataset, 
              file= paste(name.file, "/dataset", d, ".csv", sep = ""),  
              row.names = FALSE)

}

