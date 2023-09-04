########################################################################################
# Permutes intervention allocation
########################################################################################

#####  Arguments  #####
# data: dataset to analyse (data.frame)

####  Values ####
#results: dataset with permuted intervention allocation (data.frame)


allocation <- function(data){
 
  K <- max(data$cluster)
  allocation <- sample(1:K, K/2, replace = FALSE)
  
  #Re-assign the intervention and control group
  data$arm <- 0
  data[data$cluster %in% allocation, ]$arm <- 1
  
  return(data)
  
}
