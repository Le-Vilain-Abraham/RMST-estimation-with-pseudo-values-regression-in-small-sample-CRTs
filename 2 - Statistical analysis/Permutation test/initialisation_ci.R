########################################################################################
# To initialise the confidence bounds for the search procedure: 
# Estimate one permuted delta_offset from the permutation test H0: theta = theta_hat
########################################################################################

#####  Arguments  #####
# data: dataset to analyse (data.frame)
# matrix: working correlation matrix (character: "independence" for independent matrix)

####  Values ####
#results: estimation of one kappa for a permuted dataset (data.frame)


initialisation_ci <- function(data, matrix) {
  
  #Re-assign the intervention and control group
  dataset <- allocation(data)
  
  #Analysis
  log <- capture.output(suppressMessages(fit_gee_per <- gee(pv ~ arm + offset(arm.obs*theta_hat), 
                                                            data = dataset, 
                                                            id = cluster, 
                                                            family = gaussian, 
                                                            corstr = matrix)))
    
  return(coef(summary(fit_gee_per))["arm","Estimate"])
}


