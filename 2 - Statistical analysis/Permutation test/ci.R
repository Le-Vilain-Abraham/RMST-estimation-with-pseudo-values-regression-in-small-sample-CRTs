########################################################################################
# Permutation-based confidence interval construction for one dataset
########################################################################################

#####  Arguments  #####
# data: dataset to analyse (data.frame)
# matrix: working correlation matrix (character: "independence" for independent matrix)

####  Values ####
#results: estimations of the permutation-based confidence interval (data.frame)


ci <- function(data, matrix) {

  log <- capture.output(suppressMessages(fit_gee <- gee(pv ~ arm, 
                                                            data = data, 
                                                            id = cluster, 
                                                            family = gaussian, 
                                                            corstr = matrix)))
  

    ####### Initialization ####### 
    theta_hat <-  coef(summary(fit_gee))["arm","Estimate"]
    init <- initialisation(data, matrix, theta_hat)
    
    ####### Confidence interval ####### 
    results_low <- ci_permutation(data, matrix, init$m, init$k, init$low, nperm=5000, theta_hat, bound="low")
    results_up <- ci_permutation(data, matrix, init$m, init$k, init$up, nperm=5000, theta_hat, bound="up")

  
    return(data.frame(ci_low = results_low[[1]], ci_up = results_up[[1]]))
}
