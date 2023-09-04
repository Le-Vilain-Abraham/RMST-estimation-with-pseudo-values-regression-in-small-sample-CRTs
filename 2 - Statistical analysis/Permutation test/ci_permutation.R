########################################################################################
# Construction of one permutation-based confidence bound (lower or upper bound)
########################################################################################

#####  Arguments  #####
# data: dataset to analyse (data.frame)
# matrix: working correlation matrix (character: "independence" for independent matrix)
# m: magnitude of step
# k: step lenght multiplier
# bound_p: initialisation of the confidence bound (float)
# nperm: number of step of the search procedure (integer > 0)
# theta_hat: point estimate of beta1 (difference in RMST) based on original dataset (float)
# bound: confidence bound (character: "low" for the lower bound or "up" for the upper bound)

####  Values ####
#results: initialisations for the search procedure (data.frame)


ci_permutation <- function(data, matrix, m, k, bound_p, nperm, theta_hat, bound){
 
  init <- bound_p
  vec <- rep(NA, nperm)

  data$bound_p <- bound_p

 for(p in 1:nperm){
    #observé
    log <- capture.output(suppressMessages(fit_gee <- gee(pv ~ arm + offset(arm.obs * bound_p), 
                                                              data = data, 
                                                              id = cluster, 
                                                              family = gaussian, 
                                                              corstr = matrix)))
   
    #permuté
    data_permu <- allocation(data)
    log <- capture.output(suppressMessages(fit_gee_per <- gee(pv ~ arm + offset(arm.obs * bound_p), 
                                                              data = data_permu, 
                                                              id = cluster, 
                                                              family = gaussian, 
                                                              corstr = matrix)))
    
    T_obs <- coef(summary(fit_gee))["arm","Robust z"]
    T_permuted <- coef(summary(fit_gee_per))["arm","Robust z"]
    
    # update using Robbins-Monro step
    bound_p <- update_bound(bound_p, theta_hat, T_obs, T_permuted, alpha=0.05, p, m, k, bound)
    data$bound_p <- bound_p
    vec[p] <- bound_p
    }

  return(list(vec[nperm]))
}
