########################################################################################
# Update bound using Robbins-Monro search procedure
########################################################################################

#####  Arguments  #####
# bound_p: current value of the bound (float)
# theta_hat: point estimate of beta1 (difference in RMST) based on original dataset (float)
# T_obs: observed test statistic (float)
# T_permuted: permuted test statistic (float)
# alpha: significance level of the test (float between 0 and 1)
# p: iteration of the search procedure (interger >0)
# m: magnitude of step
# k: step lenght multiplier
# bound: confidence bound (character: "low" for the lower bound or "up" for the upper bound)

####  Values ####
#results: initialisations for the search procedure (data.frame)


update_bound <- function(bound_p, theta_hat, T_obs, T_permuted, alpha, p, m, k, bound){
  p <- m + p - 1
  
  if(bound=="low"){
  c <- k*(theta_hat-bound_p)
    bound_p1 <- ifelse(T_obs > T_permuted,
                       bound_p + c*(alpha/2)/p,
                       bound_p - c*(1-alpha/2)/p)
  }else{
  c <- k*(bound_p-theta_hat)
    bound_p1 <- ifelse(T_obs < T_permuted,
                       bound_p - c*(alpha/2)/p,
                       bound_p + c*(1-alpha/2)/p)
  }
  
  return(bound_p1)
}
