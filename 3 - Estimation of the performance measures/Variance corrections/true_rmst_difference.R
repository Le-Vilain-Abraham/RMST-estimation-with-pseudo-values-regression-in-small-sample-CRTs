##########################################################################################################################################################
# Computation of the true difference in RMST
##########################################################################################################################################################

#####  Arguments  #####
# t: horizon time (integer > 0)
# step: step of the right Riemann sum (float > 0)
# lambda and rho: parameters of the Weibull distribution (float > 0)
# beta: treatment effect (float)
# gamma: parameter of the gamma distribution of the frailty terms (float >0)

####  Values ####
# True difference in RMSTs (float)


true_rmst_difference <- function(t, step, lambda, rho, beta, gamma) {
  
  return((sum(survival_function(seq(1, t, by=step),lambda, rho, 1, beta, gamma))
         - sum(survival_function(seq(1, t, by=step),lambda, rho, 0, beta, gamma)))*step)
  
}
