##########################################################################################################################################################
# True survival function
###########################################################################################################################################################

#####  Arguments  #####
# t: time (float >=0)
# lambda and rho: parameters of the Weibull distribution (float > 0)
# Z: arm (0 for control, 1 for intervention)
# beta: treatment effect (float)
# gamma: parameter of the gamma distribution of the frailty terms (float >0)
# T_change: change point

####  Values ####
# survival at time t


survival_function <- function(t, lambda, rho, Z, beta, gamma, T_change){
  if(gamma == 0) {
    exp(-(lambda*T_change^rho + lambda*exp(Z*beta)*(t^rho-T_change^rho)))
    
  } else {
    1/(1+(lambda*T_change^rho + lambda*exp(Z*beta)*(t^rho-T_change^rho))/gamma)^gamma
  }
}
