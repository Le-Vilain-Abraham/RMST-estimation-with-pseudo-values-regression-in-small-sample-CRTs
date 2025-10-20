##########################################################################################################################################################
# True survival function
###########################################################################################################################################################

#####  Arguments  #####
# t: time (float >=0)
# lambda and rho: parameters of the Weibull distribution (float > 0)
# Z: arm (0 for control, 1 for intervention)
# beta: treatment effect (float)
# gamma: parameter of the gamma distribution of the frailty terms (float >0)

####  Values ####
# survival at time t


survival_function <- function(t, lambda, rho, Z, beta, gamma){
  if(gamma == 0 ) {
    exp(-lambda*(t^rho*exp(beta*Z)))
  } else {
    1/(1+(lambda*t^(rho)*exp(Z*beta))/gamma)^gamma
  }
}
