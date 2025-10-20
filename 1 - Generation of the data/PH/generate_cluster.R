##########################################################################################################################################################
# Generate time-to-event data for one cluster
##########################################################################################################################################################

#####  Arguments  #####
# m: cluster sizes (int > 0)
# lambda and rho: parameters of the Weibull distribution (floats > 0)
# u_k: frailty term of the cluster k (float > 0)
# beta : intervention effect (float)
# Z: arm (1 = intervention arm or 0 = control arm)

####  Values ####
# t_k: simulated time-to-event data for one cluster (vector) 


generate_cluster <- function(m, lambda, rho, u_k, beta, Z) {
  
  U <- runif(n = m, min = 0, max = 1) #Uniform distribution
  
  t_k <- (-log(U)/(u_k*lambda*exp(beta*Z)))^(1/rho) #TTE
  
  return(t_k)
}
