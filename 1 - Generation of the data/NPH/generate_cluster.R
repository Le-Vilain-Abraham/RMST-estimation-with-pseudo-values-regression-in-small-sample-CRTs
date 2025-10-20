##########################################################################################################################################################
# Generate time-to-event data for one cluster
##########################################################################################################################################################

#####  Arguments  #####
# m: cluster sizes (int > 0)
# lambda and rho: parameters of the Weibull distribution (floats > 0)
# u_k: frailty term of the cluster k (float > 0)
# beta : intervention effect (float)
#t_delay : change point
# Z: arm (1 = intervention arm or 0 = control arm)

####  Values ####
# t_k: simulated time-to-event data for one cluster (vector) 


generate_cluster<- function(m, lambda, rho, u_k, beta, t_delay, Z) {
  
  U <- runif(n = m, min = 0, max = 1) #Uniform distribution
  
  t_k <- ifelse(U > exp(-u_k*lambda*t_delay^rho),
                (-log(U)/(u_k*lambda))^(1/rho),
                ((((-log(U)/(lambda*u_k))-t_delay^rho)/exp(Z*beta)) + t_delay^rho) ^ (1/rho))
  
  return(t_k)
}

