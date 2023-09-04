##########################################################################################################################################################
# Measures of performance
##########################################################################################################################################################

#####  Arguments  #####
# t.star: horizon time (integer > 0)
# lambda and rho: parameters of the Weibull distribution (float > 0)
# beta: treatment effect (float)
# data: estimations on the D dataset of the  (data.frame)
# gamma: parameter of the gamma distribution of the frailty terms (float >0)

####  Values ####
# relative bias (RE), relative error (RE), the coverage rate (coverage), type I error rate (for scenario where HR=1) (rejection.rate) and the number 
#of simulation iterations which converged (D) (data.frame)


performance_measures <- function(t.star, lambda, rho, beta, data, gamma){
  D <- nrow(data)
  
  #True difference in RMSTs
  if (beta !=0){
    true.delta <- true_rmst_difference(t.star, 0.0001,lambda, rho, beta, gamma)
  } else {
    true.delta <- 0
  }
  
  #Mean difference in RMSTs
  mean.delta <- mean(data$delta.rmst)
  
  
  #Measures fo performance
  bias <-mean.delta - true.delta
  relative.bias <- bias/true.delta *100
  
  ASD <- sqrt(mean(data$var))
  ESD <- sqrt(sum((data$delta.rmst - mean.delta)^2)/(D - 1)) 

  relative.error <- (ASD-ESD)/ESD *100
    
  coverage <- mean(data$ci.low <= true.delta & data$ci.up >= true.delta) * 100 #True delta is in [ci.low,ci.up]
  
  rejection.rate <- mean(data$ci.low > 0 | data$ci.up < 0) * 100 #0 is not is [ci.low,ci.up]
  
  #Results
  return(data.frame("relative.error" = relative.error, 
                    "coverage" = coverage, 
                    "rejection.rate" = rejection.rate, 
                    "D" = D))
}
