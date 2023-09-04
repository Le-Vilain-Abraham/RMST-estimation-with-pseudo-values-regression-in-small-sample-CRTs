##########################################################################################################################################################
# Coverage rate and ype I eo ae/Powe
##########################################################################################################################################################

#####  Arguments  #####
# t_star: horizon time (integer > 0)
# lambda and rho: parameters of the Weibull distribution (float > 0)
# beta: treatment effect (float)
# data: estimations on the D dataset of the  (data.frame)
# gamma: parameter of the gamma distribution of the frailty terms (float >0)

####  Values ####
# Coverage rate and ype I eo ae/Powe (data.frame)


performance_measures <- function(t.star, lambda, rho, beta, data, gamma){
  D <- nrow(data)
  
  #True difference in RMSTs
  if (beta !=0){
    true.delta <- true_rmst_difference(t.star, 0.0001,lambda, rho, beta, gamma)
  } else {
    true.delta <- 0
  }

  coverage <- mean(data$ci_low <= true.delta & data$ci_up >= true.delta) * 100 #True delta is in [ci.low,ci.up]
  
  rejection.rate <- mean(data$ci_low > 0 | data$ci_up < 0) * 100 #0 is not is [ci.low,ci.up]
  
  #Results
  return(data.frame("rejection.rate" = rejection.rate, 
                    "coverage" = coverage, 
                    "D" = D))
}
