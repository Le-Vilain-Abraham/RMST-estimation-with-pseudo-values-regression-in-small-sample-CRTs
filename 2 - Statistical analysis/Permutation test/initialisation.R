########################################################################################
# Initialisation of the search procedure
########################################################################################

#####  Arguments  #####
# data: dataset to analyse (data.frame)
# matrix: working correlation matrix (character: "independence" for independent matrix)
# theta_hat: point estimate of beta1 (difference in RMST) based on original dataset
# alpha: significance level of the test (float between 0 and 1)

####  Values ####
#results: initialisation for the search procedure (data.frame)


initialisation <- function(data, matrix, theta_hat, alpha = 0.05) {
  
  #######
  #init m : p pour le step 1
  m <- min(c(50, ceiling(0.3 * (4 - alpha) / alpha)))
  
  #######
  #calcul de k utilisé pour le calcul de c le "step lenght constant"
  z <- qnorm(1 - (alpha / 2))
  k <- 2 * sqrt(2 * pi) * exp(z^2 / 2) / z
  
  #######
  #init lower and upper bounds 
  nperm_init <- ceiling((4 - alpha) / alpha) # nombre de pemrutation pour le test
  data$theta_hat <- theta_hat
  permuted_theta <- replicate(nperm_init, initialisation_ci(data, matrix))
  
  t1 <- sort(permuted_theta)[2] # 2nd to smallest
  t2 <- sort(permuted_theta)[nperm_init - 1] # 2nd to largest
  low <- theta_hat - ((t2 - t1) / 2)
  up  <- theta_hat + ((t2 - t1) / 2)

  return(data.frame(m = m, k = k, low = low, up = up))
}
