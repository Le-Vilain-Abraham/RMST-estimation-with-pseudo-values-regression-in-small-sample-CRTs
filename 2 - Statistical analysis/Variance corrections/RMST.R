##########################################################################################################################################################
# Analysis with the pseudo-value regression
##########################################################################################################################################################

#####  Arguments  #####
# dataset: dataset to analyse (data.frame)
# t_star: horizon time (float > 0)

####  Values ####
#results: estimations of the RMST difference, variance and 95% CI for pseudo-value regression (data.frame) with :
#         - either the standard variance, Mancl and DeRouen, Kauermann and Carroll, Fay and Graubard or Morel and al. estimators 
#         - either the normal or the Student distribution 
#         - either the independent or exchangeable correlation matrix 


RMST <- function(dataset, t_star) {
  
  ############# Computation of pseudo-values for each individual #############
  dataset$id <- dataset$cluster
  data_pseudo <-cbind(dataset, 
                      pv = pseudomean(dataset$time,
                                      dataset$status,
                                      tmax = t_star))
  
  ############################### Independent matrix ############################### 
  ############# Standard estimator #############
  log  <- capture.output(suppressMessages(fit_standard_ind <- gee(pv ~ arm, 
                                                                  data = data_pseudo, 
                                                                  id = cluster, 
                                                                  family = gaussian, 
                                                                  corstr = "independence")))

  if(fit_standard_ind$error == 0){
    ############# MD estimator #############
    log <- capture.output(suppressMessages(fit_md_ind <- GEE.var.md(pv ~ arm, 
                                                                data = data_pseudo, 
                                                                id = id, 
                                                                family = gaussian, 
                                                                corstr = "independence")))
    
    ############# KC estimator #############
    log <- capture.output(suppressMessages(fit_kc_ind <- GEE.var.kc.approx(pv ~ arm, 
                                                                       data = data_pseudo, 
                                                                       id = id, 
                                                                       family = gaussian, 
                                                                       corstr = "independence")))
    
    ############# FG estimator #############
    log <- capture.output(suppressMessages(fit_fg_ind <- GEE.var.fg(pv ~ arm, 
                                                                data = data_pseudo, 
                                                                id = id, 
                                                                family = gaussian, 
                                                                corstr = "independence")))
    
    ############# MBN estimator #############
    log <- capture.output(suppressMessages(fit_mbn_ind <- GEE.var.mbn(pv ~ arm, 
                                                                  data = data_pseudo, 
                                                                  id = id, 
                                                                  family = gaussian, 
                                                                  corstr = "independence")))
    
    
    ############# Results #############
    results_ind <- data.frame(matrix = "ind",
                              method_correction = rep(c("standard","md", "kc", "fg", "mbn"), each = 2),
                              distribution = rep(c("z","t"), 5),
                              delta.rmst = summary(fit_standard_ind)$coefficients["arm","Estimate"],
                              var = rep(c(summary(fit_standard_ind)$coefficients["arm","Robust S.E."]^2,
					  fit_md_ind$cov.beta["arm"],
                                          fit_kc_ind$cov.beta["arm"],
                                          fit_fg_ind$cov.beta["arm"],
                                          fit_mbn_ind$cov.beta["arm"]), each = 2),
                              ci.low = NA,
                              ci.up = NA,
                              t_star = 365)
    
    #z-distribution
    results_ind[results_ind$distribution=="z",]$ci.low <- results_ind[results_ind$distribution=="z","delta.rmst"]-qnorm(0.975)*sqrt(results_ind[results_ind$distribution=="z","var"])
    results_ind[results_ind$distribution=="z",]$ci.up <- results_ind[results_ind$distribution=="z","delta.rmst"]+qnorm(0.975)*sqrt(results_ind[results_ind$distribution=="z","var"])
    
    #t-distribution
    results_ind[results_ind$distribution=="t",]$ci.low <- results_ind[results_ind$distribution=="t","delta.rmst"]-qt(0.975, df=length(unique(dataset$cluster))-2)*sqrt(results_ind[results_ind$distribution=="t","var"])
    results_ind[results_ind$distribution=="t",]$ci.up <- results_ind[results_ind$distribution=="t","delta.rmst"]+qt(0.975, df=length(unique(dataset$cluster))-2)*sqrt(results_ind[results_ind$distribution=="t","var"])
    
  }else{
    results_ind <- data.frame(matrix = "ind",
                              method_correction = rep(c("standard","md", "kc", "fg", "mbn"), each = 2),
                              distribution = rep(c("z","t"), 5),
                              delta.rmst = fit_standard_ind$error,
                              var = fit_standard_ind$error,
                              ci.low = fit_standard_ind$error,
                              ci.up = fit_standard_ind$error,
                              t_star = 365)
  }

  
  
  ############################### Exchangeable matrix ############################### 
  ############# Standard estimator #############
  log  <- capture.output(suppressMessages(fit_standard_exc <- gee(pv ~ arm, 
                                                                  data = data_pseudo, 
                                                                  id = cluster, 
                                                                  family = gaussian, 
                                                                  corstr = "exchangeable")))
  if(fit_standard_exc$error == 0){
    ############# MD estimator #############
    log <- capture.output(suppressMessages(fit_md_exc <- GEE.var.md(pv ~ arm, 
                                                                data = data_pseudo, 
                                                                id = id, 
                                                                family = gaussian, 
                                                                corstr = "exchangeable")))
    
    ############# KC estimator #############
    log <- capture.output(suppressMessages(fit_kc_exc <- GEE.var.kc.approx(pv ~ arm, 
                                                                       data = data_pseudo, 
                                                                       id = id, 
                                                                       family = gaussian, 
                                                                       corstr = "exchangeable")))
    
    ############# FG estimator #############
    log <- capture.output(suppressMessages(fit_fg_exc <- GEE.var.fg(pv ~ arm, 
                                                                data = data_pseudo, 
                                                                id = id, 
                                                                family = gaussian, 
                                                                corstr = "exchangeable")))
    
    ############# MBN estimator #############
    log <- capture.output(suppressMessages(fit_mbn_exc <- GEE.var.mbn(pv ~ arm, 
                                                                  data = data_pseudo, 
                                                                  id = id, 
                                                                  family = gaussian, 
                                                                  corstr = "exchangeable")))
    
    
    ############# Results #############
    results_exc <- data.frame(matrix = "exc",
                              method_correction = rep(c("standard","md", "kc", "fg", "mbn"), each = 2),
                              distribution = rep(c("z","t"), 5),
                              delta.rmst = summary(fit_standard_exc)$coefficients["arm","Estimate"],
                              var = rep(c(summary(fit_standard_exc)$coefficients["arm","Robust S.E."]^2,
					     fit_md_exc$cov.beta["arm"],
                                            fit_kc_exc$cov.beta["arm"],
                                            fit_fg_exc$cov.beta["arm"],
                                            fit_mbn_exc$cov.beta["arm"]), each=2),
                              ci.low = NA,
                              ci.up = NA,
                              t_star = 365)
    
    #z-distribution
    results_exc[results_exc$distribution=="z",]$ci.low <- results_exc[results_exc$distribution=="z","delta.rmst"]-qnorm(0.975)*sqrt(results_exc[results_exc$distribution=="z","var"])
    results_exc[results_exc$distribution=="z",]$ci.up <- results_exc[results_exc$distribution=="z","delta.rmst"]+qnorm(0.975)*sqrt(results_exc[results_exc$distribution=="z","var"])
    
    #t-distribution
    results_exc[results_exc$distribution=="t",]$ci.low <- results_exc[results_exc$distribution=="t","delta.rmst"]-qt(0.975, df=length(unique(dataset$cluster))-2)*sqrt(results_exc[results_exc$distribution=="t","var"])
    results_exc[results_exc$distribution=="t",]$ci.up <- results_exc[results_exc$distribution=="t","delta.rmst"]+qt(0.975, df=length(unique(dataset$cluster))-2)*sqrt(results_exc[results_exc$distribution=="t","var"])
    
  }else{
    results_exc <- data.frame(matrix  = "exc",
                              method_correction = rep(c("standard","md", "kc", "fg", "mbn"), each = 2),
                              distribution = rep(c("z","t"),5),
                              delta.rmst = fit_standard_exc$error,
                              var = fit_standard_exc$error,
                              ci.low = fit_standard_exc$error,
                              ci.up = fit_standard_exc$error,
                              t_star = 365)
  }
  
  return(rbind(results_ind,results_exc))
}
