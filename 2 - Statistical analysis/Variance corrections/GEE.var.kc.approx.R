##########################################################################################################################################################
# Kauerman and Carroll variance estimator
##########################################################################################################################################################

# Function from the geesvm package modified to used the approximation suggested by Gallis et al. (2020)
# Original code available at https://github.com/cran/geesmv

#####  Arguments  #####
# formula: specify the model of interest
# id: vector which identifies the clusters vector)
# family:  family ("gaussian", "binomial" or "poisson")
# data: dataset to analyse (data.frame)
# corstr: correlation matrix ("independence" for the independent correlation matrix or "exchangeable" for the exchangeable correlation matrix)

####  Values ####
# cov.beta: estimate of robust variance for \hat{\beta} (vector)
# cov.var estimate of the variance-covariance matrix for the variance estimator


GEE.var.kc.approx <-
  function(formula, id, family = gaussian, data, corstr="independence"){

    # Delete the records with missing data in predictors or outcomes;
    if (is.null(data$id)){
      index <- which(names(data)==id)
      data$id <- data[,index]}
    
    ### na.action: only na.omit is used for gee;
    init <- model.frame(formula, data)
    init$num <- 1:length(init[,1])
    if(any(is.na(init))){
      index <- na.omit(init)$num
      data <- data[index,]
      ### Get the design matrix;
      m <- model.frame(formula, data)
      mt <- attr(m, "terms") 
      data$response <- model.response(m, "numeric")
      mat <- as.data.frame(model.matrix(formula, m))
    }else{
      ### Get the design matrix;
      m <- model.frame(formula, data)
      mt <- attr(m, "terms") 
      data$response <- model.response(m, "numeric")
      mat <- as.data.frame(model.matrix(formula, m))
    }
    
    ### Fit the GEE model to get the estimate of parameters \hat{\beta};
    gee.fit <- gee(formula,data=data,id=id,family=family,corstr=corstr)
    beta_est <- gee.fit$coefficient
    alpha <- gee.fit$working.correlation[1,2]
    len <- length(beta_est)
    len_vec <- len^2
    
    ### Estimate the robust variance for \hat{\beta}
    data$id <- gee.fit$id
    cluster<-cluster.size(data$id)
    ncluster<-max(cluster$n)
    size<-cluster$m
    mat$subj <- rep(unique(data$id), cluster$n)
    if(is.character(corstr)){
      var <- switch(corstr,
                    "independence"=cormax.ind(ncluster),
                    "exchangeable"=cormax.exch(ncluster, alpha),
                    "AR-M"=cormax.ar1(ncluster, alpha),
                    "unstructured"=summary(gee.fit)$working.correlation,)
    }else{
      print(corstr)
      stop("'working correlation structure' not recognized")
    }   
    if(is.character(family)){
      family <- switch(family,
                       "gaussian"="gaussian",
                       "binomial"="binomial",
                       "poisson"="poisson")
    }else{ 
      if(is.function(family)){
        family <- family()[[1]]
      }else{
        print(family)
        stop("'family' not recognized")
      }    
    }
    
    cov.beta<-unstr<-matrix(0,nrow=len,ncol=len)
    step11<-matrix(0, nrow=len, ncol=len)
    for (i in 1:size){
      y<-as.matrix(data$response[data$id==unique(data$id)[i]])
      covariate<-as.matrix(subset(mat[,-length(mat[1,])], mat$subj==unique(data$id)[i]))
      var_i=var[1:cluster$n[i],1:cluster$n[i]]
      if (family=="gaussian"){ 
        xx<-t(covariate)%*%solve(var_i)%*%covariate
        step11<-step11+xx  
      }else if (family=="poisson") {
        D<-mat.prod(covariate, exp(covariate%*%beta_est))
        Vi <- diag(sqrt(c(exp(covariate%*%beta_est))),cluster$n[i])%*%var_i%*%diag(sqrt(c(exp(covariate%*%beta_est))),cluster$n[i])
        xx<-t(D)%*%solve(Vi)%*%D
        step11<-step11+xx
      }else if (family=="binomial"){
        D<-mat.prod(covariate, exp(covariate%*%beta_est)/((1+exp(covariate%*%beta_est))^2))
        Vi <- diag(sqrt(c(exp(covariate%*%beta_est)/(1+exp(covariate%*%beta_est))^2)),cluster$n[i])%*%var_i%*%diag(sqrt(c(exp(covariate%*%beta_est)/(1+exp(covariate%*%beta_est))^2)),cluster$n[i])
        xx<-t(D)%*%solve(Vi)%*%D
        step11<-step11+xx 
      }
    }
    
    
    step12<-matrix(0,nrow=len,ncol=len)
    sumxy<-matrix(0,nrow=len,ncol=len)
    sumyx<-matrix(0,nrow=len,ncol=len)
   
    #step13<-matrix(0,nrow=len_vec,ncol=1)
    #step14<-matrix(0,nrow=len_vec,ncol=len_vec)
    #p<-matrix(0,nrow=len_vec,ncol=size)
   
     for (i in 1:size){
      y<-as.matrix(data$response[data$id==unique(data$id)[i]])
      covariate<-as.matrix(subset(mat[,-length(mat[1,])], mat$subj==unique(data$id)[i]))
      var_i=var[1:cluster$n[i],1:cluster$n[i]]
      if (family=="gaussian"){ 
        
        #JT 25 june 2020: replacing true KC with approximation suggested by gallis 2020 and ford 2018
        
        #xy<-t(covariate)%*%solve(var_i)%*%mat.sqrt.inv(cormax.ind(cluster$n[i])-covariate%*%solve(step11)%*%t(covariate)%*%solve(var_i))%*%(y-covariate%*%beta_est)
        #step12<-step12+xy%*%t(xy)
        #step13<-step13+vec(xy%*%t(xy))
        #p[,i]<-vec(xy%*%t(xy))
        
        xy <- t(covariate)%*%
          solve(var_i)%*%
          solve(cormax.ind(cluster$n[i])-
                  covariate%*%
                  solve(step11)%*%
                  t(covariate)%*%
                  solve(var_i))%*%
          (y-covariate%*%beta_est)%*%
          t(y-covariate%*%beta_est)%*%
          solve(var_i)%*%
          covariate
        
        yx <- t(covariate)%*%
          solve(var_i)%*%
          (y-covariate%*%beta_est)%*%
          t(y-covariate%*%beta_est)%*%
          solve(cormax.ind(cluster$n[i])-
                  t(covariate%*%
                      solve(step11)%*%
                      t(covariate)%*%
                      solve(var_i)))%*%
          solve(var_i)%*%
          covariate
        
        sumxy <- sumxy + xy
        
        sumyx <- sumyx + yx
        
      }else if (family=="poisson") {
        D<-mat.prod(covariate, exp(covariate%*%beta_est))
        Vi <- diag(sqrt(c(exp(covariate%*%beta_est))),cluster$n[i])%*%var_i%*%diag(sqrt(c(exp(covariate%*%beta_est))),cluster$n[i])
        
        #xy<-t(D)%*%solve(Vi)%*%mat.sqrt.inv(cormax.ind(cluster$n[i])-D%*%solve(step11)%*%t(D)%*%solve(Vi))%*%(y-exp(covariate%*%beta_est))
        
        #step12<-step12+xy%*%t(xy)
        #step13<-step13+vec(xy%*%t(xy))
        #p[,i]<-vec(xy%*%t(xy))
        
        xy <- t(D)%*%
          solve(Vi)%*%
          solve(cormax.ind(cluster$n[i])-
                  D%*%
                  solve(step11)%*%
                  t(D)%*%
                  solve(Vi))%*%
          (y-exp(covariate%*%beta_est))%*%
          t(y-exp(covariate%*%beta_est))%*%
          solve(Vi)%*%
          D
        
        yx <- t(D)%*%
          solve(Vi)%*%
          (y-exp(covariate%*%beta_est))%*%
          t(y-exp(covariate%*%beta_est))%*%
          solve(cormax.ind(cluster$n[i])-
                  t(D%*%
                      solve(step11)%*%
                      t(D)%*%
                      solve(Vi)))%*%
          solve(Vi)%*%
          D
        
        
        sumxy <- sumxy + xy
        
        sumyx <- sumyx + yx
        
      }else if (family=="binomial"){
        D<-mat.prod(covariate, exp(covariate%*%beta_est)/((1+exp(covariate%*%beta_est))^2))
        Vi <- diag(sqrt(c(exp(covariate%*%beta_est)/(1+exp(covariate%*%beta_est))^2)),cluster$n[i])%*%var_i%*%diag(sqrt(c(exp(covariate%*%beta_est)/(1+exp(covariate%*%beta_est))^2)),cluster$n[i])
        
        #xy<-t(D)%*%solve(Vi)%*%mat.sqrt.inv(cormax.ind(cluster$n[i])-D%*%solve(step11)%*%t(D)%*%solve(Vi))%*%(y-exp(covariate%*%beta_est)/(1+exp(covariate%*%beta_est)))
        #step12<-step12+xy%*%t(xy)
        #step13<-step13+vec(xy%*%t(xy))
        #p[,i]<-vec(xy%*%t(xy)) 
        
        xy <- t(D)%*%
          solve(Vi)%*%
          solve(cormax.ind(cluster$n[i])-
                  D%*%
                  solve(step11)%*%
                  t(D)%*%
                  solve(Vi))%*%
          (y-exp(covariate%*%beta_est)/(1+exp(covariate%*%beta_est)))%*%
          t(y-exp(covariate%*%beta_est)/(1+exp(covariate%*%beta_est)))%*%
          solve(Vi)%*%
          D
        
        yx <- t(D)%*%
          solve(Vi)%*%
          (y-exp(covariate%*%beta_est)/(1+exp(covariate%*%beta_est)))%*%
          t(y-exp(covariate%*%beta_est)/(1+exp(covariate%*%beta_est)))%*%
          solve(cormax.ind(cluster$n[i])-
                  t(D%*%
                      solve(step11)%*%
                      t(D)%*%
                      solve(Vi)))%*%
          solve(Vi)%*%
          D
        
        sumxy <- sumxy + xy
        
        sumyx <- sumyx + yx
      }    
     }
    
    step12 <- (sumxy + sumyx) / 2
    #for (i in 1:size){
    #  dif<-(p[,i]-step13/size)%*%t(p[,i]-step13/size)
    #  step14<-step14+dif
    #}
    
    cov.beta<-solve(step11)%*%(step12)%*%solve(step11)
    #cov.var<-size/(size-1)*kronecker(solve(step11), solve(step11))%*%step14%*%kronecker(solve(step11), solve(step11))
    
    return(list(cov.beta=diag(cov.beta)))
  }
