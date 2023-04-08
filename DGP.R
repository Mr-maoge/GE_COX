library(MASS)
library(Matrix)

dgp1 <- function(n_=100,p=50,q=4,
                 cor1_type=1,cor2_type=1,rho1=0.3,rho2=0.3,censoring=0.25,
                 beta=NULL,alpha=NULL,gamma=NULL,signal=0.4,
                 n_main=10,n_interactions=16,
                 c_rate=NULL,
                 n_val=0,n_tes=0,seed=0)
{
  # n_: sample size ; 
  # p: num of GE; 
  # q: num of environment features
  # cor_type: 1(AR type); 2(CS: constant covariance)
  # rho1, rho2, rho3: covariance param fo X
  # censoring: rate of censoring
  # beta: p-dimensional vector; alpha:q-dimensional; gamma: q*p matrix
  # c_rate: portion of censoring
  n = n_+n_val+n_tes
  #### initialization steps ####
  set.seed(10)  # seed for random coefficients
  p0 = 200
  # initialize coefficients
  if(is.null(beta)){
    if(T){
    idxs = sample(1:p0,n_main,replace=F)
    #beta = sparseVector(c(rep(0.3,5),rep(-0.3,5)),1:10,r)
    beta = sparseVector(runif(n_main,signal/2,signal),
                        idxs,p)
    idxs2 = sample(0:(n_main*q-1),n_interactions,replace=F)
    j = idxs[idxs2%/%q + 1]
    i = idxs2%%q+1
    gamma = sparseMatrix(i=i,j=j,
                         x=runif(n_interactions,signal/2,signal),
                         dims=c(q,p))  
    }else{
      # main-effect interaction (group structure)
      idxs = sample(1:p0,n_main,replace=F)
      beta = sparseVector(runif(n_main,signal/2,signal),
                          idxs,p)
      j = rep(idxs,q)
      i = rep(1:q,rep(n_main,q))
      gamma = sparseMatrix(i=i,j=j,
                           x=runif(n_main*q,signal/2,signal),
                           dims=c(q,p))  
    }
    
  }else{
    beta = as(beta, "sparseVector")
    gamma = as(gamma, "sparseMatrix")
  }
  
  if(is.null(alpha)){
    #alpha = rep(0.2,q)
    alpha = runif(q,signal/2,signal)
  }

  phi1 = c(beta,as.vector(t(gamma)))
  set.seed(seed) # seed for generate data
  if(cor1_type==1){
    # AR structure 
    Sigma1 = rho1^(abs(outer(1:p,1:p,"-")))
  }else if(cor1_type==2){
    # constant covariance
    Sigma1 = matrix(rho1,p,p) + diag(1-rho1,p,p)
  }else if(cor1_type==3){
    # banded 1
    Sigma1 = diag(1,p,p)+0.4*(abs(outer(1:p,1:p,"-"))==1)
  }else if(cor1_type==4){
    # banded 2
    Sigma1 = diag(1,p,p)+0.6*(abs(outer(1:p,1:p,"-"))==1)+0.2*(abs(outer(1:p,1:p,"-"))==2)
  }
  
  if(cor2_type==1){
    # AR structure 
    Sigma2 = rho2^(abs(outer(1:q,1:q,"-")))
  }else if(cor2_type==2){
    # constant covariance
    Sigma2 = matrix(rho2,q,q) + diag(1-rho2,q,q)
  }else if(cor2_type==3){
    # banded 1
    Sigma2 = diag(1,q,q)+0.4*(abs(outer(1:q,1:q,"-"))==1)
  }else if(cor2_type==4){
    # banded 2
    Sigma2 = diag(1,q,q)+0.6*(abs(outer(1:q,1:q,"-"))==1)+0.2*(abs(outer(1:q,1:q,"-"))==2)
  }
  
  #### generate data ####
  X = mvrnorm(n,mu=rep(0,p),Sigma = Sigma1)
  E = mvrnorm(n,mu=rep(0,q),Sigma = Sigma2)
  
  E[,((q+1)%/%2+1):q] = apply(E[,((q+1)%/%2+1):q],c(1,2),function(x)as.numeric(x>-0.7))
  #if(E_is_continuous!=1){
  #  E = apply(E,c(1,2),function(x)as.numeric(x>-0.7))
  #}
  interaction_term = rep(NA,n) 
  for(i in 1:n){
    interaction_term[i] = sum(gamma*(outer(E[i,],X[i,],"*")))
  }
  # lazard rate
  lambda = exp(1+(X%*%beta)[,1]+(E%*%alpha)[,1]+interaction_term)
  # survival time 
  T_ = rexp(n,rate=lambda)
  
  # generate censoring (use exponential distribution to generate censoring time C)
  # the parameter of this exponential distribution is set to met the censoring rate requirement
  if(is.null(c_rate)){
    f <- function(rate) (1-mean(exp(-rate*T_)))-censoring
    c_rate = uniroot(f,interval=c(1e-10,1e10))$root
  }
  # censoring time 
  C = rexp(n,c_rate)
  # observed survival time
  Y = pmin(T_,C)
  # censoring indicator
  N_ = as.integer(T_<=C)  # = 0 if censoring
  
  phi1_all = c(phi1,alpha)
  res = list(X=X[1:n_,],E=E[1:n_,],Y=Y[1:n_],T_=T_[1:n_],N_=N_[1:n_],
             phi1_all=phi1_all,
             beta=beta,gamma=gamma,c_rate=c_rate)
  if(n_val>0){
    res$X_val = X[(n_+1):(n_+n_val),]
    res$E_val = E[(n_+1):(n_+n_val),]
    res$Y_val = Y[(n_+1):(n_+n_val)]
    res$T_val = T_[(n_+1):(n_+n_val)]
    res$N_val = N_[(n_+1):(n_+n_val)]
  }
  
  if(n_tes>0){
    res$X_tes = X[(n_+n_val+1):(n_+n_val+n_tes),]
    res$E_tes = E[(n_+n_val+1):(n_+n_val+n_tes),]
    res$Y_tes = Y[(n_+n_val+1):(n_+n_val+n_tes)]
    res$T_tes = T_[(n_+n_val+1):(n_+n_val+n_tes)]
    res$N_tes = N_[(n_+n_val+1):(n_+n_val+n_tes)]
  }
  return(res)
}



