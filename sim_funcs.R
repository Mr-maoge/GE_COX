source("DGP.R")
source("Method.R")
source("utils.R")
library(abind)

extract_res<- function(result)
{
  return(c(result$FNR,result$FPR,result$FR,
           result$TP,result$FP,
           result$TP_main,result$FP_main,result$TP_inter,result$FP_inter,
           result$MSE1,result$MSE2,result$pl_tes,result$C_index_tes))
}

sim.combine <- function(a,b)
{
  return(list(result=abind(a$result,b$result,along=1),
              tunings=rbind(a$tunings,b$tunings)))
}

sim.example1 <- function(seed,n=100,p=50,n_main=12,n_interactions=16,E_is_continuous=1,
                         cor1_type=1,cor2_type=1,rho1=0.3,rho2=0.3,
                         signal=0.6,metric="pl_val",censoring=0.4,
                         lambda1_set=exp(seq(log(0.005),log(0.2),length.out = 10)),
                         lambda2_set=lambda1_set,
                         n.lambda2 = 3,
                         mode="normal"){
  source("DGP.R")
  source("Method.R")
  tryCatch({
  if(mode=="trend"){
    # to see the trend of metrics when n,p,r changes 
    # fix the coefficients 
    # main effects (G) 
    set.seed(42) 
    q = 5
    idxs = sample(1:200,n_main,replace=F)
    vals = runif(n_main,signal/2,signal)
    beta = sparseVector(vals,idxs,p)
    
    # interactions
    idxs2 = sample(0:(n_main*q-1),n_interactions,replace=F)
    j = idxs[idxs2%/%q + 1]
    i = idxs2%%q+1
    vals = runif(n_interactions,signal/2,signal)
    gamma = sparseMatrix(i=i,j=j,
                         x=vals,
                         dims=c(q,p))  
    
    # main effects (E) 
    alpha = runif(q,signal/2,signal)
  }else{
    beta = NULL
    alpha = NULL
    gamma = NULL
  }
  
  #### generate data ####
  data = dgp1(seed=seed,n_=n,p=p,q=5,n_val =50,n_tes=400,
              cor1_type = cor1_type,cor2_type = cor2_type,rho1=rho1,rho2=rho2,
              beta=beta,alpha=alpha,gamma=gamma,
              n_main = n_main,n_interactions = n_interactions,
              signal=signal,censoring=censoring)
  phi1_all = data$phi1
  X = data$X
  E = data$E
  Y = data$Y
  N_ = data$N_
  p = dim(X)[2]
  q = dim(E)[2]
  
  if(metric=="pl_val"){
    X_val = data$X_val
    E_val = data$E_val
    Y_val = data$Y_val
    N_val = data$N_val
  }else{
    X_val = NULL
    E_val = NULL
    Y_val = NULL
    N_val = NULL
  }
  
  X_tes = data$X_tes
  E_tes = data$E_tes
  Y_tes = data$Y_tes
  N_tes = data$N_tes

  #### simulation ####
  alt2.lasso = list()
  alt2.MCP = list()
  alt2.grMCP = list()
  alt3 = list()
  our = list() 
  
  # Lasso to generate initial value
  alt2.lasso = ALT_2(E,X,Y,N_,
                  penalty="lasso",
                  phi1_all_true=phi1_all,metric = "bic",
                  E_val = E_val,Z_val=X_val,Y_val=Y_val,N_val=N_val,
                  E_tes = E_tes,Z_tes=X_tes,Y_tes=Y_tes,N_tes=N_tes,best=T)
  # MCP
  alt2.MCP = ALT_1.path(E,X,Y,N_,lambda1_set = 1.5*lambda1_set,
                   lambda2_set=1,n.lambda2 = 3,
                   penalty="MCP",metric=metric,early_stop = F,
                   phi1_all=alt2.lasso$phi1_all,phi1_all_true=phi1_all,
                   E_val = E_val,Z_val=X_val,Y_val=Y_val,N_val=N_val,
                   E_tes = E_tes,Z_tes=X_tes,Y_tes=Y_tes,N_tes=N_tes)$best_model
  
  # grMCP
  alt2.grMCP = ALT_1.path(E,X,Y,N_,lambda1_set = 0.5*lambda1_set,#1.2*lambda1_set,
                        lambda2_set=1,n.lambda2 = 3,
                        penalty="grMCP",metric=metric,early_stop = F,
                        phi1_all=alt2.lasso$phi1_all,phi1_all_true=phi1_all,
                        E_val = E_val,Z_val=X_val,Y_val=Y_val,N_val=N_val,
                        E_tes = E_tes,Z_tes=X_tes,Y_tes=Y_tes,N_tes=N_tes)$best_model
  
  # Marginal analysis
  thresh_adj = ifelse(cor1_type==2,0.05,0.5)  # cor1_type=2: constant covariance structure, the p-values are smaller
  method_adj = ifelse(cor1_type==2,"bonferroni","fdr")
  alt3 = ALT_3(E,X,Y,N_,phi1_all, method_adj, thresh_adj)
  
  # proposed method
  our = ALT_1.path(E,X,Y,N_,lambda1_set = lambda1_set,
                   lambda2_set=lambda2_set,n.lambda2 = 3,
                    penalty="ourMCP",metric=metric,early_stop = F,
                    phi1_all=alt2.lasso$phi1_all,phi1_all_true=phi1_all,
                    E_val = E_val,Z_val=X_val,Y_val=Y_val,N_val=N_val,
                    E_tes = E_tes,Z_tes=X_tes,Y_tes=Y_tes,N_tes=N_tes)$best_model

  #### output the results ####
  result = array(NA,dim=c(1,4,14))
  result[1,1,] = c(our$auc,extract_res(our))
  result[1,2,] = c(alt2.MCP$auc,extract_res(alt2.MCP))
  result[1,3,] = c(alt2.grMCP$auc,extract_res(alt2.grMCP))
  result[1,4,1:10] = c(alt3$auc,extract_res(alt3)) 

  metrics = c("AUC","FNR","FPR","FNR+FPR",
              "TP","FP","TP_main","FP_main","TP_inter","FP_inter",
              "MSE1","MSE2","pl_test","C_index_tes")
  
  methods = c("Proposed","MCP","grMCP","Marginal")
  
  dimnames(result) = list("replicate",methods,metrics)
  
  #### output tuning parameters ####
  tunings = c(our$lambda1,our$lambda2,
              alt2.MCP$lambda1,
              alt2.grMCP$lambda1
              )
  
  names(tunings) = c("lambda1","lambda2","MCP","grMCP")
  
  return(list(result=result,tunings=tunings))
  }
  ,error=function(e){
    print(e$message)
  })
  return(list(result=NULL,tunings=NULL))
}











