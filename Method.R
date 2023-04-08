library(glmnet)
library(ncvreg)
library(grpreg)
library(survival)
library(matrixStats)
source("utils.R")

#### auxiliary functions (used to calculate C-index) ####
C_index0 <- function(risk_ord,N_ord,Y_ord){
  n = length(risk_ord)
  cnt1 = 0 # total pairs
  cnt2 = 0 # concordant pair
  for(i in 1:(n-1)){
    for(j in i:n){
      if((N_ord[j]==1)&(Y_ord[i]>=Y_ord[j])){ # i not censoring
        cnt1 = cnt1+1
        cnt2 = cnt2+((risk_ord[i]<=risk_ord[j])|(Y_ord[i]==Y_ord[j]))
      }
    }
  }
  return(cnt2/cnt1)
}

C_index <- function(phi, N_ord, A_ord, Y_ord){
  n = dim(A_ord)[1]
  risk = as.vector(A_ord%*%phi)
  cnt1 = 0 # total pairs
  cnt2 = 0 # concordant pair
  for(i in 1:(n-1)){
    for(j in i:n){
      if((N_ord[j]==1)&(Y_ord[i]>=Y_ord[j])){ # i not censoring
        cnt1 = cnt1+1
        cnt2 = cnt2+((risk[i]<=risk[j])|(Y_ord[i]==Y_ord[j]))
      }
    }
  }
  return(cnt2/cnt1)
}
#### auxiliary functions (for Cox model part) ####
# partial_loglikelihood 
partial_loglike <- function(phi,N_ord,A_ord,Y_ord_raw,part=NULL)
{
  n = length(Y_ord_raw)
  tmp = (A_ord%*%phi)[,1]
  tmp2 = cumsum(exp(tmp))
  for(i in (n-1):1){
    if(Y_ord_raw[i]==Y_ord_raw[i+1]){
      tmp2[i] = tmp2[i+1]
    }
  }
  if(is.null(part)){
    pl = mean(N_ord*(tmp-log(tmp2)))
  }else{
    pl = mean((N_ord*(tmp-log(tmp2)))[part]) 
  }
  return(pl)
}

# function used to calculate the partial-likelihood for validation set
partial_loglike_val <- function(phi,Y_train,N_train,A_train,Y_val,N_val,A_val)
{
  n_train = length(Y_train)
  n_val = length(Y_val)
  Y_all = c(Y_train,Y_val)
  N_all = c(N_train,N_val)
  A_all = rbind(A_train,A_val)
  
  Y_ord_all = order(Y_all,decreasing = T)
  N_ord_all = N_all[Y_ord_all]
  A_ord_all = A_all[Y_ord_all,]
  
  part = which(Y_ord_all>n_train)
  pl_val = partial_loglike(phi,N_ord_all,A_ord_all,Y_ord_all,part=part)
  return(pl_val)
}

# first order derivative for partial loglikelihood (gradient)
partial_loglike.first <- function(phi,N_ord,A_ord,Y_ord_raw)
{
  n = length(Y_ord_raw)
  tmp = exp((A_ord%*%phi)[,1])
  tmp2 = cumsum(tmp)
  tmp3 = colCumsums(A_ord*tmp)
  for(i in (n-1):1){
    if(Y_ord_raw[i]==Y_ord_raw[i+1]){
      tmp2[i] = tmp2[i+1]
      tmp3[i,] = tmp3[i+1,]
    }
  }
  pl.first = apply((A_ord-tmp3/tmp2)*N_ord,2,mean)
  return(pl.first)
}



#### auxiliary functions (for penalty part) ####
# soft-thresholding function
St <- function(x,lambda){
  sign(x)*pmax(abs(x)-lambda,0)
}

scad <- function(x,lambda,gamma=3.7)
{
  # proximal map func for SCAD
  x_up = (abs(x)<=2*lambda)*St(x,lambda)+
         (abs(x)>2*lambda)*(abs(x)<=gamma*lambda)*St(x,gamma*lambda/(gamma-1))/(1-1/(gamma-1))+
         (abs(x)>gamma*lambda)*x
  return(x_up)
}

mcp <- function(x,lambda,gamma=3)
{
  # proximal map func for MCP
  x_up = (abs(x)<=gamma*lambda)*St(x,lambda)/(1-1/gamma)+ (abs(x)>gamma*lambda)*x
  return(x_up)
}

# first order derivative for penalty function(non-zero region)
St.d <- function(x,lambda)
{
  return(sign(x)*lambda)
}

scad.d <- function(x,lambda,gamma=3.7)
{
  x0 = abs(x)
  d = (x0<=lambda)*lambda+
    (x0<=gamma*lambda)*(x0>lambda)*(gamma*lambda-x0)/(gamma-1)
  return(sign(x)*d)
}

mcp.d <- function(x,lambda,gamma=3)
{
  x0 = abs(x)
  d = (x0<=lambda*gamma)*(lambda-x/gamma)
  return(sign(x)*d)
}

#### function for update in second step of ADMM (proximal map) ####
# updatate step II
## penalty function for our method (sparse group penalty)
update.II <- function(a,lambda1,lambda2=lambda1,penalty="ourMCP"){
  q = length(a)-1
  st_res = St(a[2:(q+1)],lambda2)
  tmp = c(a[1],st_res)
  ## Case I: all coefficients in group is 0
  if(sqrt(sum(tmp^2))<=lambda1*sqrt(q+1)){
    return(rep(0,q+1))
  }
  ## Case II: some coefficients are non-zero
  # result for sparse-group lasso penalty
  if(penalty=="ourLasso"){
    return((1-lambda1*sqrt(q+1)/sqrt(sum(tmp^2)))*tmp)
  }
  
  # iterative to find coefficients within group (for sparse-gropu scad/mcp)
  a = c(a[1],a[2:(q+1)]*(abs(a[2:(q+1)])>lambda2))
  h = a
  for(iter in 1:100){
    h0 = h
    h_norm = sqrt(sum(h0^2))
    if(penalty=="ourSCAD"){
      d = scad.d(h_norm,sqrt(q+1)*lambda1)
      h[1] = a[1]/(1+d/h_norm)
      h[2:(q+1)] = (a[2:(q+1)]-scad.d(h0[2:(q+1)],lambda2))/(1+d/h_norm)
    }else if(penalty=="ourMCP"){
      d = mcp.d(h_norm,sqrt(q+1)*lambda1)
      h[1] = a[1]/(1+d/h_norm)
      h[2:(q+1)] = (a[2:(q+1)]-mcp.d(h0[2:(q+1)],lambda2))/(1+d/h_norm)
    }
    if(mean((h-h0)^2)<=1e-3){
      break
    }
  }
  return(h)
}

## for grMCP/grSCAD and cMCP
update.II.2 <- function(a,lambda1,penalty="grMCP"){
  q = length(a)-1
  # all group is 0
  if(sqrt(sum(a^2))<=lambda1*sqrt(q+1)){
    return(rep(0,q+1))
  }
  # iterative to find coefficients within group
  h = a
  for(iter in 1:100){
    h0 = h
    h_norm = sqrt(sum(h0^2))
    if(penalty=="grLasso"){
      d = sqrt(q+1)*lambda1
    }else if(penalty=="grSCAD"){
      d = scad.d(h_norm,sqrt(q+1)*lambda1)
    }else if(penalty=="grMCP"){
      d = mcp.d(h_norm,sqrt(q+1)*lambda1)
    }
    h = a/(1+d/h_norm)
    if(mean((h-h0)^2)<=1e-3){
      break
    }
  }
  return(h)
}




#### function for update in first step of ADMM ####
update.ALT_1.I <- function(phi1_all,N_ord,A_all_ord,Y_all_ord,rho,Tau,q,r){
  # update the parameter using quasi-newton method
  n = dim(A_all_ord)[1]
  idx1 = 1:(r+q*r)
  obj <- function(phi1_all){
    B_t = matrix(phi1_all[idx1],r,q+1)
    objective = -partial_loglike(phi1_all,N_ord,A_all_ord,Y_all_ord)+rho/2*sum((B_t-Tau)^2)
    return(objective)
  }
  obj.first <- function(phi1_all){
    B_t = matrix(phi1_all[idx1],r,q+1)
    g = as.vector(rho*(B_t-Tau))
    gr = -partial_loglike.first(phi1_all,N_ord,A_all_ord,Y_all_ord)
    gr[idx1] = gr[idx1]+g
    return(gr)
  }
  return(optim(phi1_all,fn=obj,gr=obj.first,method="BFGS")$par)
}

#### model fitting (proposed/ MCP/ grMCP) ####
# ALT_1: function for: 
#   proposed method (set penalty="ourMCP"); 
#   MCP: (set penalty="MCP")
#   grMCP: (set penalty="grMCP")
ALT_1 <- function(E,Z,Y,N_,lambda1,lambda2=lambda1,penalty="ourMCP",
                  phi1_all=NULL,
                  maxit1=100,
                  epsilon1=5e-4,epsilon2=5e-4,
                  phi1_all_true=NULL,
                  E_val=NULL,Z_val=NULL,Y_val=NULL,N_val=NULL,
                  E_tes=NULL,Z_tes=NULL,Y_tes=NULL,N_tes=NULL
                  ){
  # E: environmantal factors; Z: gene factors;
  # Y: observed survival time; 
  # N_: survival indicator for patients (N_=0 indicates censoring);
  # lambda1, lambda2: tuning parameter for penalty 
  #                   (if penalty only requires 1 tuning parameter, lambda2 is not used)
  # phi1_all: initial value
  
  # maxit1: maxit for ADMM; 
  # epsilon1_1: for primal feasibility condition; 
  # epsilon1_2: for dual feasibility condition;
  
  rho =1
  n = dim(Z)[1]
  q = dim(E)[2]
  r = dim(Z)[2]
  idx1 = 1:(r*(q+1))
  #### initialize #### 
  # initialize parameterss
  if(is.null(phi1_all)){
    phi1_all = runif(r+r*q+q,-0.5,0.5)
  }
  # vectorize variables
  A1 = Z
  for(k in 1:q){
    W = E[,k]*Z
    A1 = cbind(A1,W)
  }
  A1_all = cbind(A1,E)
  # order the data (follow decreasing order of Y)
  Y_ord = order(Y,decreasing = T)
  Y_ord_raw = Y[Y_ord]
  N_ord = N_[Y_ord]
  A1_ord_all = A1_all[Y_ord,]
  
  # Additional variable
  H = matrix(phi1_all[idx1],q+1,r,byrow = T)
  # Lagrangian Multiplier
  Pi1 = matrix(0,q+1,r)
  
  #### model fitting(ADMM) ####
  primal_trace = c()
  dual_trace = c()
  objective_trace = c()
  
  for(i1 in 1:maxit1){
    H0 = H    # these two quantities are used to calculate stoping critiria
    #### Step I ####
    Tau1 = t(H+Pi1)
    phi1_all = update.ALT_1.I(phi1_all,N_ord,A1_ord_all,Y_ord_raw,rho,Tau1,q,r)
    B = matrix(phi1_all[idx1],q+1,r,byrow = T)
    #### Step II ####
    for(j in 1:r){
      if(startsWith(penalty,"our")){
        H[,j] = update.II(B[,j]-Pi1[,j],lambda1/rho,lambda2/rho,penalty=penalty)
      }else if(penalty=="lasso"){
        H[,j] = St(B[,j]-Pi1[,j],lambda1/rho)
      }else if(penalty=="SCAD"){
        H[,j] = scad(B[,j]-Pi1[,j],lambda1/rho)
      }else if(penalty=="MCP"){
        H[,j] = mcp(B[,j]-Pi1[,j],lambda1/rho)
      }else if(penalty%in%c("grLasso","grSCAD","grMCP")){
        H[,j] = update.II.2(B[,j]-Pi1[,j],lambda1/rho,penalty=penalty)
      }
    }
    
    #### step III ####
    Pi1 = H-B+Pi1
    
    #### stoping criteria ####
    stop.primal = mean((H-B)^2)
    stop.dual = mean((H-H0)^2)

    primal_trace = c(primal_trace, stop.primal)
    dual_trace = c(dual_trace, stop.dual)
    
    
    if((stop.primal<=epsilon1)&(stop.dual<=epsilon2)){
      break
    }
  }
  if(i1==maxit1){
    warning(sprintf("ADMM desn't converge in %d iterations.",maxit1))
  }
  
  #### output ####
  phi1_all[idx1] = as.vector(t(H))
  num_nonzeros_select = sum(phi1_all!=0)-q
  #gamma = 1-1/(2*(log(r*q+r)/log(n)))
  gamma = 0.2
  pl = partial_loglike(phi1_all,N_ord,A1_ord_all,Y_ord_raw)
  aic = -2*n*pl + 2*num_nonzeros_select
  bic = -2*n*pl + log(n)*num_nonzeros_select
  ebic = -2*n*pl + log(n)*num_nonzeros_select + 2*gamma*log(choose(r+r*q, num_nonzeros_select))
  result = list(r=r,q=q,
                phi1_all=phi1_all,
                pl=pl,aic=aic,bic=bic,ebic=ebic,
                lambda1=lambda1,lambda2=lambda2,
                primal_trace=primal_trace, dual_trace=dual_trace
  )
  
  if(!is.null(phi1_all_true)){
    phi1_tmp = phi1_all[idx1]
    phi1_true_tmp = phi1_all_true[idx1]
    
    result$TP = sum((phi1_tmp!=0)*(phi1_true_tmp!=0))
    result$FP = sum((phi1_tmp!=0)*(phi1_true_tmp==0))
    # main effect
    phi1_main = phi1_tmp[1:r]
    phi1_true_main = phi1_true_tmp[1:r]
    result$TP_main = sum((phi1_main!=0)*(phi1_true_main!=0))
    result$FP_main = sum((phi1_main!=0)*(phi1_true_main==0))
    # interactions
    phi1_inter = phi1_tmp[(r+1):(r*(q+1))]
    phi1_true_inter = phi1_true_tmp[(r+1):(r*(q+1))]
    result$TP_inter = sum((phi1_inter!=0)*(phi1_true_inter!=0))
    result$FP_inter = sum((phi1_inter!=0)*(phi1_true_inter==0))
    
    result$FNR = sum((phi1_tmp==0)*(phi1_true_tmp!=0))/sum(phi1_true_tmp!=0)
    result$FPR = sum((phi1_tmp!=0)*(phi1_true_tmp==0))/sum(phi1_true_tmp==0)
    result$FR =  result$FNR+result$FPR
    result$MSE1 = sum((phi1_all-phi1_all_true)^2)
    result$MSE2 = mean((A1_ord_all%*%(phi1_all-phi1_all_true))^2)
  }
  if(!is.null(Y_val)){
    A1_val = Z_val
    for(k in 1:q){
      W = E_val[,k]*Z_val
      A1_val = cbind(A1_val,W)
    }
    A1_all_val = cbind(A1_val,E_val)
    pl_val = partial_loglike_val(phi1_all,Y,N_,A1_all,Y_val,N_val,A1_all_val)
    result$pl_val = pl_val
  }
  
  if(!is.null(Y_tes)){
    A1_tes = Z_tes
    for(k in 1:q){
      W = E_tes[,k]*Z_tes
      A1_tes = cbind(A1_tes,W)
    }
    A1_tes = cbind(A1_tes,E_tes)
    Y_tes_ord = order(Y_tes,decreasing = T)
    Y_tes_ord_raw = Y_tes[Y_tes_ord]
    N_tes_ord = N_tes[Y_tes_ord]
    A1_tes_ord = A1_tes[Y_tes_ord,]
    pl_tes = partial_loglike(phi1_all,N_tes_ord,A1_tes_ord,Y_tes_ord_raw)
    result$pl_tes = pl_tes
    C_index_tes = C_index(phi1_all,N_tes_ord,A1_tes_ord,Y_tes_ord_raw)
    result$C_index_tes = C_index_tes
  }
  return(result)
}

# model fitting with parameter selection
ALT_1.path <- function(E,Z,Y,N_,lambda1_set,lambda2_set=NULL,n.lambda2=3,
                       penalty="MCP",metric="bic",early_stop=T,
                       maxit1=300,
                       phi1_all=NULL,phi1_all_true=NULL,
                       E_val=NULL,Z_val=NULL,Y_val=NULL,N_val=NULL,
                       E_tes=NULL,Z_tes=NULL,Y_tes=NULL,N_tes=NULL
                       ){
  # lambda1_set, lambda2_set: candidate tuning parameters
  # metric: method for tuning parameter selection
  
  n.lambda1 = length(lambda1_set)
  best_model = NULL
  metric_val = Inf
  
  q = dim(E)[2]
  r = dim(Z)[2]
  
  freq_pos = rep(0,r+q+r*q)  #  freq of variables being selected
  if(is.null(lambda2_set)){
    d = n.lambda2%/%2
    n.lambda2 = 2*d+1
  }else{
    n.lambda2 = length(lambda2_set)
    d = n.lambda2%/%2
  }
  metrics_result = matrix(NA,n.lambda1,n.lambda2)
  for(k in 1:n.lambda1){
    print(k)
    for(k2 in 1:n.lambda2){
      lambda1 = lambda1_set[k]
      if(is.null(lambda2_set)){
        idx = k-d-1+k2
        if(idx>n.lambda1|idx<1){next}
        lambda2 = lambda1_set[idx]
      }else{
        lambda2 = lambda2_set[k2]
      }
      
      result = ALT_1(E,Z,Y,N_,lambda1,lambda2,penalty=penalty,
                     phi1_all=phi1_all,phi1_all_true = phi1_all_true,
                     E_val=E_val,Z_val=Z_val,Y_val=Y_val,N_val=N_val,
                     E_tes=E_tes,Z_tes=Z_tes,Y_tes=Y_tes,N_tes=N_tes)
      num_coefs = sum(result$phi1_all[1:(r*(q+1))]!=0)
      metric_tmp = result[[metric]]
      if(metric=="pl_val"){
        metric_tmp = -metric_tmp
      }
      
      metrics_result[k,k2] = metric_tmp
      if(is.null(best_model)|((metric_tmp<metric_val)&(num_coefs>0))){
        metric_val = metric_tmp
        best_model = result
      }
      phi1_all = result$phi1_all
      freq_pos = freq_pos+(phi1_all!=0)
      cat("lambda1= ",lambda1,"lambda2=",lambda2," metric= ",metric_tmp," num of variables= ",sum((phi1_all!=0)),"\n")
      if(!is.null(E_tes)){
        cat("C_index_tes = ",result$C_index_tes,"\n")
      }
    }
    if(early_stop&(k>=4)){
      if((metrics_result[k-1,d+1]>metrics_result[k-2,d+1])&(metrics_result[k,d+1]>metrics_result[k-1,d+1])){
        break
      }
    }
  }
  prob_pos = as.vector(freq_pos/(n.lambda1*n.lambda2))
  
  
  result = list(metrics_result=metrics_result,prob_pos=prob_pos)
  if(!is.null(phi1_all_true)){
    # calculate the AUC for variable selection
    true_pos = as.vector(phi1_all_true!=0)
    result$auc = cal.auc(true_pos[1:(r*(q+1))],prob_pos[1:(r*(q+1))])
  }
  #if(!is.null(Y_tes)){
  #  best_model = ALT_1(E,Z,Y,N_,best_model$lambda1,best_model$lambda2,penalty=penalty,
  #                     phi1_all=best_model$phi1_all,phi1_all_true=phi1_all_true,
  #                     E_val=E_val,Z_val=Z_val,Y_val=Y_val,N_val=N_val,
  #                     E_tes=E_tes,Z_tes=Z_tes,Y_tes=Y_tes,N_tes=N_tes)
  #}
  best_model$auc = result$auc 
  result$best_model = best_model
  return(result)
}


#### model fitting (an alternative function) #### 
# ALT_2: use `ncvsurv` package to fit the proposed method 
#        used to generate the initial value in our simulation studies
ALT_2 <- function(E,Z,Y,N_,lambda_set=NULL,penalty="MCP",metric="bic",
                  phi1_all_true=NULL,
                  E_val=NULL,Z_val=NULL,Y_val=NULL,N_val=NULL,
                  E_tes=NULL,Z_tes=NULL,Y_tes=NULL,N_tes=NULL,best=T){
  n = dim(Z)[1]
  q = dim(E)[2]
  r = dim(Z)[2]
  A1 = Z
  for(k in 1:q){
    W = E[,k]*Z
    A1 = cbind(A1,W)
  }
  A1_all = cbind(A1,E)
  y = Surv(Y, N_)
  if(penalty%in%c("Oracle","MCP","SCAD","lasso")){
    if(penalty=="Oracle"){
      penalty.factor = c(rep(0,r),rep(1,q*r),rep(0,q)) 
      penalty = "MCP"
    }else{
      penalty.factor = c(rep(1,(q+1)*r),rep(0,q))
    }
    if(is.null(lambda_set)){
      model = ncvsurv(A1_all,y,penalty = penalty, max.iter = 5000,
                      nlambda = 50,penalty.factor = penalty.factor,
                      alpha=0.99,standardize = FALSE)
    }else{
      lambda_set = sort(lambda_set,decreasing=T)
      model = ncvsurv(A1_all,y,penalty = penalty,lambda=lambda_set,max.iter = 5000,
                      nlambda = 50,penalty.factor = penalty.factor,
                      alpha=0.99,standardize = FALSE)
    }
  }else if(penalty%in%c("grLasso", "grMCP", "grSCAD","gel", "cMCP")){
    group = c(rep(1:r,q+1),rep(0,q))
    if(is.null(lambda_set)){
      model = grpsurv(A1_all,y,group,penalty=penalty,nlambda=50,max.iter = 60000,
                      alpha=0.99,standardize = FALSE)
    }else{
      model = grpsurv(A1_all,y,group,penalty=penalty,lambda=lambda_set,nlambda=50,max.iter = 60000,
                      alpha=0.99,standardize = FALSE)
    }
  }
  n.lambda = length(model$lambda)-1
  #gamma = 1-1/(2*(log(r*q+r)/log(n)))
  gamma = 0.2
  num_nonzeros = apply(model$beta!=0,2,sum)-q
  bic = (2*model$loss+num_nonzeros*log(n) )[1:(n.lambda)]
  ebic =(2*model$loss+num_nonzeros*log(n)+2*gamma*log(choose(r+r*q, num_nonzeros)))[1:(n.lambda)]
  aic = (2*model$loss+num_nonzeros*2 )[1:(n.lambda)]
  freq_pos = rep(0,r+q+r*q)  #  freq of variables being selected
  for(i in 1:n.lambda){
    phi1_all = model$beta[,i]
    freq_pos = freq_pos+(phi1_all!=0)
  }
  prob_pos = as.vector(freq_pos/n.lambda)
  
  if(metric=="bic"){
    idx = which.min(bic)
  }else if(metric=="aic"){
    idx = which.min(aic)
  }else if(metric=="ebic"){
    idx = which.min(ebic)
  }else if(metric=="pl_val"){
    A1_val = Z_val
    for(k in 1:q){
      W = E_val[,k]*Z_val
      A1_val = cbind(A1_val,W)
    }
    A1_all_val = cbind(A1_val,E_val)
    pl_val = rep(NA,n.lambda)
    for(i in 1:n.lambda){
      phi1_all=model$beta[,i]
      pl_val[i] = partial_loglike_val(phi1_all,Y,N_,A1_all,Y_val,N_val,A1_all_val)
    }
    idx = which.max(pl_val)
  }
  
  
  result = list(r=r,q=q,best_idx=idx,aic=aic,bic=bic,ebic=ebic,
                prob_pos=prob_pos)
  
  if(best){
    result$lambda=model$lambda[idx]
    result$phi1_all = model$beta[,idx]
  }else{
    result$lambda=model$lambda
    result$phi1_all = model$beta
  }
  
  if(!is.null(phi1_all_true)){
    true_pos = as.vector(phi1_all_true!=0)
    if(mean(true_pos[1:(r*(q+1))])!=1){
      result$auc = cal.auc(true_pos[1:(r*(q+1))],prob_pos[1:(r*(q+1))])
    }else{
      result$auc = NA
    }
    
    phi1_true_tmp = phi1_all_true[1:(r*(q+1))]
    if(best){
      phi1 = model$beta[,idx]
      phi1_tmp = phi1_all[1:(r*(q+1))]
      
      result$TP = sum((phi1_tmp!=0)*(phi1_true_tmp!=0))
      result$FP = sum((phi1_tmp!=0)*(phi1_true_tmp==0))
      # main effect
      phi1_main = phi1_tmp[1:r]
      phi1_true_main = phi1_true_tmp[1:r]
      result$TP_main = sum((phi1_main!=0)*(phi1_true_main!=0))
      result$FP_main = sum((phi1_main!=0)*(phi1_true_main==0))
      # interactions
      phi1_inter = phi1_tmp[(r+1):(r*(q+1))]
      phi1_true_inter = phi1_true_tmp[(r+1):(r*(q+1))]
      result$TP_inter = sum((phi1_inter!=0)*(phi1_true_inter!=0))
      result$FP_inter = sum((phi1_inter!=0)*(phi1_true_inter==0))
      
      result$FNR = sum((phi1_tmp==0)*(phi1_true_tmp!=0))/sum(phi1_true_tmp!=0)
      result$FPR = sum((phi1_tmp!=0)*(phi1_true_tmp==0))/sum(phi1_true_tmp==0)
      result$FR = result$FNR+result$FPR
      result$MSE1 = sum((phi1-phi1_all_true)^2)
      result$MSE2 = mean((A1_all%*%(phi1-phi1_all_true))^2)
    }else{
      result$FNR = rep(NA,n.lambda)
      result$FPR = rep(NA,n.lambda)
      result$FR = rep(NA,n.lambda)
      result$MSE1 = rep(NA,n.lambda)
      result$MSE2 = rep(NA,n.lambda)
      result$TP = rep(NA,n.lambda)
      result$FP = rep(NA,n.lambda)
      result$TP_main = rep(NA,n.lambda)
      result$FP_main = rep(NA,n.lambda)
      result$TP_inter = rep(NA,n.lambda)
      result$FP_inter = rep(NA,n.lambda)
      for(k in 1:n.lambda){
        phi1 = model$beta[,k]
        phi1_tmp = phi1[1:(r*(q+1))]
        
        result$TP[k] = sum((phi1_tmp!=0)*(phi1_true_tmp!=0))
        result$FP[k] = sum((phi1_tmp!=0)*(phi1_true_tmp==0))
        # main effect
        phi1_main = phi1_tmp[1:r]
        phi1_true_main = phi1_true_tmp[1:r]
        result$TP_main[k] = sum((phi1_main!=0)*(phi1_true_main!=0))
        result$FP_main[k] = sum((phi1_main!=0)*(phi1_true_main==0))
        # interactions
        phi1_inter = phi1_tmp[(r+1):(r*(q+1))]
        phi1_true_inter = phi1_true_tmp[(r+1):(r*(q+1))]
        result$TP_inter[k] = sum((phi1_inter!=0)*(phi1_true_inter!=0))
        result$FP_inter[k] = sum((phi1_inter!=0)*(phi1_true_inter==0))
        
        result$FNR[k] = sum((phi1_tmp==0)*(phi1_true_tmp!=0))/sum(phi1_true_tmp!=0)
        result$FPR[k] = sum((phi1_tmp!=0)*(phi1_true_tmp==0))/sum(phi1_true_tmp==0)
        result$FR[k] = result$FNR[k]+result$FPR[k]
        
        result$MSE1[k] = sum((phi1-phi1_all_true)^2)
        result$MSE2[k] = mean((A1_all%*%(phi1-phi1_all_true))^2)
      }
    }
  }

  if(!is.null(Y_tes)){
    A1_tes = Z_tes
    for(k in 1:q){
      W = E_tes[,k]*Z_tes
      A1_tes = cbind(A1_tes,W)
    }
    A1_tes_all = cbind(A1_tes,E_tes)
    Y_tes_ord = order(Y_tes,decreasing = T)
    Y_tes_ord_raw = Y_tes[Y_tes_ord]
    N_tes_ord = N_tes[Y_tes_ord]
    A1_tes_ord = A1_tes_all[Y_tes_ord,]
    if(best){
      phi1_all = model$beta[,idx]
      #pl_tes = partial_loglike_val(phi1,Y,N_,A1,Y_tes,N_tes,A1_tes)
      pl_tes = partial_loglike(phi1_all,N_tes_ord,A1_tes_ord,Y_tes_ord_raw)
      C_index_tes = C_index(phi1_all,N_tes_ord,A1_tes_ord,Y_tes_ord_raw)
    }else{
      pl_tes = rep(NA,n.lambda)
      C_index_tes = rep(NA,n.lambda)
      for(k in 1:n.lambda){
        phi1_all = model$beta[,k]
        pl_tes[k] = partial_loglike(phi1_all,N_tes_ord,A1_tes_ord,Y_tes_ord_raw)
        C_index_tes[k] = C_index(phi1_all,N_tes_ord,A1_tes_ord,Y_tes_ord_raw)
      }
    }
    result$pl_tes = pl_tes
    result$C_index_tes = C_index_tes
  }
  return(result)
}

#### model fitting (marginal analysis) ####
# ALT_3: marginal analysis  
ALT_3 <- function(E,Z,Y,N_,phi1_all_true=NULL,method="bonferroni",thresh=0.5)
{
  n = dim(Z)[1]
  q = dim(E)[2]
  r = dim(Z)[2]
  response = Surv(Y, N_)  
  #H = matrix(NA,q+1,r)
  H_pvalue = matrix(NA,q+1,r)
  for(i in 1:r){
    A = matrix(NA,n,2*q+1)
    A[,1] = Z[,i]
    A[,2:(q+1)] = E*Z[,i]
    A[,(q+2):(2*q+1)] = E
    model = summary(coxph(response ~ A))
    #H[,i] = model$coefficients[1:(q+1),1]
    H_pvalue[,i] = model$coefficients[1:(q+1),5]
  }
  phi1_pvalue = as.vector(t(H_pvalue))
  phi1_pvalue_adj = p.adjust(phi1_pvalue,method=method) # multiple test adjust
  phi1_all_sig = (phi1_pvalue_adj<=thresh)  # selected or not
  phi1_all = c(phi1_all_sig, rep(T,q))
  
  result = list()
  result$phi1_all = phi1_all
  if(!is.null(phi1_all_true)){
    phi1_true_tmp = phi1_all_true[1:((q+1)*r)]
    phi1_true_sig = as.vector(phi1_true_tmp!=0)
    result$auc = cal.auc(phi1_true_sig,1-phi1_pvalue)
    result$TP = sum(phi1_all_sig*(phi1_true_tmp!=0))
    result$FP = sum(phi1_all_sig*(phi1_true_tmp==0))
    # main effect
    phi1_sig_main = phi1_all_sig[1:r]
    phi1_true_main = phi1_true_tmp[1:r]
    result$TP_main = sum(phi1_sig_main*(phi1_true_main!=0))
    result$FP_main = sum(phi1_sig_main*(phi1_true_main==0))
    # interactions
    phi1_sig_inter = phi1_all_sig[(r+1):(r*(q+1))]
    phi1_true_inter = phi1_true_tmp[(r+1):(r*(q+1))]
    result$TP_inter = sum(phi1_sig_inter*(phi1_true_inter!=0))
    result$FP_inter = sum(phi1_sig_inter*(phi1_true_inter==0))
    
    result$FNR = sum((phi1_all_sig==0)*(phi1_true_tmp!=0))/sum(phi1_true_tmp!=0)
    result$FPR = sum(phi1_all_sig*(phi1_true_tmp==0))/sum(phi1_true_tmp==0)
    result$FR =  result$FNR+result$FPR
  }
  return(result)
}




