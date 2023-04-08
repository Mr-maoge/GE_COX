library(pROC)



cal.auc <- function(y,y_hat)
{
  as.double(pROC::auc(roc(y,y_hat,quiet = T)))
}

result.table <- function(r.mean,r.std)
{
  n = dim(r.mean)[1]
  p = dim(r.mean)[2]
  r = matrix(" ",n,p)
  for(i in 1:n){
    for(j in 1:p){
      r[i,j] = sprintf("%.3f(%.3f)",r.mean[i,j],r.std[i,j])
    }
  }
  r
}

seq_log <- function(minn,maxn,length.out)
{
  return(exp(seq(log(minn),log(maxn),length.out = length.out)))
}

runif_section <- function(n,locs)
{
  num_sec = length(locs)%/%2
  sep = rep(0,num_sec+1)
  for(i in 1:num_sec){
    sep[i+1] = sep[i]+locs[2*i]-locs[2*i-1]
  }
  res = runif(n,min=0,max=sep[num_sec+1])
  for(k in 1:n){
    val = res[k]
    for(i in 1:num_sec){
      if(val>=sep[i]&val<sep[i+1]){
        res[k] = locs[2*i-1]+val-sep[i]
      }
    }
  }
  return(res)
}

summary.coefs <- function(coef_model,p,q)
{
  tab3 = matrix(coef_model[1:((q+1)*p)],p,q+1)
  tab3 = rbind(c(0,coef_model[((q+1)*p+1):((q+1)*p+q)]),tab3)
  row.names(tab3) = c("Main.E",colnames(X))
  colnames(tab3) = c("Main",colnames(E))
  #tab3 = tab3[which(apply(tab3,1,sum)>0),]
  return(tab3)
}

