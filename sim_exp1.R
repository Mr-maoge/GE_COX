#setwd("Z:\\ARMI\\code_topic1\\code_submit")
library(foreach)
library(doParallel)
source("sim_funcs.R")


#### settings ####
sim_scheme = read.csv("0_sim_scheme_topic1_5.csv")
result_folder = "./result/"

#### parallel simulation ####
cl=makeCluster(23, outfile="")
registerDoParallel(cl)
repeat_times = 20

for(i in 1:24){
  p = 600
  n = 350
  signal = 0.8
  censoring = 0.4
  cor1_type = sim_scheme$cor1_type[i]
  cor2_type = sim_scheme$cor2_type[i]
  rho1 = sim_scheme$rho1[i]
  rho2 = sim_scheme$rho2[i]
  n_main = sim_scheme$n_main[i]
  n_interaction = sim_scheme$n_interaction[i]
  
  lambda1_set = seq_log(sim_scheme$lambda1_min[i]*1.5,sim_scheme$lambda1_max[i]*1.5,length.out = 8)
  lambda2_set = seq_log(sim_scheme$lambda1_min[i],sim_scheme$lambda1_max[i],length.out = 6)
  
  set.seed(1)
  seed_each = sample.int(100001, repeat_times)
  start.time = Sys.time()
  cat("==== i=",i," ====","\n")
  print(start.time)
  
  result_all = foreach(seed=seed_each,.combine=sim.combine) %dopar%
    sim.example1(seed = seed,n=n,p=p,n_main=n_main,n_interactions=n_interaction,
                 E_is_continuous=E_is_continuous,
                 rho1=rho1,rho2=rho2,cor1_type=cor1_type,cor2_type=cor2_type,signal=signal,
                 metric = "ebic",censoring=censoring,
                 lambda1_set = lambda1_set,lambda2_set=lambda2_set,
                 mode="normal"
    )
  
  end.time = Sys.time()
  print(end.time)
  print(end.time-start.time)
  
  result_mean = apply(result_all$result,2:3,mean,na.rm=T)
  result_std = apply(result_all$result,2:3,sd,na.rm=T)
  
  metrics = c("AUC","FNR","FPR","FNR+FPR",
              "TP","FP","TP_main","FP_main","TP_inter","FP_inter",
              "MSE1","MSE2","pl_test","C_index_tes")
  methods = c("Proposed","MCP","grMCP","Marginal")
  
  result_out0 = result.table(result_mean,result_std)
  result_out1 = data.frame(result_out0)
  row.names(result_out1) = methods
  colnames(result_out1) = metrics
  filename = sprintf("%s/sim1_p%d_%d.csv",result_folder,p,i)
  write.csv(result_out1,file=filename)
  filename2 = sprintf("%s/sim1_p%d_%d_tuning.csv",result_folder,p,i)
  write.csv(result_all$tunings,file=filename2)
  print(result_out1)
}

stopCluster(cl)
