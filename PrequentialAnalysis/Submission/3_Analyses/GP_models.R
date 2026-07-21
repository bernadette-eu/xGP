# Load required libraries
library(scoringRules)  # for CRPS and log score
library(Metrics)       # for RMSE and MAE
library(bayesplot)
library(rstan)
library(ggplot2)
library(loo)
library(MASS)
library(dplyr)

df  <- read.csv("...//1_Data//rt1.csv")
df2 <- read.csv("...//1_Data//rt2.csv")

df$date <- as.Date(df$date)

GR<-df[df$country=='Greece',]
UK<-df[df$country=='United Kingdom',]
BEL<-df2[df2$country=='Belgium',]
ITA<-df2[df2$country=='Italy',]
SPAIN<-df2[df2$country=='Spain',]
GER<-df[df$country=='Germany',]
POL <-df[df$country=='Poland',]
DEN <-df[df$country=='Denmark',]
FIN<-df[df$country=='Finland',]
NOR<-df[df$country=='Norway',]
LUX<-df[df$country=='Luxembourg',]
LITH<-df[df$country=='Lithuania',]


min_date<-min(c(GR$date,UK$date,POL$date,DEN$date))
max_date<-max(c(GR$date,UK$date,POL$date,DEN$date))

start_date<-as.Date('2020-03-16')
dates_train<-list()
dates_test<-list()
logRt_train<-list()
logRt_test<-list()
Times_train<-list()
Times_test<-list()
dim_Y<-list()
N_preq<-8
horizon_of_study<-180
M<-3

data_xGPs<-list()
data_iGP <-list()

for(i in 1:N_preq){
  dates_train[[i]]<-seq(from = start_date, to = start_date+horizon_of_study+7*(i-1), by = "day")
  dates_test[[i]]<-seq(from = start_date+horizon_of_study+7*(i-1)+1, to = start_date+horizon_of_study+7*(i), by = "day")
  logRt_train[[i]]<-c(log(GR$median[which(GR$date %in% dates_train[[i]])]),log(UK$median[which(UK$date %in% dates_train[[i]])]),log(GER$median[which(GER$date %in% dates_train[[i]])]))
  logRt_test[[i]]<-c(log(GR$median[which(GR$date %in% dates_test[[i]])]),log(UK$median[which(UK$date %in% dates_test[[i]])]),log(GER$median[which(GER$date %in% dates_test[[i]])]))
  Times_train[[i]]<-rep(1:length(log(GR$median[which(GR$date %in% dates_train[[i]])])),M)
  Times_test[[i]]<-rep((length(log(GR$median[which(GR$date %in% dates_train[[i]])]))+1):(length(log(GR$median[which(GR$date %in% dates_train[[i]])]))+7),M)
  dim_Y[[i]]<-c(length(log(GR$median[which(GR$date %in% dates_train[[i]])])),length(log(UK$median[which(UK$date %in% dates_train[[i]])])),length(log(GER$median[which(GER$date %in% dates_train[[i]])])))
  print(paste0(i,' Length logRt: ',length(logRt_train[[i]]),' sum(dim_Y): ',sum(dim_Y[[i]])))
  data_xGPs[[i]]<-list(M=M,N=length(logRt_train[[i]]),dim_Y=dim_Y[[i]],Y=logRt_train[[i]],Times=Times_train[[i]])
  data_iGP[[i]]<-list()
  data_iGP[[i]][['GR']]<-list(N=length(log(GR$median[which(GR$date %in% dates_train[[i]])])),Y=log(GR$median[which(GR$date %in% dates_train[[i]])]),Times=1:length(log(GR$median[which(GR$date %in% dates_train[[i]])])))
  data_iGP[[i]][['UK']]<-list(N=length(log(GR$median[which(GR$date %in% dates_train[[i]])])),Y=log(UK$median[which(UK$date %in% dates_train[[i]])]),Times=1:length(log(UK$median[which(UK$date %in% dates_train[[i]])])))
  data_iGP[[i]][['GER']]<-list(N=length(log(GER$median[which(GER$date %in% dates_train[[i]])])),Y=log(GER$median[which(GER$date %in% dates_train[[i]])]),Times=1:length(log(GER$median[which(GR$date %in% dates_train[[i]])])))
}


iterations<-20000
warm_up_iter<-iterations/2
thin<-10
fit_xGPs<-list()
fit_iGP<-list()
fit_iGPs<-list()
fit_mxGPs<-list()
fit_xGPs_without_M<-list()
xGPs_mdl<- stan_model(file = "//2_Stan_files//exchangeableGPs_squared_exponential_marginal - without kronecker.stan")
xGPs_without_M_mdl<- stan_model(file = "//2_Stan_files//exchangeableGPs_withoutM_squared_exponential_marginal - without kronecker.stan")
iGP_mdl<-stan_model(file = "//2_Stan_files//indepedentGP_squared_exponential_marginal.stan")
iGPs_mdl<-stan_model(file = "//2_Stan_files//indepedentGPs_squared_exponential_marginal - without kronecker.stan")
mxGPs_mdl<- stan_model(file = "//2_Stan_files//mexchangeableGPs_squared_exponential_marginal - without kronecker.stan")


#do from week 6 onwards for all the models
for(i in 1:N_preq){
  print(i)
  fit_xGPs_without_M[[i]]<-sampling(
    object=xGPs_without_M_mdl,
    data = data_xGPs[[i]],    # named list of data
    chains = 3,             # number of Markov chains
    warmup = warm_up_iter,          # number of warmup iterations per chain
    iter = iterations,            # total number of iterations per chain
    thin=thin,           # total number of iterations per chain
    cores  = 3,              # number of cores (could use one per chain)
    refresh = 10,             # per how many iterations progress is shown
    control=list(max_treedepth =15,adapt_delta=0.8),
    pars=c('length_scale_A','sigma_A','sigma_e'),
    init='random' )
  save(fit_xGPs_without_M,file=paste0("//2_Stan_files//Rt multiple countries\\saved fits\\preq_analysis_fit_xGPs_without_M_",iterations,"iter.RData"))
  print('xGPs without M Done!')
  print(max(apply(get_elapsed_time(fit_xGPs_without_M[[i]]),1,sum))/3600)
}

times_xGPS_without_M<-matrix(NA,N_preq,1)
for(i in 1:N_preq){
  print(i)
  print('xGPS without M')
  print(max(apply(get_elapsed_time(fit_xGPs_without_M[[i]]),1,sum))/3600)
  times_xGPS_without_M[i,1]<-max(apply(get_elapsed_time(fit_xGPs_without_M[[i]]),1,sum))/3600
}

times_mxGPS<-matrix(NA,N_preq,1)
times_xGPS<-matrix(NA,N_preq,1)
times_iGPS<-matrix(NA,N_preq,3)
times_iGPS_alltogether<-matrix(NA,N_preq,1)
colnames(times_iGPS)<-c('GR','UK','GER')

for(i in 1:N_preq){
  print(i)
  print('mxGPS')
  print(max(apply(get_elapsed_time(fit_mxGPs[[i]]),1,sum))/3600)
  times_mxGPS[i,1]<-max(apply(get_elapsed_time(fit_mxGPs[[i]]),1,sum))/3600
  print('xGPS')
  print(max(apply(get_elapsed_time(fit_xGPs[[i]]),1,sum))/3600)
  times_xGPS[i,1]<-max(apply(get_elapsed_time(fit_xGPs[[i]]),1,sum))/3600
  print('iGPS')
  print(max(apply(get_elapsed_time(fit_iGPs[[i]]),1,sum))/3600)
  times_iGPS_alltogether[i,1]<-max(apply(get_elapsed_time(fit_iGPs[[i]]),1,sum))/3600
  print('GR')
  print(max(apply(get_elapsed_time(fit_iGP[[i]][['GR']]),1,sum))/3600)
  times_iGPS[i,1]<-max(apply(get_elapsed_time(fit_iGP[[i]][['GR']]),1,sum))/3600
  print('UK')
  times_iGPS[i,2]<-max(apply(get_elapsed_time(fit_iGP[[i]][['UK']]),1,sum))/3600
  
  print(max(apply(get_elapsed_time(fit_iGP[[i]][['UK']]),1,sum))/3600)
  print('GER')
  print(max(apply(get_elapsed_time(fit_iGP[[i]][['GER']]),1,sum))/3600)
  times_iGPS[i,3]<-max(apply(get_elapsed_time(fit_iGP[[i]][['GER']]),1,sum))/3600
  
  
}

#COMPARE WAIC AND LOO
set.seed(1)

#xGPs predictions and generated quantities
length_scale_A_mx<-list()
sigma_A_mx<-list()
length_scale_M_mx<-list()
sigma_M_mx<-list()
sigma_e_mx<-list()
f_mxGPs<-list()
predictions_mxGPs<-list()

start_time<-Sys.time()
for(i in 1:N_preq){
  print('mxGPs')
  print(i)
  matrix_of_draws_mx <- as.data.frame(fit_mxGPs[[i]])
  
  length_scale_A_mx[[i]] <-matrix_of_draws_mx[,c('length_scale_A[1]','length_scale_A[2]','length_scale_A[3]')]
  sigma_A_mx[[i]]     <-matrix_of_draws_mx[,c('sigma_A[1]','sigma_A[2]','sigma_A[3]')]     
  length_scale_M_mx[[i]]<-matrix_of_draws_mx[,'length_scale_M']
  sigma_M_mx[[i]]      <-matrix_of_draws_mx[,'sigma_M']
  sigma_e_mx[[i]]   <-matrix_of_draws_mx[,'sigma_e']
  rm(matrix_of_draws)
  K_mxGPs<-compute_K_mxGP(data_xGPs[[i]]$dim_Y,data_xGPs[[i]]$Times,sigma_A_mx[[i]],length_scale_A_mx[[i]],sigma_M_mx[[i]],length_scale_M_mx[[i]],NULL)
  print('K calculated')
  print(object.size(K_mxGPs))
  
  f_mxGPs[[i]]<-compute_f_GPs(K_mxGPs,sigma_e_mx[[i]],data_xGPs[[i]]$Y,NULL)
  print('f calculated')
  
  predictions_mxGPs[[i]]<-matrix(NA,length(Times_test[[i]]),length(sigma_e_mx[[i]]))
  times_M<-cbind(Times_test[[i]],rep(1:M, each = length(Times_test[[i]])/M))
  for(j in 1:dim(times_M)[1]){
    predictions_mxGPs[[i]][j,]<-compute_Y_pred_mxGP(times_M[j,1],times_M[j,2],data_xGPs[[i]]$dim_Y,data_xGPs[[i]]$Times,f_mxGPs[[i]],K_mxGPs,data_xGPs[[i]]$Y,sigma_A_mx[[i]], length_scale_A_mx[[i]], 
                                                  sigma_M_mx[[i]], length_scale_M_mx[[i]],sigma_e_mx[[i]],fun=NULL)
  }
  print('predictions done')
  rm(K_mxGPs)
  generated_quantites_mxGPs<-list(f_mxGPs)
  save(generated_quantites_mxGPs,file=paste0("//2_Stan_files//Rt multiple countries\\saved fits\\preq_analysis_generated_quantities_mxGPs_",iterations,"iter.RData"))
  
  save(predictions_mxGPs,file=paste0("//2_Stan_files//Rt multiple countries\\saved fits\\preq_analysis_predictions_mxGPs_",iterations,"iter.RData"))
  
}
end_time<-Sys.time()
print(end_time-start_time)

rm(length_scale_A_mx)
rm(sigma_A_mx)
rm(length_scale_M_mx)
rm(sigma_M_mx)
rm(sigma_e_mx)
          
#xGPs predictions and generated quantities
length_scale_A_x<-list()
sigma_A_x<-list()
length_scale_M_x<-list()
sigma_M_x<-list()
sigma_e_x<-list()
f_xGPs<-list()
predictions_xGPs<-list()

start_time<-Sys.time()
for(i in 1:N_preq){
  print('xGPs')
  print(i)

  matrix_of_draws_x <- as.data.frame(fit_xGPs[[i]])
  
  length_scale_A_x[[i]] <-matrix_of_draws_x[,'length_scale_A']
  sigma_A_x[[i]]     <-matrix_of_draws_x[,'sigma_A']   
  length_scale_M_x[[i]]<-matrix_of_draws_x[,'length_scale_M']
  sigma_M_x[[i]]      <-matrix_of_draws_x[,'sigma_M']
  sigma_e_x[[i]]   <-matrix_of_draws_x[,'sigma_e']
  rm(matrix_of_draws)
  K_xGPs<-compute_K_xGP(data_xGPs[[i]]$dim_Y,data_xGPs[[i]]$Times,sigma_A_x[[i]],length_scale_A_x[[i]],sigma_M_x[[i]],length_scale_M_x[[i]],NULL)
  print('K calculated')
  print(object.size(K_xGPs))

  f_xGPs[[i]]<-compute_f_GPs(K_xGPs,sigma_e_x[[i]],data_xGPs[[i]]$Y,NULL)
  print('f calculated')
  
  predictions_xGPs[[i]]<-matrix(NA,length(Times_test[[i]]),length(sigma_e_x[[i]]))
  times_M<-cbind(Times_test[[i]],rep(1:M, each = length(Times_test[[i]])/M))
  for(j in 1:dim(times_M)[1]){
    predictions_xGPs[[i]][j,]<-compute_Y_pred_xGP(times_M[j,1],times_M[j,2],data_xGPs[[i]]$dim_Y,data_xGPs[[i]]$Times,f_xGPs[[i]],K_xGPs,data_xGPs[[i]]$Y,sigma_A_x[[i]], length_scale_A_x[[i]], 
                                                  sigma_M_x[[i]], length_scale_M_x[[i]],sigma_e_x[[i]],fun=NULL)
  }
  print('predictions done')
  rm(K_xGPs)
  generated_quantites_xGPs<-list(f_xGPs)
  save(generated_quantites_xGPs,file=paste0("//2_Stan_files//Rt multiple countries\\saved fits\\preq_analysis_generated_quantities_xGPs_",iterations,"iter.RData"))
  
  save(predictions_xGPs,file=paste0("//2_Stan_files//Rt multiple countries\\saved fits\\preq_analysis_predictions_xGPs_",iterations,"iter.RData"))
  
}
end_time<-Sys.time()
print(end_time-start_time)

rm(length_scale_A_x)
rm(sigma_A_x)
rm(length_scale_M_x)
rm(sigma_M_x)
rm(sigma_e_x)


######################################################################
#
# iGPs predictions and generated quantities
#
######################################################################
length_scale<-list()
sigma<-list()

sigma_e<-list()
f_iGP<-list()
predictions_iGP<-list()

start_time<-Sys.time()
for(i in 1:N_preq){
  print(i)
  length_scale[[i]]<-list()
  sigma[[i]]<-list()
  sigma_e[[i]]<-list()
  predictions_iGP[[i]]<-list()
  matrix_of_draws1 <- as.data.frame(fit_iGP[[i]][['GR']])
  matrix_of_draws2 <- as.data.frame(fit_iGP[[i]][['UK']])
  matrix_of_draws3<- as.data.frame(fit_iGP[[i]][['GER']])
  length_scale[[i]][['GR']] <-matrix_of_draws1[,'length_scale']
  sigma[[i]][['GR']] <-matrix_of_draws1[,'sigma']
  sigma_e[[i]][['GR']] <-matrix_of_draws1[,'sigma_e']
  length_scale[[i]][['UK']] <-matrix_of_draws2[,'length_scale']
  sigma[[i]][['UK']] <-matrix_of_draws2[,'sigma']
  sigma_e[[i]][['UK']] <-matrix_of_draws2[,'sigma_e']
  length_scale[[i]][['GER']] <-matrix_of_draws3[,'length_scale']
  sigma[[i]][['GER']] <-matrix_of_draws3[,'sigma']
  sigma_e[[i]][['GER']] <-matrix_of_draws3[,'sigma_e']
  rm(matrix_of_draws1)
  rm(matrix_of_draws2)
  rm(matrix_of_draws3)
  K_iGP_GR<-compute_K_iGP(data_iGP[[i]][['GR']]$Times,sigma[[i]][['GR']],length_scale[[i]][['GR']],NULL)
  K_iGP_UK<-compute_K_iGP(data_iGP[[i]][['UK']]$Times,sigma[[i]][['UK']],length_scale[[i]][['UK']],NULL)
  K_iGP_GER<-compute_K_iGP(data_iGP[[i]][['GER']]$Times,sigma[[i]][['GER']],length_scale[[i]][['GER']],NULL)
  
  print('Ks calculated')

  f_iGP[[i]]<-list()
  start_time3<-Sys.time()
  f_iGP[[i]][['GR']]<-compute_f_iGP(K_iGP_GR,sigma_e[[i]][['GR']],data_iGP[[i]][['GR']]$Y,NULL)
  f_iGP[[i]][['UK']]<-compute_f_iGP(K_iGP_UK,sigma_e[[i]][['UK']],data_iGP[[i]][['UK']]$Y,NULL)
  f_iGP[[i]][['GER']]<-compute_f_iGP(K_iGP_GER,sigma_e[[i]][['GER']],data_iGP[[i]][['GER']]$Y,NULL)
  
  print('fs calculated')
  
  generated_quantites_iGP<-list(f_iGP)
  save(generated_quantites_iGP,file=paste0("//2_Stan_files//Rt multiple countries\\saved fits\\preq_analysis_generated_quantities_iGP_",iterations,"iter.RData"))
  
  predictions_iGP[[i]][['GR']]<-matrix(NA,length(Times_test[[i]]),length(sigma_e[[i]][['GR']]))
  predictions_iGP[[i]][['UK']]<-matrix(NA,length(Times_test[[i]]),length(sigma_e[[i]][['UK']]))
  predictions_iGP[[i]][['GER']]<-matrix(NA,length(Times_test[[i]]),length(sigma_e[[i]][['GER']]))
  #times_M<-cbind(Times_test[[i]],rep(1:M, each = length(Times_test[[i]])/M))
  for(j in 1:(length(Times_test[[i]])/M)){
    predictions_iGP[[i]][['GR']][j,]<-compute_Y_pred_iGP(Times_test[[i]][j],data_iGP[[i]][['GR']]$Times,f_iGP[[i]][['GR']],K_iGP_GR,data_iGP[[i]][['GR']]$Y,sigma[[i]][['GR']], length_scale[[i]][['GR']],sigma_e[[i]][['GR']],fun=NULL)
    predictions_iGP[[i]][['UK']][j,]<-compute_Y_pred_iGP(Times_test[[i]][j],data_iGP[[i]][['UK']]$Times,f_iGP[[i]][['UK']],K_iGP_UK,data_iGP[[i]][['UK']]$Y,sigma[[i]][['UK']], length_scale[[i]][['UK']],sigma_e[[i]][['UK']],fun=NULL)
    predictions_iGP[[i]][['GER']][j,]<-compute_Y_pred_iGP(Times_test[[i]][j],data_iGP[[i]][['GER']]$Times,f_iGP[[i]][['GER']],K_iGP_GER,data_iGP[[i]][['GER']]$Y,sigma[[i]][['GER']], length_scale[[i]][['GER']],sigma_e[[i]][['GER']],fun=NULL)
  }
  
  save(predictions_iGP,file=paste0("//2_Stan_files//Rt multiple countries\\saved fits\\preq_analysis_predictions_iGP_",iterations,"iter.RData"))
  
  
  rm(K_iGP_GR)
  rm(K_iGP_UK)
  rm(K_iGP_GER)
}
end_time<-Sys.time()
print(end_time-start_time)

######################################################################
#
# xGPs without M predictions and generated quantities
#
######################################################################
length_scale_A_x_without_M<-list()
sigma_A_x_without_M<-list()
sigma_e_x_without_M<-list()
f_xGPs_without_M<-list()
predictions_xGPs_without_M<-list()

start_time<-Sys.time()
for(i in 1:N_preq){
  print('xGPs')
  print(i)
  
  matrix_of_draws_x_without_M <- as.data.frame(fit_xGPs_without_M[[i]])
  
  length_scale_A_x_without_M[[i]] <-matrix_of_draws_x_without_M[,'length_scale_A']
  sigma_A_x_without_M[[i]]     <-matrix_of_draws_x_without_M[,'sigma_A']   
  sigma_e_x_without_M[[i]]   <-matrix_of_draws_x_without_M[,'sigma_e']
  rm(matrix_of_draws_without_M)
  K_xGPs_without_M<-compute_K_xGP_without_M(data_xGPs[[i]]$dim_Y,data_xGPs[[i]]$Times,sigma_A_x_without_M[[i]],length_scale_A_x_without_M[[i]],NULL)
  print('K calculated')
  print(object.size(K_xGPs_without_M))
  
  f_xGPs_without_M[[i]]<-compute_f_GPs(K_xGPs_without_M,sigma_e_x_without_M[[i]],data_xGPs[[i]]$Y,NULL)
  print('f calculated')
  
  predictions_xGPs_without_M[[i]]<-matrix(NA,length(Times_test[[i]]),length(sigma_e_x_without_M[[i]]))
  times_M<-cbind(Times_test[[i]],rep(1:M, each = length(Times_test[[i]])/M))
  for(j in 1:dim(times_M)[1]){
    predictions_xGPs_without_M[[i]][j,]<-compute_Y_pred_xGP_without_M(times_M[j,1],times_M[j,2],data_xGPs[[i]]$dim_Y,data_xGPs[[i]]$Times,f_xGPs_without_M[[i]],K_xGPs_without_M,data_xGPs[[i]]$Y,sigma_A_x_without_M[[i]], length_scale_A_x_without_M[[i]], 
                                                                      sigma_e_x_without_M[[i]],fun=NULL)
  }
  print('predictions done')
  rm(K_xGPs_without_M)
  generated_quantites_xGPs_without_M<-list(f_xGPs_without_M)
  save(generated_quantites_xGPs_without_M,file=paste0("//2_Stan_files//Rt multiple countries\\saved fits\\preq_analysis_generated_quantities_xGPs_without_M_",iterations,"iter.RData"))
  
  save(predictions_xGPs_without_M,file=paste0("//2_Stan_files//Rt multiple countries\\saved fits\\preq_analysis_predictions_xGPs_without_M_",iterations,"iter.RData"))
  
}
end_time<-Sys.time()
print(end_time-start_time)

rm(length_scale_A_x_without_M)
rm(sigma_A_x_without_M)
rm(sigma_e_x_without_M)
