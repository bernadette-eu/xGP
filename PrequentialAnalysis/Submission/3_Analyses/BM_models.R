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

start_date  <-as.Date('2020-03-16')
dates_train <-list()
dates_test  <-list()
logRt_train <-list()
logRt_test  <-list()
Times_train <-list()
Times_test  <- list()
dim_Y       <-list()
N_preq      <- 8
horizon_of_study <- 180

M <- 3

data_xGPs <- list()
data_iGP  <- list()

for(i in 1:N_preq){
  dates_train[[i]]<-seq(from = start_date, to = start_date+horizon_of_study+7*(i-1), by = "day")
  dates_test[[i]]<-seq(from = start_date+horizon_of_study+7*(i-1)+1, to = start_date+horizon_of_study+7*(i), by = "day")
  logRt_train[[i]]<-c(log(GR$median[which(GR$date %in% dates_train[[i]])]),
                      log(UK$median[which(UK$date %in% dates_train[[i]])]),
                      log(GER$median[which(GER$date %in% dates_train[[i]])]))
  logRt_test[[i]]<-c(log(GR$median[which(GR$date %in% dates_test[[i]])]),log(UK$median[which(UK$date %in% dates_test[[i]])]),log(GER$median[which(GER$date %in% dates_test[[i]])]))
  Times_train[[i]]<-rep(1:length(log(GR$median[which(GR$date %in% dates_train[[i]])])),M)
  Times_test[[i]]<-rep((length(log(GR$median[which(GR$date %in% dates_train[[i]])]))+1):(length(log(GR$median[which(GR$date %in% dates_train[[i]])]))+7),M)
  dim_Y[[i]]<-c(length(log(GR$median[which(GR$date %in% dates_train[[i]])])),length(log(UK$median[which(UK$date %in% dates_train[[i]])])),length(log(GER$median[which(GER$date %in% dates_train[[i]])])))
  print(paste0(i,' Length logRt: ',length(logRt_train[[i]]),' sum(dim_Y): ',sum(dim_Y[[i]])))
  data_xGPs[[i]]<-list(M=M,N=length(logRt_train[[i]]),dim_Y=dim_Y[[i]],Y=logRt_train[[i]],Times=Times_train[[i]])
  data_iGP[[i]]<-list()
  data_iGP[[i]][['GR']]<-list(N=length(log(GR$median[which(GR$date %in% dates_train[[i]])])),
                              Y=log(GR$median[which(GR$date %in% dates_train[[i]])]),
                              Times=1:length(log(GR$median[which(GR$date %in% dates_train[[i]])])))
  data_iGP[[i]][['UK']]<-list(N=length(log(GR$median[which(GR$date %in% dates_train[[i]])])),
                              Y=log(UK$median[which(UK$date %in% dates_train[[i]])]),
                              Times=1:length(log(UK$median[which(UK$date %in% dates_train[[i]])])))
  data_iGP[[i]][['GER']]<-list(N=length(log(GER$median[which(GER$date %in% dates_train[[i]])])),
                               Y=log(GER$median[which(GER$date %in% dates_train[[i]])]),
                               Times=1:length(log(GER$median[which(GR$date %in% dates_train[[i]])])))
}


iterations<-20000
warm_up_iter<-iterations/2
thin<-10
fit_xBMs<-list()
fit_iBM<-list()
#fit_iGPs<-list()
fit_mxBMs<-list()
fit_xBMs_without_M<-list()
xBMs_mdl<- stan_model(file = "C:\\Users\\barbounakis\\OneDrive - aueb.gr\\Post Doc AUEB\\BOURANIS LAMPROS - xGP_project\\petros work\\GP rstan\\exchangeableGPs_brownian_motion_marginal - without kronecker.stan")
xBMs_without_M_mdl<- stan_model(file = "C:\\Users\\barbounakis\\OneDrive - aueb.gr\\Post Doc AUEB\\BOURANIS LAMPROS - xGP_project\\petros work\\GP rstan\\exchangeableGPs_withoutM_brownian_motion_marginal - without kronecker.stan")
iBM_mdl<-stan_model(file = "C:\\Users\\barbounakis\\OneDrive - aueb.gr\\Post Doc AUEB\\BOURANIS LAMPROS - xGP_project\\petros work\\GP rstan\\indepedentGP_brownian_motion_marginal.stan")
#iGPs_mdl<-stan_model(file = "C:\\Users\\barbounakis\\OneDrive - aueb.gr\\Post Doc AUEB\\BOURANIS LAMPROS - xGP_project\\petros work\\GP rstan\\indepedentGPs_brownian_motion_marginal - without kronecker.stan")
mxBMs_mdl<- stan_model(file = "C:\\Users\\barbounakis\\OneDrive - aueb.gr\\Post Doc AUEB\\BOURANIS LAMPROS - xGP_project\\petros work\\GP rstan\\mexchangeableGPs_brownian_motion_marginal - without kronecker.stan")

fit_xBMs_without_M[[6]]<-sampling(
  object=xBMs_without_M_mdl,
  data = data_xGPs[[6]],    # named list of data
  chains = 3,             # number of Markov chains
  warmup = warm_up_iter,          # number of warmup iterations per chain
  iter = iterations,            # total number of iterations per chain
  thin=thin,           # total number of iterations per chain
  cores  = 3,              # number of cores (could use one per chain)
  refresh = 10,             # per how many iterations progress is shown
  control=list(max_treedepth =15,adapt_delta=0.8),
  pars=c('sigma_A','sigma_e'),
  init='random' )
save(fit_xBMs_without_M,file=paste0("...//preq_analysis_fit_xBMs_without_M_",iterations,"iter.RData"))
print('xBMs without M Done!')
print(max(apply(get_elapsed_time(fit_xBMs_without_M[[6]]),1,sum))/3600)


fit_iBM[[6]]<-list()

fit_iBM[[6]][['GR']]<-sampling(
  object=iBM_mdl,
  data = data_iGP[[6]][['GR']], #    named list of data
  chains = 3,      #        number of Markov chains
  warmup = warm_up_iter,   #        number of warmup iterations per chain
  iter = iterations,    #         total number of iterations per chain
  thin=thin,        #    total number of iterations per chain
  cores  = 3,       #        number of cores (could use one per chain)
  refresh = 10,       #       per how many iterations progress is shown
  control=list(max_treedepth =15,adapt_delta=0.8),
  pars=c('sigma','sigma_e'),
  init='random' )
save(fit_iBM,file=paste0("...//preq_analysis_fit_iBM_",iterations,"iter.RData"))
print('iBM GR Done!')
print(max(apply(get_elapsed_time(fit_iBM[[6]][['GR']]),1,sum))/3600)


fit_iBM[[6]][['UK']]<-sampling(
  object=iBM_mdl,
  data = data_iGP[[6]][['UK']],    # named list of data
  chains = 3,       #       number of Markov chains
  warmup = warm_up_iter,    #       number of warmup iterations per chain
  iter = iterations,    #         total number of iterations per chain
  thin=thin,    #        total number of iterations per chain
  cores  = 3,    #           number of cores (could use one per chain)
  refresh = 10,     #         per how many iterations progress is shown
  control=list(max_treedepth =15,adapt_delta=0.8),
  pars=c('sigma','sigma_e'),
  init='random' )
save(fit_iBM,file=paste0("...//preq_analysis_fit_iBM_",iterations,"iter.RData"))
print('iBM UK Done!')
print(max(apply(get_elapsed_time(fit_iBM[[6]][['UK']]),1,sum))/3600)

fit_iBM[[6]][['GER']]<-sampling(
  object=iBM_mdl,
  data = data_iGP[[6]][['GER']], #    named list of data
  chains = 3,       #       number of Markov chains
  warmup = warm_up_iter,      #     number of warmup iterations per chain
  iter = iterations,    #         total number of iterations per chain
  thin=thin,       #     total number of iterations per chain
  cores  = 3,     #          number of cores (could use one per chain)
  refresh = 10,     #         per how many iterations progress is shown
  control=list(max_treedepth =15,adapt_delta=0.8),
  pars=c('sigma','sigma_e'),
  init='random' )
save(fit_iBM,file=paste0("...//preq_analysis_fit_iBM_",iterations,"iter.RData"))
print('iBM GER Done!')
print(max(apply(get_elapsed_time(fit_iBM[[6]][['GER']]),1,sum))/3600)

print('iBMs Done!')


#do from week 7 onwards for all the models
for(i in 7:N_preq){
  
  fit_mxBMs[[i]]<-sampling(
    object=mxBMs_mdl,
    data = data_xGPs[[i]],    # named list of data
     chains = 3,             # number of Markov chains
     warmup = warm_up_iter,          # number of warmup iterations per chain
     iter = iterations,            # total number of iterations per chain
     thin=thin,           # total number of iterations per chain
     cores  = 3,              # number of cores (could use one per chain)
     refresh = 10,             # per how many iterations progress is shown
     control=list(max_treedepth =15,adapt_delta=0.8),
     pars=c('sigma_A','sigma_M','sigma_e'),
     init='random' )
   save(fit_mxBMs,file=paste0("...//preq_analysis_fit_mxBMs_",iterations,"iter.RData"))
   print('mxBMs Done!')
   print(max(apply(get_elapsed_time(fit_mxBMs[[i]]),1,sum))/3600)
   

    fit_xBMs[[i]]<-sampling(
     object=xBMs_mdl,
      data = data_xGPs[[i]],    # named list of data
     chains = 3,             # number of Markov chains
     warmup = warm_up_iter,          # number of warmup iterations per chain
     iter = iterations,            # total number of iterations per chain
     thin=thin,           # total number of iterations per chain
     cores  = 3,              # number of cores (could use one per chain)
     refresh = 10,             # per how many iterations progress is shown
     control=list(max_treedepth =15,adapt_delta=0.8),
     pars=c('sigma_A','sigma_M','sigma_e'),
     init='random' )
  save(fit_xBMs,file=paste0("...//preq_analysis_fit_xBMs_",iterations,"iter.RData"))
   print('xBMs Done!')
    print(max(apply(get_elapsed_time(fit_xBMs[[i]]),1,sum))/3600)
    
    fit_xBMs_without_M[[i]]<-sampling(
      object=xBMs_without_M_mdl,
      data = data_xGPs[[i]],    # named list of data
      chains = 3,             # number of Markov chains
      warmup = warm_up_iter,          # number of warmup iterations per chain
      iter = iterations,            # total number of iterations per chain
      thin=thin,           # total number of iterations per chain
      cores  = 3,              # number of cores (could use one per chain)
      refresh = 10,             # per how many iterations progress is shown
      control=list(max_treedepth =15,adapt_delta=0.8),
      pars=c('sigma_A','sigma_e'),
      init='random' )
    save(fit_xBMs_without_M,file=paste0("...//preq_analysis_fit_xBMs_without_M_",iterations,"iter.RData"))
    print('xBMs without M Done!')
    print(max(apply(get_elapsed_time(fit_xBMs_without_M[[i]]),1,sum))/3600)
    
    
    
    fit_iBM[[i]]<-list()
   
    fit_iBM[[i]][['GR']]<-sampling(
      object=iBM_mdl,
      data = data_iGP[[i]][['GR']], #    named list of data
      chains = 3,      #        number of Markov chains
      warmup = warm_up_iter,   #        number of warmup iterations per chain
      iter = iterations,    #         total number of iterations per chain
      thin=thin,        #    total number of iterations per chain
      cores  = 3,       #        number of cores (could use one per chain)
      refresh = 10,       #       per how many iterations progress is shown
      control=list(max_treedepth =15,adapt_delta=0.8),
      pars=c('sigma','sigma_e'),
      init='random' )
    save(fit_iBM,file=paste0("...//preq_analysis_fit_iBM_",iterations,"iter.RData"))
    print('iBM GR Done!')
    print(max(apply(get_elapsed_time(fit_iBM[[i]][['GR']]),1,sum))/3600)
    
    
    fit_iBM[[i]][['UK']]<-sampling(
      object=iBM_mdl,
      data = data_iGP[[i]][['UK']],    # named list of data
      chains = 3,       #       number of Markov chains
      warmup = warm_up_iter,    #       number of warmup iterations per chain
      iter = iterations,    #         total number of iterations per chain
      thin=thin,    #        total number of iterations per chain
      cores  = 3,    #           number of cores (could use one per chain)
      refresh = 10,     #         per how many iterations progress is shown
      control=list(max_treedepth =15,adapt_delta=0.8),
      pars=c('sigma','sigma_e'),
      init='random' )
    save(fit_iBM,file=paste0("...//preq_analysis_fit_iBM_",iterations,"iter.RData"))
    print('iBM UK Done!')
    print(max(apply(get_elapsed_time(fit_iBM[[i]][['UK']]),1,sum))/3600)
    
    fit_iBM[[i]][['GER']]<-sampling(
      object=iBM_mdl,
      data = data_iGP[[i]][['GER']], #    named list of data
      chains = 3,       #       number of Markov chains
      warmup = warm_up_iter,      #     number of warmup iterations per chain
      iter = iterations,    #         total number of iterations per chain
      thin=thin,       #     total number of iterations per chain
      cores  = 3,     #          number of cores (could use one per chain)
      refresh = 10,     #         per how many iterations progress is shown
      control=list(max_treedepth =15,adapt_delta=0.8),
      pars=c('sigma','sigma_e'),
      init='random' )
    save(fit_iBM,file=paste0("...//preq_analysis_fit_iBM_",iterations,"iter.RData"))
    print('iBM GER Done!')
    print(max(apply(get_elapsed_time(fit_iBM[[i]][['GER']]),1,sum))/3600)
    
    print('iBMs Done!')
   
}

times_mxBMS<-matrix(NA,N_preq,1)
times_xBMS<-matrix(NA,N_preq,1)
times_xBMS_without_M<-matrix(NA,N_preq,1)
times_iBMS<-matrix(NA,N_preq,3)
colnames(times_iBMS)<-c('GR','UK','GER')

for(i in 1:N_preq){
  print(i)
  print('mxBMS')
  print(max(apply(get_elapsed_time(fit_mxBMs[[i]]),1,sum))/3600)
  times_mxBMS[i,1]<-max(apply(get_elapsed_time(fit_mxBMs[[i]]),1,sum))/3600
  print('xBMS')
  print(max(apply(get_elapsed_time(fit_xBMs[[i]]),1,sum))/3600)
  times_xBMS[i,1]<-max(apply(get_elapsed_time(fit_xBMs[[i]]),1,sum))/3600
  print('xBMS without M')
  print(max(apply(get_elapsed_time(fit_xBMs_without_M[[i]]),1,sum))/3600)
  times_xBMS_without_M[i,1]<-max(apply(get_elapsed_time(fit_xBMs_without_M[[i]]),1,sum))/3600
  print('iBMS')
  print('GR')
  print(max(apply(get_elapsed_time(fit_iBM[[i]][['GR']]),1,sum))/3600)
  times_iBMS[i,1]<-max(apply(get_elapsed_time(fit_iBM[[i]][['GR']]),1,sum))/3600
  print('UK')
  times_iBMS[i,2]<-max(apply(get_elapsed_time(fit_iBM[[i]][['UK']]),1,sum))/3600
  
  print(max(apply(get_elapsed_time(fit_iBM[[i]][['UK']]),1,sum))/3600)
  print('GER')
  print(max(apply(get_elapsed_time(fit_iBM[[i]][['GER']]),1,sum))/3600)
  times_iBMS[i,3]<-max(apply(get_elapsed_time(fit_iBM[[i]][['GER']]),1,sum))/3600
  
  
}

#COMPARE WAIC AND LOO
set.seed(1)

#xBMs predictions and generated quantities
sigma_A_mx<-list()
sigma_M_mx<-list()
sigma_e_mx<-list()
f_mxBMs<-list()
predictions_mxBMs<-list()

start_time<-Sys.time()
for(i in 1:N_preq){
  print('mxBMs')
  print(i)
  matrix_of_draws_mx <- as.data.frame(fit_mxBMs[[i]])
  
  sigma_A_mx[[i]] <- matrix_of_draws_mx[,c('sigma_A[1]','sigma_A[2]','sigma_A[3]')]     
  sigma_M_mx[[i]] <- matrix_of_draws_mx[,'sigma_M']
  sigma_e_mx[[i]] <- matrix_of_draws_mx[,'sigma_e']
  
  rm(matrix_of_draws_mx)
  K_mxBMs<-compute_K_mxBM(data_xGPs[[i]]$dim_Y,data_xGPs[[i]]$Times,sigma_A_mx[[i]],sigma_M_mx[[i]],NULL)
  print('K calculated')
  print(object.size(K_mxBMs))
  
  f_mxBMs[[i]]<-compute_f_BMs(K_mxBMs,sigma_e_mx[[i]],data_xGPs[[i]]$Y,NULL)
  print('f calculated')
  
  predictions_mxBMs[[i]]<-matrix(NA,length(Times_test[[i]]),length(sigma_e_mx[[i]]))
  times_M<-cbind(Times_test[[i]],rep(1:M, each = length(Times_test[[i]])/M))
  for(j in 1:dim(times_M)[1]){
    predictions_mxBMs[[i]][j,]<-compute_Y_pred_mxBM(times_M[j,1],times_M[j,2],data_xGPs[[i]]$dim_Y,data_xGPs[[i]]$Times,f_mxBMs[[i]],K_mxBMs,data_xGPs[[i]]$Y,sigma_A_mx[[i]],
                                                  sigma_M_mx[[i]],sigma_e_mx[[i]],fun=NULL)
  }
  print('predictions done')
  rm(K_mxBMs)
  generated_quantites_mxBMs<-list(f_mxBMs)
  save(generated_quantites_mxBMs,file=paste0("...//preq_analysis_generated_quantities_mxBMs_",iterations,"iter.RData"))
  
  save(predictions_mxBMs,file=paste0("...//preq_analysis_predictions_mxBMs_",iterations,"iter.RData"))
  
}
end_time<-Sys.time()
print(end_time-start_time)

rm(sigma_A_mx)
rm(sigma_M_mx)
rm(sigma_e_mx)
          
#xBMs predictions and generated quantities
sigma_A_x <- list()
sigma_M_x <- list()
sigma_e_x <- list()
f_xBMs    <- list()
predictions_xBMs <- list()
# f_M_xBMs<-list()

start_time <- Sys.time()

for(i in 1:N_preq){
  print('xBMs')
  print(i)

  matrix_of_draws_x <- as.data.frame(fit_xBMs[[i]])
  
  sigma_A_x[[i]]     <-matrix_of_draws_x[,'sigma_A']   
  sigma_M_x[[i]]      <-matrix_of_draws_x[,'sigma_M']
  sigma_e_x[[i]]   <-matrix_of_draws_x[,'sigma_e']
  rm(matrix_of_draws_x)
  K_xBMs<-compute_K_xBM(data_xGPs[[i]]$dim_Y,data_xGPs[[i]]$Times,sigma_A_x[[i]],sigma_M_x[[i]],NULL)
  print('K calculated')
  print(object.size(K_xBMs))

  f_xBMs[[i]]<-compute_f_BMs(K_xBMs,sigma_e_x[[i]],data_xGPs[[i]]$Y,NULL)
  print('f calculated')
  
  predictions_xBMs[[i]]<-matrix(NA,length(Times_test[[i]]),length(sigma_e_x[[i]]))
  times_M<-cbind(Times_test[[i]],rep(1:M, each = length(Times_test[[i]])/M))
  for(j in 1:dim(times_M)[1]){
    predictions_xBMs[[i]][j,]<-compute_Y_pred_xBM(times_M[j,1],times_M[j,2],data_xGPs[[i]]$dim_Y,data_xGPs[[i]]$Times,f_xBMs[[i]],K_xBMs,data_xGPs[[i]]$Y,sigma_A_x[[i]], 
                                                  sigma_M_x[[i]],sigma_e_x[[i]],fun=NULL)
  }
  print('predictions done')
  rm(K_xBMs)
  generated_quantites_xBMs<-list(f_xBMs)
  save(generated_quantites_xBMs,file=paste0("...//preq_analysis_generated_quantities_xBMs_",iterations,"iter.RData"))
  
  save(predictions_xBMs,file=paste0("...//preq_analysis_predictions_xBMs_",iterations,"iter.RData"))
  
}
end_time<-Sys.time()
print(end_time-start_time)

rm(sigma_A_x)
rm(sigma_M_x)
rm(sigma_e_x)



#iBMs predictions and generated quantities
sigma<-list()

sigma_e<-list()
f_iBM<-list()
predictions_iBM<-list()

start_time<-Sys.time()
for(i in 1:N_preq){
  print(i)
  sigma[[i]]<-list()
  sigma_e[[i]]<-list()
  predictions_iBM[[i]]<-list()
  matrix_of_draws1 <- as.data.frame(fit_iBM[[i]][['GR']])
  matrix_of_draws2 <- as.data.frame(fit_iBM[[i]][['UK']])
  matrix_of_draws3<- as.data.frame(fit_iBM[[i]][['GER']])
  sigma[[i]][['GR']] <-matrix_of_draws1[,'sigma']
  sigma_e[[i]][['GR']] <-matrix_of_draws1[,'sigma_e']
  sigma[[i]][['UK']] <-matrix_of_draws2[,'sigma']
  sigma_e[[i]][['UK']] <-matrix_of_draws2[,'sigma_e']
  sigma[[i]][['GER']] <-matrix_of_draws3[,'sigma']
  sigma_e[[i]][['GER']] <-matrix_of_draws3[,'sigma_e']
  rm(matrix_of_draws1)
  rm(matrix_of_draws2)
  rm(matrix_of_draws3)
  K_iBM_GR<-compute_K_iBM(data_iGP[[i]][['GR']]$Times,sigma[[i]][['GR']],NULL)
  K_iBM_UK<-compute_K_iBM(data_iGP[[i]][['UK']]$Times,sigma[[i]][['UK']],NULL)
  K_iBM_GER<-compute_K_iBM(data_iGP[[i]][['GER']]$Times,sigma[[i]][['GER']],NULL)
  
  print('Ks calculated')
  #print(object.size(K_xBMs))
  
  f_iBM[[i]]<-list()
  start_time3<-Sys.time()
  f_iBM[[i]][['GR']]<-compute_f_iBM(K_iBM_GR,sigma_e[[i]][['GR']],data_iGP[[i]][['GR']]$Y,NULL)
  f_iBM[[i]][['UK']]<-compute_f_iBM(K_iBM_UK,sigma_e[[i]][['UK']],data_iGP[[i]][['UK']]$Y,NULL)
  f_iBM[[i]][['GER']]<-compute_f_iBM(K_iBM_GER,sigma_e[[i]][['GER']],data_iGP[[i]][['GER']]$Y,NULL)
  
  print('fs calculated')
  
  generated_quantites_iBM<-list(f_iBM)
  save(generated_quantites_iBM,file=paste0("...//preq_analysis_generated_quantities_iBM_",iterations,"iter.RData"))
  
  predictions_iBM[[i]][['GR']]<-matrix(NA,length(Times_test[[i]]),length(sigma_e[[i]][['GR']]))
  predictions_iBM[[i]][['UK']]<-matrix(NA,length(Times_test[[i]]),length(sigma_e[[i]][['UK']]))
  predictions_iBM[[i]][['GER']]<-matrix(NA,length(Times_test[[i]]),length(sigma_e[[i]][['GER']]))
  #times_M<-cbind(Times_test[[i]],rep(1:M, each = length(Times_test[[i]])/M))
  for(j in 1:(length(Times_test[[i]])/M)){
    predictions_iBM[[i]][['GR']][j,]<-compute_Y_pred_iBM(Times_test[[i]][j],data_iGP[[i]][['GR']]$Times,f_iBM[[i]][['GR']],K_iBM_GR,data_iGP[[i]][['GR']]$Y,sigma[[i]][['GR']], sigma_e[[i]][['GR']],fun=NULL)
    predictions_iBM[[i]][['UK']][j,]<-compute_Y_pred_iBM(Times_test[[i]][j],data_iGP[[i]][['UK']]$Times,f_iBM[[i]][['UK']],K_iBM_UK,data_iGP[[i]][['UK']]$Y,sigma[[i]][['UK']], sigma_e[[i]][['UK']],fun=NULL)
    predictions_iBM[[i]][['GER']][j,]<-compute_Y_pred_iBM(Times_test[[i]][j],data_iGP[[i]][['GER']]$Times,f_iBM[[i]][['GER']],K_iBM_GER,data_iGP[[i]][['GER']]$Y,sigma[[i]][['GER']], sigma_e[[i]][['GER']],fun=NULL)
  }
  
  save(predictions_iBM,file=paste0("...//preq_analysis_predictions_iBM_",iterations,"iter.RData"))
  
  
  rm(K_iBM_GR)
  rm(K_iBM_UK)
  rm(K_iBM_GER)
}
end_time<-Sys.time()
print(end_time-start_time)

#xBMs without M predictions and generated quantities
sigma_A_x_without_M<-list()
sigma_e_x_without_M<-list()
f_xBMs_without_M<-list()
predictions_xBMs_without_M<-list()
# f_M_xBMs<-list()

start_time<-Sys.time()
for(i in 1:N_preq){
  print('xBMs')
  print(i)
  
  matrix_of_draws_x_without_M <- as.data.frame(fit_xBMs_without_M[[i]])
  
  sigma_A_x_without_M[[i]]     <-matrix_of_draws_x_without_M[,'sigma_A']   
  sigma_e_x_without_M[[i]]   <-matrix_of_draws_x_without_M[,'sigma_e']
  rm(matrix_of_draws_without_M)
  K_xBMs_without_M<-compute_K_xBM_without_M(data_xGPs[[i]]$dim_Y,data_xGPs[[i]]$Times,sigma_A_x_without_M[[i]],NULL)
  print('K calculated')
  print(object.size(K_xBMs_without_M))
  
  f_xBMs_without_M[[i]]<-compute_f_BMs(K_xBMs_without_M,sigma_e_x_without_M[[i]],data_xGPs[[i]]$Y,NULL)
  print('f calculated')
  
  predictions_xBMs_without_M[[i]]<-matrix(NA,length(Times_test[[i]]),length(sigma_e_x_without_M[[i]]))
  times_M<-cbind(Times_test[[i]],rep(1:M, each = length(Times_test[[i]])/M))
  for(j in 1:dim(times_M)[1]){
    predictions_xBMs_without_M[[i]][j,]<-compute_Y_pred_xBM_without_M(times_M[j,1],times_M[j,2],data_xGPs[[i]]$dim_Y,data_xGPs[[i]]$Times,f_xBMs_without_M[[i]],K_xBMs_without_M,data_xGPs[[i]]$Y,sigma_A_x_without_M[[i]], 
                                                                      sigma_e_x_without_M[[i]],fun=NULL)
  }
  print('predictions done')
  rm(K_xBMs_without_M)
  generated_quantites_xBMs_without_M<-list(f_xBMs_without_M)
  save(generated_quantites_xBMs_without_M,file=paste0("...//preq_analysis_generated_quantities_xBMs_without_M_",iterations,"iter.RData"))
  
  save(predictions_xBMs_without_M,file=paste0("...//preq_analysis_predictions_xBMs_without_M_",iterations,"iter.RData"))
  
}
end_time<-Sys.time()
print(end_time-start_time)

rm(sigma_A_x_without_M)
rm(sigma_e_x_without_M)

