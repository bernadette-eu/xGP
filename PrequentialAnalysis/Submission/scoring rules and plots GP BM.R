lib <- c("MASS",
         "dplyr",
         "lubridate",
         "ggplot2",
         "tidyr",
         "rstan",
         "loo",
         "bayesplot",
         "scoringRules", # for CRPS and log score
         "Metrics",       # for RMSE and MAE
         "xtable"
         )
lapply(lib, require, character.only = TRUE)

output_path   <- "...//Outputs//"
plot_filepath <- output_path

#load BM
load(paste0(output_path, "preq_analysis_predictions_mxBMs_20000iter.RData"))
load(paste0(output_path, "preq_analysis_predictions_iBM_20000iter.RData"))
load(paste0(output_path, "preq_analysis_predictions_xBMs_20000iter.RData"))
#load GP
load(paste0(output_path, "preq_analysis_predictions_xGPs_20000iter.RData"))
load(paste0(output_path, "preq_analysis_predictions_mxGPs_20000iter.RData"))
load(paste0(output_path, "preq_analysis_predictions_iGP_20000iter.RData"))

#---- Set system locale to English:
Sys.setlocale("LC_ALL", "English")

df  <- read.csv("...//1_Data//rt1.csv")
df2 <- read.csv("...//1_Data//rt2.csv")

df$date<-as.Date(df$date)

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
LUX <-df[df$country=='Luxembourg',]
LITH<-df[df$country=='Lithuania',]


min_date <- min(c(GR$date,UK$date,POL$date,DEN$date))
max_date <- max(c(GR$date,UK$date,POL$date,DEN$date))

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

data_xGPs <- list()
data_iGP  <- list()
for(i in 1:N_preq){
  dates_train[[i]]  <- seq(from = start_date, to = start_date+horizon_of_study+7*(i-1), by = "day")
  dates_test[[i]]   <- seq(from = start_date+horizon_of_study+7*(i-1)+1, to = start_date+horizon_of_study+7*(i), by = "day")
  logRt_train[[i]]  <- c(log(GR$median[which(GR$date %in% dates_train[[i]])]),log(UK$median[which(UK$date %in% dates_train[[i]])]),log(GER$median[which(GER$date %in% dates_train[[i]])]))
  logRt_test[[i]]   <- c(log(GR$median[which(GR$date %in% dates_test[[i]])]),log(UK$median[which(UK$date %in% dates_test[[i]])]),log(GER$median[which(GER$date %in% dates_test[[i]])]))
  Times_train[[i]]  <- rep(1:length(log(GR$median[which(GR$date %in% dates_train[[i]])])),M)
  Times_test[[i]]   <- rep((length(log(GR$median[which(GR$date %in% dates_train[[i]])]))+1):(length(log(GR$median[which(GR$date %in% dates_train[[i]])]))+7),M)
  dim_Y[[i]]        <- c(length(log(GR$median[which(GR$date %in% dates_train[[i]])])),length(log(UK$median[which(UK$date %in% dates_train[[i]])])),length(log(GER$median[which(GER$date %in% dates_train[[i]])])))
  
  print(paste0(i,' Length logRt: ',length(logRt_train[[i]]),' sum(dim_Y): ',sum(dim_Y[[i]])))
 
  data_xGPs[[i]]         <- list(M=M,N=length(logRt_train[[i]]),dim_Y=dim_Y[[i]],Y=logRt_train[[i]],Times=Times_train[[i]])
  data_iGP[[i]]          <- list()
  data_iGP[[i]][['GR']]  <- list(N=length(log(GR$median[which(GR$date %in% dates_train[[i]])])),Y=log(GR$median[which(GR$date %in% dates_train[[i]])]),Times=1:length(log(GR$median[which(GR$date %in% dates_train[[i]])])))
  data_iGP[[i]][['UK']]  <- list(N=length(log(GR$median[which(GR$date %in% dates_train[[i]])])),Y=log(UK$median[which(UK$date %in% dates_train[[i]])]),Times=1:length(log(UK$median[which(UK$date %in% dates_train[[i]])])))
  data_iGP[[i]][['GER']] <- list(N=length(log(GER$median[which(GER$date %in% dates_train[[i]])])),Y=log(GER$median[which(GER$date %in% dates_train[[i]])]),Times=1:length(log(GER$median[which(GR$date %in% dates_train[[i]])])))
}

##################################################################
# 
# SCORING RULES
#
##################################################################
evaluate_predictive_scores <- function(y_true, y_pred_samples) {
  # Load necessary packages
  if (!requireNamespace("scoringRules", quietly = TRUE)) {
    stop("Please install the 'scoringRules' package.")
  }
  if (!requireNamespace("Metrics", quietly = TRUE)) {
    stop("Please install the 'Metrics' package.")
  }
  
  # Input checks
  if (length(y_true) != nrow(y_pred_samples)) {
    stop("Length of y_true must match number of rows in y_pred_samples.")
  }
  
  # Compute scores
  crps_vals      <- scoringRules::crps_sample(y = y_true, dat = y_pred_samples)
  log_score_vals <- scoringRules::logs_sample(y = y_true, dat = y_pred_samples)
  
  y_pred_mean <- rowMeans(y_pred_samples)
  rmse_val    <- Metrics::rmse(y_true, y_pred_mean)
  mae_val     <- Metrics::mae(y_true, y_pred_mean)
  
  # Return as named list
  return(list(
    crps           = crps_vals,
    log_score      = log_score_vals,
    rmse           = rmse_val,
    mae            = mae_val,
    mean_crps      = mean(crps_vals),
    mean_log_score = mean(log_score_vals)
  ))
}

evaluate_predictive_scores_by_group <- function(y_true, y_pred_samples, group_labels) {
  # Load required packages
  if (!requireNamespace("scoringRules", quietly = TRUE)) stop("Please install 'scoringRules'")
  if (!requireNamespace("Metrics", quietly = TRUE)) stop("Please install 'Metrics'")
  
  # Input checks
  N <- length(y_true)
  if (nrow(y_pred_samples) != N) stop("y_true and y_pred_samples must have the same number of rows.")
  if (N %% length(group_labels) != 0) stop("y_true must be divisible by the number of groups.")
  
  # Determine group structure
  group_size <- N / length(group_labels)
  group_names <- unique(group_labels)
  
  # Store results
  results <- list()
  
  for (i in seq_along(group_labels)) {
    group <- group_labels[i]
    idx <- ((i - 1) * group_size + 1):(i * group_size)
    
    y_sub <- y_true[idx]
    pred_sub <- y_pred_samples[idx, , drop = FALSE]
    
    crps_vals <- scoringRules::crps_sample(y = y_sub, dat = pred_sub)
    log_score_vals <- scoringRules::logs_sample(y = y_sub, dat = pred_sub)
    y_mean <- rowMeans(pred_sub)
    
    results[[group]] <- list(
      mean_crps = mean(crps_vals),
      mean_log_score = mean(log_score_vals),
      rmse = Metrics::rmse(y_sub, y_mean),
      mae = Metrics::mae(y_sub, y_mean)
    )
  }
  
  # Total (overall) scores
  total_crps <- scoringRules::crps_sample(y = y_true, dat = y_pred_samples)
  total_log <- scoringRules::logs_sample(y = y_true, dat = y_pred_samples)
  total_mean <- rowMeans(y_pred_samples)
  
  results[["Total"]] <- list(
    mean_crps = mean(total_crps),
    mean_log_score = mean(total_log),
    rmse = Metrics::rmse(y_true, total_mean),
    mae = Metrics::mae(y_true, total_mean)
  )
  
  return(results)
}

evaluate_predictive_scores_multiple_times <- function(y_true, y_pred_samples, group_labels, n_times) {
  
  if (!requireNamespace("scoringRules", quietly = TRUE)) stop("Please install 'scoringRules'")
  if (!requireNamespace("Metrics", quietly = TRUE)) stop("Please install 'Metrics'")
  
  # Input checks
  total_N <- length(y_true)
  n_groups <- length(group_labels)
  if (nrow(y_pred_samples) != total_N) stop("y_true and y_pred_samples must have the same number of rows.")
  if (total_N %% (n_times * n_groups) != 0) stop("Mismatch: data length must be divisible by n_times * n_groups.")
  
  group_size <- total_N / (n_times * n_groups)
  
  results <- list()
  
  for (t in 1:n_times) {
    time_result <- list()
    time_start <- (t - 1) * group_size * n_groups + 1
    for (g in seq_along(group_labels)) {
      group <- group_labels[g]
      idx_start <- time_start + (g - 1) * group_size
      idx_end <- idx_start + group_size - 1
      
      y_sub <- y_true[idx_start:idx_end]
      pred_sub <- y_pred_samples[idx_start:idx_end, , drop = FALSE]
      
      crps_vals <- scoringRules::crps_sample(y = y_sub, dat = pred_sub)
      log_score_vals <- scoringRules::logs_sample(y = y_sub, dat = pred_sub)
      y_mean <- rowMeans(pred_sub)
      
      time_result[[group]] <- list(
        mean_crps = mean(crps_vals),
        mean_log_score = mean(log_score_vals),
        rmse = Metrics::rmse(y_sub, y_mean),
        mae = Metrics::mae(y_sub, y_mean)
      )
    }
    
    # Overall scores for this time point
    time_idx <- time_start:(time_start + group_size * n_groups - 1)
    y_time <- y_true[time_idx]
    pred_time <- y_pred_samples[time_idx, , drop = FALSE]
    
    crps_vals <- scoringRules::crps_sample(y = y_time, dat = pred_time)
    log_score_vals <- scoringRules::logs_sample(y = y_time, dat = pred_time)
    y_mean <- rowMeans(pred_time)
    
    time_result[["Total"]] <- list(
      mean_crps = mean(crps_vals),
      mean_log_score = mean(log_score_vals),
      rmse = Metrics::rmse(y_time, y_mean),
      mae = Metrics::mae(y_time, y_mean)
    )
    
    results[[paste0("Time_", t)]] <- time_result
  }
  
  # Overall across all times
  y_all <- y_true
  pred_all <- y_pred_samples
  crps_vals <- scoringRules::crps_sample(y = y_all, dat = pred_all)
  log_score_vals <- scoringRules::logs_sample(y = y_all, dat = pred_all)
  y_mean <- rowMeans(pred_all)
  
  results[["Overall"]] <- list(
    mean_crps = mean(crps_vals),
    mean_log_score = mean(log_score_vals),
    rmse = Metrics::rmse(y_all, y_mean),
    mae = Metrics::mae(y_all, y_mean)
  )
  
  return(results)
}

evaluate_predictive_scores_multiple_times_df <- function(y_true, y_pred_samples, group_labels, n_times) {
  if (!requireNamespace("scoringRules", quietly = TRUE)) stop("Please install 'scoringRules'")
  if (!requireNamespace("Metrics", quietly = TRUE)) stop("Please install 'Metrics'")
  
  total_N <- length(y_true)
  n_groups <- length(group_labels)
  if (nrow(y_pred_samples) != total_N) stop("y_true and y_pred_samples must have the same number of rows.")
  if (total_N %% (n_times * n_groups) != 0) stop("Mismatch: data length must be divisible by n_times * n_groups.")
  
  group_size <- total_N / (n_times * n_groups)
  
  results_df <- data.frame()
  
  for (t in 1:n_times) {
    time_start <- (t - 1) * group_size * n_groups + 1
    for (g in seq_along(group_labels)) {
      group <- group_labels[g]
      idx_start <- time_start + (g - 1) * group_size
      idx_end <- idx_start + group_size - 1
      
      y_sub <- y_true[idx_start:idx_end]
      pred_sub <- y_pred_samples[idx_start:idx_end, , drop = FALSE]
      
      crps_vals <- scoringRules::crps_sample(y = y_sub, dat = pred_sub)
      log_score_vals <- scoringRules::logs_sample(y = y_sub, dat = pred_sub)
      y_mean <- rowMeans(pred_sub)
      
      results_df <- rbind(results_df, data.frame(
        time = paste0("Week ", t),
        group = group,
        mean_crps = mean(crps_vals),
        mean_log_score = mean(log_score_vals),
        rmse = Metrics::rmse(y_sub, y_mean),
        mae = Metrics::mae(y_sub, y_mean)
      ))
    }
    
    # Time-level total
    time_idx <- time_start:(time_start + group_size * n_groups - 1)
    y_time <- y_true[time_idx]
    pred_time <- y_pred_samples[time_idx, , drop = FALSE]
    
    crps_vals <- scoringRules::crps_sample(y = y_time, dat = pred_time)
    log_score_vals <- scoringRules::logs_sample(y = y_time, dat = pred_time)
    y_mean <- rowMeans(pred_time)
    
    results_df <- rbind(results_df, data.frame(
      time = paste0("Week ", t),
      group = "Total",
      mean_crps = mean(crps_vals),
      mean_log_score = mean(log_score_vals),
      rmse = Metrics::rmse(y_time, y_mean),
      mae = Metrics::mae(y_time, y_mean)
    ))
  }
  
  # Overall across all times
  crps_vals <- scoringRules::crps_sample(y = y_true, dat = y_pred_samples)
  log_score_vals <- scoringRules::logs_sample(y = y_true, dat = y_pred_samples)
  y_mean <- rowMeans(y_pred_samples)
  
  results_df <- rbind(results_df, data.frame(
    time = "Overall",
    group = "Total",
    mean_crps = mean(crps_vals),
    mean_log_score = mean(log_score_vals),
    rmse = Metrics::rmse(y_true, y_mean),
    mae = Metrics::mae(y_true, y_mean)
  ))
  
  return(results_df)
}

predictions_iGP_conc<-list()
for(i in 1:N_preq){
  predictions_iGP_conc[[i]]<-rbind(na.omit(predictions_iGP[[i]][['GR']]),na.omit(predictions_iGP[[i]][['UK']]),na.omit(predictions_iGP[[i]][['GER']]))
}


predictions_iBM_conc<-list()
for(i in 1:N_preq){
  predictions_iBM_conc[[i]]<-rbind(na.omit(predictions_iBM[[i]][['GR']]),na.omit(predictions_iBM[[i]][['UK']]),na.omit(predictions_iBM[[i]][['GER']]))
}

scores_all_iBM  <- evaluate_predictive_scores_multiple_times_df(unlist(logRt_test),do.call(rbind, predictions_iBM_conc), group_labels = c('GR','UK','GER'),8) %>% group_by(group)
scores_all_xBM  <- evaluate_predictive_scores_multiple_times_df(unlist(logRt_test),do.call(rbind, predictions_xBMs),     group_labels = c('GR','UK','GER'),8) %>% group_by(group)
scores_all_mxBM <- evaluate_predictive_scores_multiple_times_df(unlist(logRt_test),do.call(rbind, predictions_mxBMs),    group_labels = c('GR','UK','GER'),8) %>% group_by(group)

scores_all_iGP  <- evaluate_predictive_scores_multiple_times_df(unlist(logRt_test),
                                                                do.call(rbind, predictions_iGP_conc), 
                                                                group_labels = c('GR','UK','GER'),8) %>% 
                   group_by(group)

scores_all_xGP <- evaluate_predictive_scores_multiple_times_df(unlist(logRt_test),
                                                               do.call(rbind, predictions_xGPs),
                                                               group_labels = c('GR','UK','GER'),8) %>% 
                  group_by(group)

scores_all_mxGP<- evaluate_predictive_scores_multiple_times_df(unlist(logRt_test),
                                                               do.call(rbind, predictions_mxGPs),
                                                               group_labels = c('GR','UK','GER'),8) %>% 
                  group_by(group)

crps_table <- list()
log_table  <- list()
rmse_table <- list()
mae_table  <- list()

for(group_idx in c('GR','UK','GER','Total')){
  crps_table[[group_idx]] <- cbind(scores_all_mxGP[scores_all_mxGP$group==group_idx,c(1,3)],scores_all_xGP[scores_all_xGP$group==group_idx,3],scores_all_iGP[scores_all_iGP$group==group_idx,3],scores_all_xBM[scores_all_xBM$group==group_idx,3],scores_all_iBM[scores_all_iBM$group==group_idx,3],scores_all_iBM[scores_all_iBM$group==group_idx,3])
  colnames(crps_table[[group_idx]]) <- c('time','mxGP','xGP','iGP','mxBM','xBM','iBM')
  log_table[[group_idx]]  <- cbind(scores_all_mxGP[scores_all_mxGP$group==group_idx,c(1,4)],scores_all_xGP[scores_all_xGP$group==group_idx,4],scores_all_iGP[scores_all_iGP$group==group_idx,4],scores_all_xBM[scores_all_xBM$group==group_idx,4],scores_all_iBM[scores_all_iBM$group==group_idx,4],scores_all_iBM[scores_all_iBM$group==group_idx,4])
  colnames(log_table[[group_idx]]) <- c('time','mxGP','xGP','iGP','mxBM','xBM','iBM')
  rmse_table[[group_idx]] <- cbind(scores_all_mxGP[scores_all_mxGP$group==group_idx,c(1,5)],scores_all_xGP[scores_all_xGP$group==group_idx,5],scores_all_iGP[scores_all_iGP$group==group_idx,5],scores_all_xBM[scores_all_xBM$group==group_idx,5],scores_all_iBM[scores_all_iBM$group==group_idx,5],scores_all_iBM[scores_all_iBM$group==group_idx,5])
  colnames(rmse_table[[group_idx]]) <- c('time','mxGP','xGP','iGP','mxBM','xBM','iBM')
  mae_table[[group_idx]]  <- cbind(scores_all_mxGP[scores_all_mxGP$group==group_idx,c(1,6)],scores_all_xGP[scores_all_xGP$group==group_idx,6],scores_all_iGP[scores_all_iGP$group==group_idx,6],scores_all_xBM[scores_all_xBM$group==group_idx,6],scores_all_iBM[scores_all_iBM$group==group_idx,6],scores_all_iBM[scores_all_iBM$group==group_idx,6])
  colnames(mae_table[[group_idx]]) <- c('time','mxGP','xGP','iGP','mxBM','xBM','iBM')
}

cols_to_check <- c('mxGP','xGP','iGP','mxBM','xBM','iBM')

#-----
crps_table_all         <- do.call(rbind, crps_table)%>%tibble::rownames_to_column( var = "Country")
crps_table_all$Country <- rep(c('GR','UK','GER','Total'),c(8,8,8,9))
crps_table_all         <- crps_table_all%>% arrange(time)
crps_table_all$best_model <- apply(
  crps_table_all[cols_to_check], 1,
  function(row) names(row)[which.min(row)]
)

crps_table_all <- crps_table_all %>%
                  mutate(across(where(is.numeric), ~ round(., 5)))

colnames(crps_table_all)[2] <- "Week"
ordered_names <- c("Week", "Country", "iBM", "xBM", "mxBM", "iGP", "xGP", "mxGP", "best_model")

crps_table_all <- crps_table_all[ordered_names]

#-----
log_table_all <- do.call(rbind, log_table)%>%tibble::rownames_to_column( var = "Country")
log_table_all$Country<-rep(c('GR','UK','GER','Total'),c(8,8,8,9))
log_table_all <- log_table_all%>% arrange(time)
log_table_all$best_model <- apply(
  log_table_all[cols_to_check], 1,
  function(row) names(row)[which.min(row)]
)

log_table_all <- log_table_all %>%
                mutate(across(where(is.numeric), ~ round(., 3)))

colnames(log_table_all)[2] <- "Week"
log_table_all <- log_table_all[ordered_names]

#-----
rmse_table_all<-do.call(rbind, rmse_table)%>%tibble::rownames_to_column( var = "Country")
rmse_table_all$Country<-rep(c('GR','UK','GER','Total'),c(8,8,8,9))
rmse_table_all<-rmse_table_all%>% arrange(time)
rmse_table_all$best_model <- apply(
  rmse_table_all[cols_to_check], 1,
  function(row) names(row)[which.min(row)]
)

rmse_table_all <- rmse_table_all %>%
  mutate(across(where(is.numeric), ~ round(., 5)))

colnames(rmse_table_all)[2] <- "Week"
rmse_table_all <- rmse_table_all[ordered_names]

#-----
mae_table_all<-do.call(rbind, mae_table)%>%tibble::rownames_to_column( var = "Country")
mae_table_all$Country<-rep(c('GR','UK','GER','Total'),c(8,8,8,9))
mae_table_all<-mae_table_all%>% arrange(time)
mae_table_all$best_model <- apply(
  mae_table_all[cols_to_check], 1,
  function(row) names(row)[which.min(row)]
)
mae_table_all <- mae_table_all %>%
                 mutate(across(where(is.numeric), ~ round(., 5)))

colnames(mae_table_all)[2] <- "Week"
mae_table_all <- mae_table_all[ordered_names]

#---- Export tables in latex format:
print(xtable(crps_table_all, digits = 3), include.rownames=FALSE)
print(xtable(log_table_all,  digits = 3), include.rownames=FALSE)
print(xtable(rmse_table_all, digits = 3), include.rownames=FALSE)
print(xtable(mae_table_all,  digits = 3), include.rownames=FALSE)

# Main text table:
predictive_ability_overall_per_criterion <- 
  rbind(crps_table_all[1,],
        log_table_all[1,],
        rmse_table_all[1,],
        mae_table_all[1,])
colnames(predictive_ability_overall_per_criterion)[1] <- "Criterion"
predictive_ability_overall_per_criterion$Country      <- NULL
predictive_ability_overall_per_criterion[,1] <- c("CRPS", "Log score", "RMSE", "MAE") 

print(xtable(predictive_ability_overall_per_criterion,  digits = 3), include.rownames = FALSE)

##################################################################
# 
# Graphical outputs
#
##################################################################
Sys.setlocale("LC_ALL", "English")

plot_Rt_by_country_compare_models_all_together <- function(logRt_test, 
                                                           predictions_xGPs, 
                                                           predictions_mxGPs, 
                                                           predictions_iGPs,
                                                           predictions_xBMs, 
                                                           predictions_mxBMs, 
                                                           predictions_iBMs, 
                                                           dates_test,
                                                           ncol = 3) {
  countries   <- c("Greece", "United Kingdom", "Germany")
  n_weeks     <- length(logRt_test)
  n_countries <- length(countries)
  timepoints_per_country <- 7
  
  model_list <- list(
    xEQ   = lapply(predictions_xGPs, exp),
    mxEQ  = lapply(predictions_mxGPs,exp),
    iEQ   = lapply(predictions_iGPs, exp),
    xBM  = lapply(predictions_xBMs, exp),
    mxBM = lapply(predictions_mxBMs,exp),
    iBM  = lapply(predictions_iBMs, exp)
  )
  Rt_test <- lapply(logRt_test,exp)
  
  all_data <- list()
  
  for (model_name in names(model_list)) {
    model_preds <- model_list[[model_name]]
    
    for (week in seq_len(n_weeks)) {
      obs         <- Rt_test[[week]]
      date_seq    <- as.Date(dates_test[[week]])
      pred_matrix <- matrix(model_preds[[week]], ncol = 3000, byrow = FALSE)
      
      pred_mean  <- rowMeans(pred_matrix)
      pred_lower <- apply(pred_matrix, 1, quantile, probs = 0.025)
      pred_upper <- apply(pred_matrix, 1, quantile, probs = 0.975)
      
      for (j in seq_len(n_countries)) {
        idx_start <- (j - 1) * timepoints_per_country + 1
        idx_end   <- j * timepoints_per_country
        
        df <- data.frame(
          Date     = date_seq,  # <- use full 7-day vector for each country
          Country  = countries[j],
          Observed = obs[idx_start:idx_end],
          Mean     = pred_mean[idx_start:idx_end],
          Lower    = pred_lower[idx_start:idx_end],
          Upper    = pred_upper[idx_start:idx_end],
          Model    = model_name
        )
        
        all_data[[length(all_data) + 1]] <- df
      }
    }
  }
  
  full_df <- bind_rows(all_data)
  
  # Plot all countries in one figure with facets
  p <- ggplot(full_df, 
              aes(x = Date, 
                  color = Model, 
                  fill = Model)) +
    geom_ribbon(aes(ymin = Lower, 
                    ymax = Upper), 
                alpha = 0.15, linetype = 0) +
    geom_line(aes(y = Mean), linewidth = 1) +
    geom_point(aes(y = Observed), color = "black", size = 1.5) +
    #facet_wrap(~ Country, nrow = 1) +
    facet_wrap(. ~ Country, 
               ncol           = ncol,
               strip.position = "top",
               shrink         = T) + 
    labs(x = "Epidemiological Date",
         y = "Effective reproduction number") +
    scale_x_date(date_breaks = "1 week", date_labels = "%d-%b-%y") +
    #theme_minimal(base_size = 16) +  # increase base font size
    #theme_bw(base_size = 16) +
    theme_bw() +
    theme(
      axis.title   = element_text(size = 18, face = "bold"),
      axis.text    = element_text(size = 12, face = "bold"),
      strip.text   = element_text(size = 16, face = "bold"),
      legend.title = element_text(size = 16, face = "bold"),
      legend.text  = element_text(size = 14, face = "bold"),
      plot.title   = element_text(size = 20, face = "bold", hjust = 0.5),
      axis.text.x  = element_text(angle = 45, hjust = 1),
      legend.position = "bottom"
    )
  
  print(p)
}

predictionPlot <- 
plot_Rt_by_country_compare_models_all_together(
  logRt_test,
  predictions_xGPs,
  predictions_mxGPs,
  predictions_iGP_conc,
  predictions_xBMs,
  predictions_mxBMs,
  predictions_iBM_conc,
  dates_test,
  ncol = 1
)

# Preserve transparency when saving to .eps
# Source: https://www.sthda.com/english/wiki/saving-high-resolution-ggplots-how-to-preserve-semi-transparency
grDevices::cairo_ps(filename = paste0(plot_filepath,'prequential_analysis_predicted_rt.eps'),
                    width     = 12, 
                    height    = 12, 
                    pointsize = 12,
                    fallback_resolution = 300)
print(predictionPlot)
dev.off()
