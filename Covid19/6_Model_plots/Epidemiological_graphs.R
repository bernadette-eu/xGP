###########################################################################
#
# Load libraries:
#
###########################################################################
lib <- c("readr",
         "dplyr",
         "tidyr",
         "rstan",
         "bayesplot",
         "extraDistr",
         "MASS",
         "gridExtra",
         "compiler"
)

lapply(lib, require, character.only = TRUE)

###########################################################################
#
# Setup:
#
###########################################################################

#--- Paths:
pc_path          <- "C://Users//bouranis//"
project_path     <- paste0(pc_path,"OneDrive - aueb.gr/BERNADETTE/")
experiments_path <- "10_Stan_Project_code/13_Exchange_BM"
path_to_source   <- paste0(pc_path, "OneDrive - aueb.gr//BERNADETTE//10_Stan_Project_code//11_SIR_GR_Age_model4//Github_files//")

#---- Paths to model outputs:
store_model_out_model1 <- paste0(project_path, 
                                 experiments_path, 
                                 "/3_Output/", 
                                 #"SEEIIR_renewal_multitype_XGBM_v4E2", 
                                 "SEEIIR_renewal_multitype_XGBM_v5C",
                                 ".RData")

rt_path <- paste0(project_path, 
                  experiments_path, 
                  "/3_Output/", 
                  "SEEIIR_renewal_multitype_XGBM_v5C_rt_renewal", 
                  ".RData")  

# singletype_model_output_path <- paste0(pc_path,
#                                        "OneDrive - aueb.gr//BERNADETTE//10_Stan_Project_code//Flaxman_GBM/5_Output//",  
#                                        "SIR_Flaxman_GBM_v5", 
#                                        ".RData")

#---- Functions:
`%nin%` <- Negate(`%in%`)

#-- Data engineering; daily time series of new infections and deaths for a given period:
source(paste0(path_to_source, "2_Data_Engineering", ".R"))

#-- Calculation of the effective reproduction number with the NGM method:
source(paste0(project_path, experiments_path, 
              "//2_Code//3_Multi_SEEIIR//3_Model_comparison//", 
              "7_rt_renewal", ".R"))

# mortality_data  = dt_mortality_analysis_period
# infections_data = dt_cases_analysis_period
# nuts_fit        = nuts_fit_1

posterior_age_specific_counts <- function(mortality_data, 
                                          infections_data,
                                          nuts_fit,
                                          cov_data){
  
  #---- Checks:
  if(class(nuts_fit)[1] != "stanfit") stop("Provide an object of class 'stanfit' using rstan::sampling()")
  
  #---- Posterior draws:
  posterior_draws <- rstan::extract(nuts_fit)
  
  #---- Placeholder for output tables:
  data_cases_cols      <- c("Date", "Group", "New_Cases", "median", "low", "high")
  data_cases           <- data.frame(matrix(ncol = length(data_cases_cols), nrow = 0))
  colnames(data_cases) <- data_cases_cols
  
  data_deaths_cols      <- c("Date", "Group", "New_Deaths", "median", "low", "high")
  data_deaths           <- data.frame(matrix(ncol = length(data_deaths_cols), nrow = 0))
  colnames(data_deaths) <- data_deaths_cols
  
  if( colnames(infections_data)[5] == "Total_Cases") colnames(infections_data)[5] <- "New_Cases"
  
  if( "y_deaths" %in% names(cov_data)) {
    group_id <- colnames(cov_data$y_deaths)
    
  } else { 
    group_id <- colnames(cov_data$y_data)
  
  }# End if
  
  # Loop over each age group:
  for (i in 1:cov_data$A){
    
    word <- c("E_casesByAge", "E_deathsByAge")

    if( any( word %in% names(posterior_draws)) ){
      
      fit_cases  <- posterior_draws$E_casesByAge[,,i]  #E_casesByAge[,,i]
      fit_deaths <- posterior_draws$E_deathsByAge[,,i] #E_deathsByAge[,,i]  #
      
    } else {
      fit_cases  <- posterior_draws$E_casesAge[,,i]  #E_casesByAge[,,i]
      fit_deaths <- posterior_draws$E_deathsAge[,,i] #E_deathsByAge[,,i]  #
    }

    dt_cases_age_grp  <- data.frame(Date  = unique(mortality_data$Date),
                                    Group = rep(group_id[i], length(unique(mortality_data$Date)) ))
    dt_deaths_age_grp <- dt_cases_age_grp
    
    
    filter_out_infections_cols <- c("Index", "Right", "Week_ID", "New_Cases")
    filter_out_deaths_cols     <- c("Index", "Right", "Week_ID", "New_Deaths")
    
    
    dt_cases_analysis_period_melt <- infections_data[, colnames(infections_data) %nin% 
                                                       filter_out_infections_cols] %>% 
      tidyr::gather(Group, New_Cases, -Date)
    
    dt_mortality_analysis_period_melt <- mortality_data[, colnames(mortality_data) %nin% 
                                                          filter_out_deaths_cols] %>% 
      tidyr::gather(key = "Group", value = "New_Deaths", -Date)
    
    dt_cases_age_grp <- dt_cases_age_grp %>%
      dplyr::left_join(dt_cases_analysis_period_melt,
                       by = c("Date"  = "Date",
                              "Group" = "Group"))
    
    dt_deaths_age_grp <- dt_deaths_age_grp %>%
      dplyr::left_join(dt_mortality_analysis_period_melt,
                       by = c("Date"  = "Date",
                              "Group" = "Group"))
    
    # Add quantiles from the model outputs:
    dt_cases_age_grp$median <- apply(fit_cases, 2, median)
    dt_cases_age_grp$low    <- apply(fit_cases, 2, quantile, probs = c(0.025))
    dt_cases_age_grp$high   <- apply(fit_cases, 2, quantile, probs = c(0.975))
    dt_cases_age_grp$low25  <- apply(fit_cases, 2, quantile, probs = c(0.25))
    dt_cases_age_grp$high75 <- apply(fit_cases, 2, quantile, probs = c(0.75))
    
    dt_deaths_age_grp$median <- apply(fit_deaths, 2, median)
    dt_deaths_age_grp$low    <- apply(fit_deaths, 2, quantile, probs = c(0.025))
    dt_deaths_age_grp$high   <- apply(fit_deaths, 2, quantile, probs = c(0.975))
    dt_deaths_age_grp$low25  <- apply(fit_deaths, 2, quantile, probs = c(0.25))
    dt_deaths_age_grp$high75 <- apply(fit_deaths, 2, quantile, probs = c(0.75))
    
    data_cases  <- rbind(data_cases, dt_cases_age_grp)
    data_deaths <- rbind(data_deaths, dt_deaths_age_grp)
    
  }# End for
  
  output <- list(Deaths     = data_deaths,
                 Infections = data_cases
  )
  
  return(output)
  
}# End function

#---- Posterior random draws of age-specific mortality counts:
# model_out = nuts_fit_1
# cov_data  = stan_data
# reported_data = data_rt_wide
deaths_random_draws <- function(cov_data,
                                model_out,
                                age_specific_deaths,
                                aggregated_deaths){
  
  set.seed(1)
  
  #---- Checks:
  if(class(model_out)[1] != "stanfit") stop("Provide an object of class 'stanfit' using rstan::sampling()")
  
  #---- Posterior draws:
  posterior_draws <- rstan::extract(model_out)
  E_deathsByAge   <- posterior_draws[["E_deathsAge"]]
  E_deaths        <- posterior_draws[["E_deaths"]]
  phi             <- posterior_draws[["phiD"]]
  
  mcmc_length     <- dim(E_deathsByAge)[1]
  ts_length       <- dim(E_deathsByAge)[2]
  age_grps        <- dim(E_deathsByAge)[3]
  dates           <- age_specific_deaths$Date
  
  death_draws           <- array(NA, c(mcmc_length, ts_length, age_grps))
  data_deaths_cols      <- c("Date",
                             "Group",
                             "rng_low",
                             "rng_low25",
                             "rng_high75",
                             "rng_high")
  age_rng_draws           <- data.frame(matrix(ncol = length(data_deaths_cols),
                                               nrow = 0))
  colnames(age_rng_draws) <- data_deaths_cols
  aggregated_rng_draws    <- age_rng_draws
  
  message(" > Estimation for age-specific deaths")
  
  for (k in 1:age_grps) {
    message(paste0(" > Estimation in group ", k))
    
    for (i in 1:mcmc_length) {
      for (j in 1:ts_length) {
        death_draws[i,j,k] <- rnbinom(1,
                                      mu   = E_deathsByAge[i,j,k],
                                      size = phi[i])
      }# End for
    }# End for
  }# End for
  
  for (k in 1:age_grps) {
    rng_draws_age_grp  <- data.frame(Date  = dates,
                                     Group = rep( colnames(cov_data$y_data)[k],
                                                  length(dates)
                                     ))
    rng_draws_age_grp$rng_low    <- apply(death_draws[,,k], 2, quantile, probs = c(0.025))
    rng_draws_age_grp$rng_low25  <- apply(death_draws[,,k], 2, quantile, probs = c(0.25))
    rng_draws_age_grp$rng_high75 <- apply(death_draws[,,k], 2, quantile, probs = c(0.75))
    rng_draws_age_grp$rng_high   <- apply(death_draws[,,k], 2, quantile, probs = c(0.975))
    
    age_rng_draws <- rbind(age_rng_draws, rng_draws_age_grp)
  }# End for
  
  message(" > Estimation for aggregated deaths")
  
  aggregated_death_draws <- apply(death_draws, MARGIN = c(1,2), sum)
  aggregated_rng_draws   <- data.frame(Date  = dates)
  
  aggregated_rng_draws$rng_low    <- apply(aggregated_death_draws, 2, quantile, probs = c(0.025)) 
  aggregated_rng_draws$rng_low25  <- apply(aggregated_death_draws, 2, quantile, probs = c(0.25))  
  aggregated_rng_draws$rng_high75 <- apply(aggregated_death_draws, 2, quantile, probs = c(0.75))  
  aggregated_rng_draws$rng_high   <- apply(aggregated_death_draws, 2, quantile, probs = c(0.975))
  
  ###################################
  age_rng_draws        <- age_rng_draws %>% distinct(Date, Group, .keep_all = TRUE)
  aggregated_rng_draws <- aggregated_rng_draws %>% distinct(Date, .keep_all = TRUE)
  
  output <- list(age_specific_deaths = age_rng_draws,
                 aggregated_deaths   = aggregated_rng_draws)
  
  return(output)
  
}# End function
deaths_random_draws <- compiler::cmpfun(deaths_random_draws)

###########################################################################
#
# Data management:
#
###########################################################################

#---- Country and analysis period:
country    <- "England"
start_date <- "2020-03-02"
end_date   <- "2020-09-27"

#---- Age mappings:
# The age groups with mortality data are: "0-19"  "20-39" "40-59" "60-79" "80+"
# The NHS groups are "0-19"  "20-39" "40-59" "60-79" "80+":
age_mapping_deaths <- c("0-19",
                        "20-39",
                        "40-59",
                        rep("60+",2) )

# The NHS groups are:
# [1] "0 to 4"      "5 to 9"      "10 to 14"    "15 to 19"    "20 to 24"    "25 to 29"    "30 to 34"   
# [8] "35 to 39"    "40 to 44"    "45 to 49"    "50 to 54"    "55 to 59"    "60 to 64"    "65 to 69"   
# [15] "70 to 74"    "75 to 79"    "80 to 84"    "85 to 89"    "90 and over"
age_mapping_cases <- c(rep("0-19",  4),
                       rep("20-39", 4),
                       rep("40-59", 4),
                       rep("60+",   7))

age_mapping_contacts <- c(rep("0-19",  4),
                          rep("20-39", 4),
                          rep("40-59", 4),
                          rep("60+",   4))

age_mapping_ifr_react <- age_mapping_cases

data_all <- data_engineering_observations_ENG(start_date         = start_date,
                                              end_date           = end_date,
                                              age_mapping_deaths = age_mapping_deaths,
                                              age_mapping_cases  = age_mapping_cases)

dt_mortality_analysis_period <- data_all$Deaths
dt_cases_analysis_period     <- data_all$Cases
aggr_age                     <- data_all$Age_distribution_aggr

###########################################################################
#
# Population infection rate (single-type model)
#
###########################################################################
# load(file = singletype_model_output_path)
# 
# singletype_posts_1 <- rstan::extract(nuts_fit_1)
# singletype_posts_1 <- singletype_posts_1$Rt_eff
# dI                 <- 7 # days
# 
# transm_rate_singletype_model <- singletype_posts_1 / matrix(data = dI, nrow = nrow(singletype_posts_1), ncol = ncol(singletype_posts_1) )
# 
# # Checks of the priors used in Monod et al (2021):
# # test <- rnorm(1, -0.6833129, 0.242431)
# # exp(test)
# # 
# # test2 <- rnorm(1000, 0.3828, 0.1638)
# # exp(test2)
# # 
# # test3 <- rnorm(1, -1.0702, 0.2170)
# # exp(test3)
# 
# median_fit_transmrate <- apply(transm_rate_singletype_model, 2, median)
# low_fit_transmrate    <- apply(transm_rate_singletype_model, 2, quantile, probs = c(0.025)) 
# high_fit_transmrate   <- apply(transm_rate_singletype_model, 2, quantile, probs = c(0.975))
# low25_fit_transmrate  <- apply(transm_rate_singletype_model, 2, quantile, probs = c(0.25))
# high75_fit_transmrate <- apply(transm_rate_singletype_model, 2, quantile, probs = c(0.75))
# 
# transmrate_plot_data        <- dt_mortality_analysis_period[,3, drop = FALSE]
# transmrate_plot_data$median <- median_fit_transmrate
# transmrate_plot_data$low    <- low_fit_transmrate
# transmrate_plot_data$high   <- high_fit_transmrate
# transmrate_plot_data$low25  <- low25_fit_transmrate
# transmrate_plot_data$high75 <- high75_fit_transmrate
# 
# transmrate_plot_data <- transmrate_plot_data[seq(1, nrow(transmrate_plot_data), cov_data$gbm_no_days),]
# transmrate_plot_data$Group <- rep("Population", nrow(transmrate_plot_data))
# transmrate_plot_data <- transmrate_plot_data[c("Date", "Group", "median", "low", "high", "low25", "high75")]
# 
# rm(nuts_fit_1, cov_data)
# 
# #---- Plot of posterior estimates:
# ggplot(transmrate_plot_data,
#        aes(Date   = Date,
#            Median = median)) +
#   geom_line(aes(x     = Date,
#                 y     = median,
#                 color = "Median"),
#             size = 1.3) +
#   geom_ribbon(aes(x    = Date,
#                   ymin = low25,
#                   ymax = high75,
#                   fill = "50% CrI"),
#               alpha = 0.5) +
#   # geom_ribbon(aes(x    = Date,
#   #                 ymin = low,
#   #                 ymax = high,
#   #                 fill = "95% CrI"),
#   #             alpha = 0.5) +
#   # labs(x = "Epidemiological Date",
#   #      y = "Effective reproduction number") +
#   labs(x = "Epidemiological Date",
#        y = "Infection rate (single-type model)") +
#   scale_x_date(date_breaks       = "2 weeks",
#                date_minor_breaks = "2 week") +
#   #geom_hline(yintercept = 1, color = "black") +
#   # scale_y_continuous(
#   #   limits = c( min(transmrate_plot_data$low25)*0.8,   max(transmrate_plot_data$high75)*1.1),    # max(death_plot_data$high)*1.2
#   #   breaks = seq( 0.2, max(transmrate_plot_data$high75)*1.1, 0.2) #max(death_plot_data$high)*1.2
#   # ) +
#   scale_colour_manual(name   = '',
#                       values = c('Median' = "black")) +
#   scale_fill_manual(values = c("50% CrI" = "grey30"#,  #"steelblue3"
#                                #"95% CrI" = "steelblue1" #"steelblue1"
#   )
#   ) +
#   
#   theme_bw() +
#   theme(strip.placement  = "outside",
#         strip.background = element_rect(fill = NA, colour="grey50"),
#         axis.text.x      = element_text(angle = 45, hjust = 1),
#         axis.title.x     = element_text(size = 14, face = "bold"),
#         axis.title.y     = element_text(size = 14, face = "bold"),
#         panel.spacing    = unit(0,"cm"),
#         legend.position  = "bottom",
#         legend.title     = element_blank(),
#         legend.box       = "vertical",
#         legend.margin    = margin() )

###########################################################################
#
# Load the posterior draws:
#
###########################################################################
load(file = store_model_out_model1)

posts_1   <- rstan::extract(nuts_fit_1)

if( "y_deaths" %in% names(cov_data)) {
  age_bands <- colnames(cov_data$y_deaths)
} else { 
  age_bands <- colnames(cov_data$y_data)
}# End if

###########################################################################
#
# Estimation of the effective reproduction number
#
###########################################################################
time_start_Rt <- Sys.time()
post_Rt <- rt_renewal(mortality_data  = dt_mortality_analysis_period,
                      model_out       = nuts_fit_1,
                      cov_data        = cov_data,
                      progress_bar    = TRUE)
time_end_Rt <- Sys.time()
duration_Rt <- time_end_Rt - time_start_Rt

save(post_Rt,
     duration_Rt,
     file = rt_path)

# Load the estimates of the effective reproduction number:
# load(rt_path)

###########################################################################
#
# Estimated Mortality (Aggregated)
#
###########################################################################
# mortality_data  = dt_mortality_analysis_period
# infections_data = dt_cases_analysis_period
# nuts_fit        = nuts_fit_1
# cov_data        = cov_data
# group_id        = unique(age_mapping_deaths)

posterior_counts <- posterior_age_specific_counts(mortality_data  = dt_mortality_analysis_period, 
                                                  infections_data = dt_cases_analysis_period,
                                                  nuts_fit        = nuts_fit_1,
                                                  cov_data        = cov_data)

#---- MAE:
MAE <- rep(0, length(age_bands))

for (i in 1:length(age_bands)) 
  MAE[i] <- posterior_counts$Deaths %>%
  dplyr::filter(Group %in% age_bands[i]) %>% 
  dplyr::mutate(MAE_median = abs(median - New_Deaths)) %>%
  dplyr::summarise(MAE_median = mean(MAE_median))
MAE <- unlist(MAE)
names(MAE) <- age_bands
print(MAE)

###########################################################################
#
# mxGBM goodness of fit:
#
###########################################################################
mortality_data     <- data_all$Deaths
infections_data    <- data_all$Cases

age_specific_fits  <- posterior_age_specific_counts(mortality_data, 
                                                    infections_data,
                                                    nuts_fit = nuts_fit_1,
                                                    cov_data)

aggregated_fits    <- posterior_aggregated_counts(mortality_data,
                                                  infections_data,
                                                  nuts_fit = nuts_fit_1)

aggregated_deaths_data   <- aggregated_fits$Deaths
aggregated_cases_data    <- aggregated_fits$Infections
age_specific_cases_data  <- age_specific_fits$Infections
age_specific_deaths_data <- age_specific_fits$Deaths

#---- Generate random samples of age-stratified deaths from the model:
deaths_rng      <- deaths_random_draws(cov_data,
                                       model_out           = nuts_fit_1,
                                       age_specific_deaths = age_specific_deaths_data,
                                       aggregated_deaths   = aggregated_deaths_data)

age_deaths_rng  <- deaths_rng$age_specific_deaths
aggr_deaths_rng <- deaths_rng$aggregated_deaths

age_specific_deaths_data <- age_specific_deaths_data %>% left_join(age_deaths_rng,  by = c("Date", "Group"))
aggregated_deaths_data   <- aggregated_deaths_data   %>% left_join(aggr_deaths_rng, by = c("Date"))

fit_deaths_all <- 
  ggplot(aggregated_deaths_data,
         aes(Date       = Date,
             New_Deaths = New_Deaths)) +
  geom_point(aes(x   = Date, 
                 y   = New_Deaths,
                 fill= "Reported")) +
  geom_line(aes(x     = Date,
                y     = median,
                color = "Median"),
            size = 1.3) +
  geom_ribbon(aes(x    = Date,
                  ymin = rng_low,
                  ymax = rng_high,
                  fill = "95% CrI"),
              alpha = 0.5) +
  labs(x = "Epidemiological Date",
       y = "New daily deaths") +
  scale_x_date(date_breaks       = "1 month", 
               date_minor_breaks = "1 month") + 
  scale_y_continuous(
    limits = c(0,   max(aggregated_deaths_data$high)*1.4),    
    breaks = seq(0, max(aggregated_deaths_data$high)*1.4, 100) 
  ) +
  scale_fill_manual(values = c('Reported' = "black")) +
  scale_colour_manual(name = '', 
                      values = c('Median'  = "black", 
                                 "95% CrI" = "gray40"       
                      )
  ) +
  theme_bw() +
  theme(strip.placement  = "outside",
        strip.background = element_rect(fill = NA, colour="grey50"),
        axis.text.x      = element_text(angle = 45, hjust = 1),
        axis.title.x     = element_text(size = 14, face = "bold"),
        axis.title.y     = element_text(size = 14, face = "bold"),
        panel.spacing    = unit(0,"cm"),
        legend.position  = "bottom",
        legend.title     = element_blank(),
        legend.box       = "vertical", 
        legend.margin    = margin()) 

###########################################################################
#
# Goodness of fit to new age-stratified deaths:
#
###########################################################################
fit_deaths_age <- 
  ggplot(age_specific_deaths_data,
         aes(Date       = Date,
             New_Deaths = New_Deaths)) +
  facet_wrap(. ~ Group, 
             scales = "free_y", 
             ncol   = 1,
             strip.position = "right") +
  geom_point(aes(x   = Date, 
                 y   = New_Deaths,
                 fill= "Reported")) +
  geom_line(aes(x     = Date,
                y     = median,
                color = "Median"),
            size = 1.3) +
  geom_ribbon(aes(x    = Date,
                  ymin = rng_low,
                  ymax = rng_high,
                  fill = "95% CrI"),
              alpha = 0.5) +
  labs(x = "Epidemiological Date",
       y = "New daily deaths") +
  scale_x_date(date_breaks       = "1 month", 
               date_minor_breaks = "1 month") + 
  scale_fill_manual(values = c('Reported' = "black")) +
  scale_colour_manual(name = '', 
                      values = c('Median'  = "black", 
                                 "95% CrI" = "gray40"       
                      )
  ) +
  theme_bw() +
  theme(panel.spacing    = unit(0.2,"cm"),
        axis.text.x      = element_text(angle = 45, hjust = 1),
        axis.title.x     = element_text(size = 14, face = "bold"),
        axis.title.y     = element_text(size = 14, face = "bold"),
        legend.position  = "bottom",
        legend.title     = element_blank(),
        legend.box       = "vertical", 
        legend.margin    = margin())

# Size: 7x12
#---- Model fit to observed deaths:
(fit_deaths_age | fit_deaths_all) + plot_annotation(tag_levels =  "A")

#---- Posterior age-stratified transmission rate:
library(ggplot2)

main_path <- "C://Users//lbour//OneDrive//Desktop1//xgp_project//GP multitype epidemic England//4_Stan_Outputs//"
model_ids <- c("SEEIIR_renewal_multitype_XGBM")

store_model_out_model1 <- paste0(main_path, 
                                 model_ids[1],
                                 ".RData")
load(file = store_model_out_model1)
start_date <- "2020-03-02"
end_date   <- "2020-09-27"
dates     <- seq(as.Date(start_date), as.Date(end_date), by = "3 days")

transRatePosterior <- function(cov_data,
                               model_out,
                               dates){
  
  
  # model_out <- nuts_fit_1
  # data_cases <- dt_mortality_analysis_period
  # data_rt <- post_Rt
  # rm(aggr_age, model_out, nuts_fit_1)
  
  `%nin%` <- Negate(`%in%`)
  
  #---- Checks:
  if(class(model_out)[1] != "stanfit") stop("Provide an object of class 'stanfit' using rstan::sampling()")
  
  #---- Posterior draws:
  posterior_draws <- rstan::extract(model_out)
  
  #---- Checks:
  age_grps     <- cov_data$A
  beta_draws   <- posterior_draws$rho_weekly
  chain_length <- nrow(beta_draws)
  ts_length    <- dim(beta_draws)[2]
  
  #---- Initiate output table for the transmission rate
  #     per age group over time:
  dataTransmissionRateCols      <- c("Date", 
                                     "Group", 
                                     "median", 
                                     "low0025", 
                                     "low25", 
                                     "high75", 
                                     "high975")
  dataTransmissionRate           <- data.frame(matrix(ncol = length(dataTransmissionRateCols), nrow = 0))
  colnames(dataTransmissionRate) <- dataTransmissionRateCols
  
  message(" > Estimation of age-specific transmission rate")
  
  for (k in 1:cov_data$A){
    
    #message(paste(" > Age group", k))
    
    trans_rate_temp          <- matrix(0L, nrow = chain_length, ncol = ts_length)
    data_trans_rate_age_grp  <- data.frame(Date  = dates,
                                           Group = rep( colnames(cov_data$y_data)[k],
                                                        length(dates)))
    
    for (j in 1:ts_length) trans_rate_temp[,j] <- beta_draws[,j,k] * cov_data$contact_matrix[k,k]
    
    data_trans_rate_age_grp$median  <- apply(trans_rate_temp, 2, median)
    data_trans_rate_age_grp$low0025 <- apply(trans_rate_temp, 2, quantile, probs = c(0.025)) # c(0.025)
    data_trans_rate_age_grp$low25   <- apply(trans_rate_temp, 2, quantile, probs = c(0.25))  # c(0.025)
    data_trans_rate_age_grp$high75  <- apply(trans_rate_temp, 2, quantile, probs = c(0.75))  # c(0.975)
    data_trans_rate_age_grp$high975 <- apply(trans_rate_temp, 2, quantile, probs = c(0.975)) # c(0.975)
    
    dataTransmissionRate <- rbind(dataTransmissionRate,
                                  data_trans_rate_age_grp)
    
  }# End for
  
  #---- Output:
  output <- list(Transmission_rate = dataTransmissionRate)
  
  return(output)
  
}# End function
transRatePosteriorOut <- transRatePosterior(cov_data, model_out = nuts_fit_1, dates)

Sys.setlocale("LC_ALL", "English")
plotTransRate <- 
  ggplot(transRatePosteriorOut$Transmission_rate) +
  facet_wrap(. ~ Group, 
             scales = "free_y", 
             ncol   = 1,
             strip.position = "top",
             shrink = T) +
  geom_line(aes(x     = Date,
                y     = median,
                color = "Median"),
            size  = 1.3) +
  geom_ribbon(aes(x    = Date,
                  ymin = low25,
                  ymax = high75,
                  fill = "50% CrI"),
              alpha = 0.6) +
  geom_ribbon(aes(x    = Date,
                  ymin = low0025,
                  ymax = high975,
                  fill = "95% CrI"),
              alpha = 0.6) +
  labs(x = "Epidemiological date",
       y = "Transmission rate") +
  scale_x_date(date_breaks       = "1 month", 
               date_minor_breaks = "1 month") +
  scale_fill_manual(values = c("50% CrI" = "gray40",
                               "95% CrI" = "gray70"
  )) +
  scale_colour_manual(name   = '',
                      values = c('Median' = "black")) +
  theme_bw() +
  theme(panel.spacing    = unit(0.2,"cm"),
        axis.text.x      = element_text(size = 12, angle = 45, hjust = 1),
        axis.text.y      = element_text(size = 12),
        axis.title.x     = element_text(size = 16, face = "bold"),
        axis.title.y     = element_text(size = 16, face = "bold"),
        legend.position  = "bottom",
        legend.title     = element_blank(),
        #legend.box       = "vertical", 
        legend.text      = element_text(size = 14),
        legend.margin    = margin(),
        strip.text.x     = element_text(size = 16)
  )

# Preserve transparency when saving to .eps
# Source: https://www.sthda.com/english/wiki/saving-high-resolution-ggplots-how-to-preserve-semi-transparency
grDevices::cairo_ps(filename = "C://Users//lbour//OneDrive//Desktop1//xgp_project//GP multitype epidemic England//6_Model_plots//COVID19_xGBM_Transrate.eps",
                    width     = 10, 
                    height    = 10, 
                    pointsize = 12,
                    fallback_resolution = 300)
print(plotTransRate)
dev.off()

###########################################################################
#
# Posterior median new aggregated deaths
#
###########################################################################
fit_dead        <- posts_1$E_deaths
median_fit_dead <- apply(fit_dead, 2, median)
low_fit_dead    <- apply(fit_dead, 2, quantile, probs = c(0.025)) #c(0.025)
high_fit_dead   <- apply(fit_dead, 2, quantile, probs = c(0.975))#c(0.975)
low25_fit_dead  <- apply(fit_dead, 2, quantile, probs = c(0.25)) #c(0.025)
high75_fit_dead <- apply(fit_dead, 2, quantile, probs = c(0.75))#c(0.975)

death_plot_data        <- dt_mortality_analysis_period #cd
death_plot_data$median <- median_fit_dead
death_plot_data$low    <- low_fit_dead
death_plot_data$high   <- high_fit_dead
death_plot_data$low25  <- low25_fit_dead
death_plot_data$high75 <- high75_fit_dead

ggplot(death_plot_data,
       aes(Date       = Date,
           New_Deaths = New_Deaths)) +
  geom_point(aes(x   = Date, 
                 y   = New_Deaths,
                 fill= "Reported deaths")) +
  geom_line(aes(x     = Date,
                y     = median,
                color = "Median Fitted deaths"),
            size = 1.3) +
  geom_ribbon(aes(x    = Date,
                  ymin = low,
                  ymax = high,
                  fill = "95% CrI Estimated deaths"),
              alpha = 0.5) +
  geom_ribbon(aes(x    = Date,
                  ymin = low25,
                  ymax = high75,
                  fill = "50% CrI Estimated deaths"),
              alpha = 0.5) +
  labs(x = "Epidemiological Date",
       y = "New daily deaths") +
  scale_x_date(date_breaks = "2 weeks", 
               date_minor_breaks = "2 weeks") + 
  scale_y_continuous(
    limits = c(0,   1200),    # max(death_plot_data$high)*1.2
    breaks = seq(0, 1200, 200) #max(death_plot_data$high)*1.2
  ) +
  scale_fill_manual(values = c('Reported deaths' = "black")) +
  scale_colour_manual(name = '', 
                      values = c('Median Fitted deaths'    = "black", # "#0072B2"
                                 "50% CrI Estimated deaths" = "gray70",  #"steelblue3"
                                 "95% CrI Estimated deaths" = "gray40" #"steelblue1"
                      )
  ) +
  theme_bw() +
  theme(strip.placement  = "outside",
        strip.background = element_rect(fill = NA, colour="grey50"),
        axis.text.x      = element_text(angle = 45, hjust = 1),
        panel.spacing    = unit(0,"cm"),
        legend.position  = "bottom",
        legend.title     = element_blank() )

###########################################################################
#
# Estimated group-specific infections
#
###########################################################################
fit_infections_age <- 
  ggplot(posterior_counts$Infections,
         aes(Date       = Date,
             median = median)) +
  facet_wrap(. ~ Group, 
             scales = "free_y", 
             ncol   = 1,
             strip.position = "right") +
  # geom_point(aes(x   = Date, 
  #                y   = New_Cases,
  #                fill= "Reported")) +
  geom_line(aes(x     = Date,
                y     = median,
                color = "Median"),
            size = 1.3) +
  geom_ribbon(aes(x    = Date,
                  ymin = low25,
                  ymax = high75,
                  fill = "50% CrI"),
              alpha = 0.5) +
  # geom_ribbon(aes(x    = Date,
  #                 ymin = low,
  #                 ymax = high,
  #                 fill = "95% CrI median"),
  #             alpha = 0.5) +
  labs(x = "Epidemiological Date",
       y = "New daily infections") +
  scale_x_date(date_breaks       = "1 month", 
               date_minor_breaks = "1 month") + 
  #scale_fill_manual(values = c('Reported' = "black")) +
  scale_fill_manual(values = c("50% CrI" = "gray40")) +
  
  scale_colour_manual(name = '',
                      values = c('Median'         = "black"#, # "#0072B2"
                                 #"50% CrI" = "gray70"  #"steelblue3"
                                 #"95% CrI median" = "gray40", #"steelblue1"
                      )
  ) +
  theme_bw() +
  theme(panel.spacing    = unit(0.2,"cm"),
        axis.text.x      = element_text(angle = 45, hjust = 1),
        axis.title.x     = element_text(size = 14, face = "bold"),
        axis.title.y     = element_text(size = 14, face = "bold"),
        legend.position  = "bottom",
        legend.title     = element_blank(),
        legend.box       = "horizontal", 
        legend.margin    = margin())

###########################################################################
#
# Estimated new aggregated infections
#
###########################################################################
fit_cases        <- posts_1$E_cases
median_fit_cases <- apply(fit_cases, 2, median)
low_fit_cases    <- apply(fit_cases, 2, quantile, probs = c(0.025)) # c(0.025)
high_fit_cases   <- apply(fit_cases, 2, quantile, probs = c(0.975)) # c(0.975)
low25_fit_cases  <- apply(fit_cases, 2, quantile, probs = c(0.25)) #c(0.025)
high75_fit_cases <- apply(fit_cases, 2, quantile, probs = c(0.75))#c(0.975)

cases_plot_data        <- dt_cases_analysis_period#death_plot_data
cases_plot_data$median <- median_fit_cases
cases_plot_data$low    <- low_fit_cases
cases_plot_data$high   <- high_fit_cases
cases_plot_data$low25  <- low25_fit_cases
cases_plot_data$high75 <- high75_fit_cases

fit_infections_all <- 
ggplot(cases_plot_data,
       aes(Date       = Date,
           New_Cases = Total_Cases)) +
  # geom_point(aes(x   = Date, 
  #                y   = Total_Cases,
  #                fill= "Reported cases")) +
  geom_line(aes(x     = Date,
                y     = median,
                color = "Median"),
            size = 1.3) +
  geom_ribbon(aes(x    = Date,
                  ymin = low25,
                  ymax = high75,
                  fill = "50% CrI"),
              alpha = 0.5) +
  # geom_ribbon(aes(x    = Date,
  #                 ymin = low,
  #                 ymax = high,
  #                 fill = "95% CrI Estimated cases"),
  #             alpha = 0.5) +
  labs(x = "Epidemiological Date",
       y = "New daily infections") +
  scale_x_date(date_breaks = "1 month", 
               date_minor_breaks = "1 month") + 
  scale_y_continuous(
    limits = c(0, 180000)  , #max(cases_plot_data$high75)*1.08)   # max(cases_plot_data$median)*1.2
    breaks = seq(0, 180000, 30000) # max(cases_plot_data$median)*1.2
  ) +
  #scale_fill_manual(values = c('Reported cases' = "black")) +
  scale_fill_manual(values = c("50% CrI" = "gray40")) +
  scale_colour_manual(name = '', 
                      values = c('Median'  = "black"#,
                                 #"50% CrI" = "gray70"
                                 #"95% CrI Estimated cases" = "gray40",
                      )
  ) +
  theme_bw() +
  theme(panel.spacing    = unit(0.2,"cm"),
        axis.text.x      = element_text(angle = 45, hjust = 1),
        axis.title.x     = element_text(size = 14, face = "bold"),
        axis.title.y     = element_text(size = 14, face = "bold"),
        legend.position  = "bottom",
        legend.title     = element_blank(),
        legend.box       = "horizontal", 
        legend.margin    = margin())

# Size: 7x12
#---- Model fit to observed deaths:
(fit_infections_age | fit_infections_all) + plot_annotation(tag_levels =  "A")

###########################################################################
#
# Population Probability of infection trajectory (discrete renewal model equations)
#
###########################################################################
# fit_rho        <- posts_1$beta_mu
# median_fit_rho <- apply(fit_rho, 2, median)
# low_fit_rho    <- apply(fit_rho, 2, quantile, probs = c(0.025)) 
# high_fit_rho   <- apply(fit_rho, 2, quantile, probs = c(0.975))
# low25_fit_rho  <- apply(fit_rho, 2, quantile, probs = c(0.25))
# high75_fit_rho <- apply(fit_rho, 2, quantile, probs = c(0.75))
# 
# rho_plot_data        <- dt_mortality_analysis_period[,3, drop = FALSE]
# rho_plot_data$median <- median_fit_rho
# rho_plot_data$low    <- low_fit_rho
# rho_plot_data$high   <- high_fit_rho
# rho_plot_data$low25  <- low25_fit_rho
# rho_plot_data$high75 <- high75_fit_rho
# 
# #---- Plot of posterior estimates:
# ggplot(rho_plot_data,
#        aes(Date   = Date,
#            Median = median)) +
#   geom_line(aes(x     = Date,
#                 y     = median,
#                 color = "Median"),
#             size = 1.3) +
#   geom_ribbon(aes(x    = Date,
#                   ymin = low25,
#                   ymax = high75,
#                   fill = "50% CrI"),
#               alpha = 0.5) +
#   # geom_ribbon(aes(x    = Date,
#   #                 ymin = low,
#   #                 ymax = high,
#   #                 fill = "95% CrI"),
#   #             alpha = 0.5) +
#   # labs(x = "Epidemiological Date",
#   #      y = "Effective reproduction number") +
#   labs(x = "Epidemiological Date",
#        y = "Effective contact rate") +
#   scale_x_date(date_breaks       = "2 weeks",
#                date_minor_breaks = "2 week") +
#   #geom_hline(yintercept = 1, color = "black") +
#   # scale_y_continuous(
#   #   limits = c( min(rho_plot_data$low25)*0.8,   max(rho_plot_data$high75)*1.1),    # max(death_plot_data$high)*1.2
#   #   breaks = seq( 0.2, max(rho_plot_data$high75)*1.1, 0.2) #max(death_plot_data$high)*1.2
#   # ) +
#   scale_colour_manual(name   = '',
#                       values = c('Median' = "black")) +
#   scale_fill_manual(values = c("50% CrI" = "grey30"#,  #"steelblue3"
#                                #"95% CrI" = "steelblue1" #"steelblue1"
#   )
#   ) +
# 
#   theme_bw() +
#   theme(strip.placement  = "outside",
#         strip.background = element_rect(fill = NA, colour="grey50"),
#         axis.text.x      = element_text(angle = 45, hjust = 1),
#         axis.title.x     = element_text(size = 14, face = "bold"),
#         axis.title.y     = element_text(size = 14, face = "bold"),
#         panel.spacing    = unit(0,"cm"),
#         legend.position  = "bottom",
#         legend.title     = element_blank(),
#         legend.box       = "vertical",
#         legend.margin    = margin() )

###########################################################################
#
# Age-specific probability of infection (Rho_t) trajectories vs population
#
###########################################################################
# data_rho_cols     <- c("Date", "Group", "median", "low", "high", "low25", "high75")
# data_rho           <- data.frame(matrix(ncol = length(data_rho_cols), nrow = 0))
# colnames(data_rho) <- data_rho_cols
# 
# for (i in 1:cov_data$A){
#   
#   fit_rho  <- posts_1$rho_weekly[,,i] #-c(1:rm_dates)
#   #fit_rho  <- posts_1$Rt[,,i] #-c(1:rm_dates)
#   
#   dt_rho_country  <- data.frame(#Date     = dt_mortality_analysis_period$Date,
#                                Date   = dt_mortality_analysis_period$Date[seq(1, length(dt_mortality_analysis_period$Date), cov_data$gbm_no_days)], # Keep every 3rd element, due to having changes every 3 days
#                                Group  = rep(age_bands[i], dim(fit_rho)[2]) )
#   
#   # Add quantiles from the model outputs:
#   dt_rho_country$median <- apply(fit_rho, 2, median)
#   dt_rho_country$low    <- apply(fit_rho, 2, quantile, probs = c(0.025))
#   dt_rho_country$high   <- apply(fit_rho, 2, quantile, probs = c(0.975)) 
#   dt_rho_country$low25  <- apply(fit_rho, 2, quantile, probs = c(0.25)) 
#   dt_rho_country$high75 <- apply(fit_rho, 2, quantile, probs = c(0.75))
#   
#   data_rho  <- rbind(data_rho, dt_rho_country)
#   
# }# End for
# 
# #---- Check the transformed model outputs:
# table(data_rho$Group)
# 
# #---- Remove redundant objects:
# rm(i, dt_rho_country)
# 
# #---- Merge the age-specific dataset with the population dataste:
# df_rho_age_vs_population <- rbind(data_rho, transmrate_plot_data)
# 
# #---- Posterior median infection probability per group:
# ggplot(df_rho_age_vs_population,
#        aes(x      = Date,
#            y      = median,
#            colour = Group)) +
#   geom_line(size = 1.1) +
#   # labs(x = "Date",
#   #      y = "Estimated rho") +
#   labs(x = "Epidemiological Date",
#        y = "Effective contact rate") +
#   scale_x_date(date_breaks       = "2 weeks", 
#                date_minor_breaks = "2 weeks") + 
#   theme_bw() +
#   theme(panel.spacing    = unit(0.4,"cm"),
#         axis.text.x      = element_text(angle = 45, hjust = 1),
#         axis.title.x     = element_text(size = 14, face = "bold"),
#         axis.title.y     = element_text(size = 14, face = "bold"),
#         legend.position  = "bottom",
#         legend.title     = element_blank(),
#         legend.box       = "vertical", 
#         legend.margin    = margin())
# 
# #---- Plot of posterior estimates:
# 
# # Check that the estimated probabilities differs between groups at a given time-point:
# data_rho_subset <- subset(data_rho, Date == "2020-03-17")
#                   subset(data_rho, Date == "2020-06-21")
#                   subset(data_rho, Date == "2020-08-26")
# 
# ggplot(data_rho,
#        aes(Date       = Date,
#            New_Deaths = median)) +
#   facet_wrap(. ~ Group, 
#              scales = "free_y", 
#              ncol   = 1,
#              strip.position = "right") +
#   geom_line(aes(x     = Date,
#                 y     = median,
#                 color = "Posterior Median rho"),
#             size = 1.3) +
#   geom_ribbon(aes(x    = Date,
#                   ymin = low25,
#                   ymax = high75,
#                   fill = "50% CrI Estimated rho"),
#               alpha = 0.5) +
#   # geom_ribbon(aes(x    = Date,
#   #                 ymin = low,
#   #                 ymax = high,
#   #                 fill = "95% CrI Estimated rho"),
#   #             alpha = 0.5) +
#   labs(x = "Date",
#        y = "Estimated infection probability") +
#   # labs(x = "Epidemiological Date",
#   #      y = "Effective contact rate") +
#   scale_x_date(date_breaks       = "2 weeks", 
#                date_minor_breaks = "2 weeks") + 
#   scale_colour_manual(name   = '',
#                       values = c('Posterior Median rho' = "black")) +
#   scale_fill_manual(values = c("50% CrI Estimated rho" = "grey30"#,  
#                                #"95% CrI Estimated rho" = "grey50"
#   )
#   ) +
#   theme_bw() +
#   theme(panel.spacing    = unit(0.4,"cm"),
#         axis.text.x      = element_text(angle = 45, hjust = 1),
#         axis.title.x     = element_text(size = 14, face = "bold"),
#         axis.title.y     = element_text(size = 14, face = "bold"),
#         legend.position  = "bottom",
#         legend.title     = element_blank(),
#         legend.box       = "vertical", 
#         legend.margin    = margin())

###########################################################################
#
# Effective reproduction number trajectory from discrete renewal process
#
###########################################################################

#---- Plot of non-daily posterior estimates:
ggplot(post_Rt$Rt_nondaily_summary,
       aes(Date   = Date,
           Median = median)) +
  geom_line(aes(x     = Date,
                y     = median#,
                #color = "Median"
                ),
            size = 1.3) +
  # geom_ribbon(aes(x    = Date,
  #                 ymin = low25,
  #                 ymax = high75,
  #                 fill = "50% CrI"),
  #             alpha = 0.5) +
  # geom_ribbon(aes(x    = Date,
  #                 ymin = low,
  #                 ymax = high,
  #                 fill = "95% CrI"),
  #             alpha = 0.5) +
  geom_hline(yintercept = 1, color = "black") +
  
  labs(x = "Epidemiological Date",
       y = "Effective reproduction number") +
  scale_x_date(date_breaks       = "2 weeks",
               date_minor_breaks = "2 weeks") +
  # scale_y_continuous(
  #   limits = c( min(post_Rt$Rt_daily_summary$low25)*0.8,   max(post_Rt$Rt_daily_summary$high75)*1.1),    # max(death_plot_data$high)*1.2
  #   breaks = seq( 0, max(post_Rt$Rt_daily_summary$high75)*1.1, 0.5) #max(death_plot_data$high)*1.2
  # ) +
  scale_y_continuous(
    limits = c( 0,   6),    
    breaks = seq( 0, 6, 0.5)
  ) +
  # scale_fill_manual(values = c("50% CrI" = "gray40")) +
  # scale_colour_manual(name = '',
  #                     values = c('Median'  = "black"#,
  #                                #"50% CrI" = "gray70"
  #                                #"95% CrI Estimated cases" = "gray40",
  #                     )
  # ) +
  theme_bw() +
  theme(strip.placement  = "outside",
        strip.background = element_rect(fill = NA, colour="grey50"),
        axis.text.x      = element_text(angle = 45, hjust = 1),
        axis.title.x     = element_text(size = 14, face = "bold"),
        axis.title.y     = element_text(size = 14, face = "bold"),
        panel.spacing    = unit(0,"cm"),
        legend.position  = "bottom",
        legend.title     = element_blank(),
        legend.box       = "horizontal", 
        legend.margin    = margin() )

#---- Plot of daily posterior estimates:
ggplot(post_Rt$Rt_daily_summary,
       aes(Date   = Date,
           Median = median)) +
  geom_line(aes(x     = Date,
                y     = median#,
                #color = "Median"
  ),
  size = 1.3) +
  # geom_ribbon(aes(x    = Date,
  #                 ymin = low25,
  #                 ymax = high75,
  #                 fill = "50% CrI"),
  #             alpha = 0.5) +
  # geom_ribbon(aes(x    = Date,
  #                 ymin = low,
  #                 ymax = high,
  #                 fill = "95% CrI"),
  #             alpha = 0.5) +
  geom_hline(yintercept = 1, color = "black") +
  
  labs(x = "Epidemiological Date",
       y = "Effective reproduction number") +
  scale_x_date(date_breaks       = "2 weeks",
               date_minor_breaks = "2 weeks") +
  # scale_y_continuous(
  #   limits = c( min(post_Rt$Rt_daily_summary$low25)*0.8,   max(post_Rt$Rt_daily_summary$high75)*1.1),    # max(death_plot_data$high)*1.2
  #   breaks = seq( 0, max(post_Rt$Rt_daily_summary$high75)*1.1, 0.5) #max(death_plot_data$high)*1.2
  # ) +
  scale_y_continuous(
    limits = c( 0,   6),    
    breaks = seq( 0, 6, 0.5)
  ) +
  # scale_fill_manual(values = c("50% CrI" = "gray40")) +
  # scale_colour_manual(name = '',
  #                     values = c('Median'  = "black"#,
  #                                #"50% CrI" = "gray70"
  #                                #"95% CrI Estimated cases" = "gray40",
  #                     )
  # ) +
  theme_bw() +
  theme(strip.placement  = "outside",
        strip.background = element_rect(fill = NA, colour="grey50"),
        axis.text.x      = element_text(angle = 45, hjust = 1),
        axis.title.x     = element_text(size = 14, face = "bold"),
        axis.title.y     = element_text(size = 14, face = "bold"),
        panel.spacing    = unit(0,"cm"),
        legend.position  = "bottom",
        legend.title     = element_blank(),
        legend.box       = "horizontal", 
        legend.margin    = margin() )



###########################################################################
#
# Effective reproduction number trajectory from NGM
#
###########################################################################

#---- Plot of posterior estimates:
# ggplot(rt_ngm_plot_data,
#        aes(Date   = Date,
#            Median = median)) +
#   geom_line(aes(x     = Date,
#                 y     = median,
#                 color = "Median"),
#             size = 1.3) +
#   # geom_ribbon(aes(x    = Date,
#   #                 ymin = low25,
#   #                 ymax = high75,
#   #                 fill = "50% CrI"),
#   #             alpha = 0.5) +
#   # geom_ribbon(aes(x    = Date,
#   #                 ymin = low,
#   #                 ymax = high,
#   #                 fill = "95% CrI"),
#   #             alpha = 0.5) +
#   geom_hline(yintercept = 1, color = "black") +
#   
#   labs(x = "Epidemiological Date",
#        y = "Effective reproduction number (NGM)") +
#   scale_x_date(date_breaks       = "2 weeks",
#                date_minor_breaks = "2 week") +
#   #geom_hline(yintercept = 1, color = "black") +
#   # scale_y_continuous(
#   #   limits = c( min(rt_ngmplot_data$low25)*0.8,   max(rt_ngmplot_data$high75)*1.1),    # max(death_plot_data$high)*1.2
#   #   breaks = seq( 0.2, max(rt_ngmplot_data$high75)*1.1, 0.2) #max(death_plot_data$high)*1.2
#   # ) +
#   scale_colour_manual(name   = '',
#                       values = c('Median' = "black")) +
#   # scale_fill_manual(values = c("50% CrI" = "grey30"#,  #"steelblue3"
#   #                              #"95% CrI" = "steelblue1" #"steelblue1"
#   # )
#   # ) +
#   
#   theme_bw() +
#   theme(strip.placement  = "outside",
#         strip.background = element_rect(fill = NA, colour="grey50"),
#         axis.text.x      = element_text(angle = 45, hjust = 1),
#         axis.title.x     = element_text(size = 14, face = "bold"),
#         axis.title.y     = element_text(size = 14, face = "bold"),
#         panel.spacing    = unit(0,"cm"),
#         legend.position  = "bottom",
#         legend.title     = element_blank(),
#         legend.box       = "vertical",
#         legend.margin    = margin() )

###########################################################################
#
# Age-specific eff reproduction number trajectories from discrete renewal process
#
###########################################################################
data_rt_age_cols      <- c("Date", "Group", "median", "low", "high")
data_rt_age           <- data.frame(matrix(ncol = length(data_rt_age_cols), nrow = 0))
colnames(data_rt_age) <- data_rt_age_cols

for (i in 1:cov_data$A){
  
  fit_rt  <- posts_1$Rt_Age[,,i]
  
  dt_rt_country  <- data.frame(Date   = dt_mortality_analysis_period$Date, # Keep every 3rd element, due to having changes every 3 days
                               Group  = rep(age_bands[i], dim(fit_rt)[2]) )
  
  # Add quantiles from the model outputs:
  dt_rt_country$median <- apply(fit_rt, 2, median)
  dt_rt_country$low    <- apply(fit_rt, 2, quantile, probs = c(0.025))
  dt_rt_country$high   <- apply(fit_rt, 2, quantile, probs = c(0.975)) 
  dt_rt_country$low25  <- apply(fit_rt, 2, quantile, probs = c(0.25)) 
  dt_rt_country$high75 <- apply(fit_rt, 2, quantile, probs = c(0.75))
  
  # Store estimates every cov_data$gbm_no_days days:
  dt_rt_country <- dt_rt_country[seq(1, nrow(dt_rt_country), cov_data$gbm_no_days),]
  data_rt_age   <- rbind(data_rt_age, dt_rt_country)
  
}# End for

#---- Check the transformed model outputs:
table(data_rt_age$Group)

#---- Remove redundant objects:
rm(i, dt_rt_country)

#---- Plot of posterior estimates:

# Check that the estimated probabilities differs between groups at a given time-point:
data_rt_subset <- subset(data_rt_age, Date == "2020-03-17")
subset(data_rt, Date == "2020-06-21")
subset(data_rt, Date == "2020-08-26")

ggplot(data_rt_age,
       aes(Date       = Date,
           New_Deaths = median)) +
  facet_wrap(. ~ Group, 
             scales = "free_y", 
             ncol   = 2,
             strip.position = "right") +
  geom_line(aes(x     = Date,
                y     = median,
                color = "Posterior Median"),
            size = 1.3) +
  geom_ribbon(aes(x    = Date,
                  ymin = low25,
                  ymax = high75,
                  fill = "50% CrI"),
              alpha = 0.5) +
  # geom_ribbon(aes(x    = Date,
  #                 ymin = low,
  #                 ymax = high,
  #                 fill = "95% CrI Estimated rho"),
  #             alpha = 0.5) +
  # labs(x = "Date",
  #      y = "Estimated rho") +
  labs(x = "Epidemiological Date",
       y = "Effective reproduction number") +
  scale_x_date(date_breaks       = "2 weeks", 
               date_minor_breaks = "2 weeks") + 
  scale_colour_manual(name   = '',
                      values = c('Posterior Median' = "black")) +
  scale_fill_manual(values = c("50% CrI Estimated" = "grey30"#,  
                               #"95% CrI Estimated rho" = "grey50"
  )
  ) +
  theme_bw() +
  theme(panel.spacing    = unit(0.4,"cm"),
        axis.text.x      = element_text(angle = 45, hjust = 1),
        axis.title.x     = element_text(size = 14, face = "bold"),
        axis.title.y     = element_text(size = 14, face = "bold"),
        legend.position  = "bottom",
        legend.title     = element_blank(),
        legend.box       = "vertical", 
        legend.margin    = margin())

#---- Posterior median R_t per group:
ggplot(data_rt_age,
       aes(x      = Date,
           y      = median,
           colour = Group)) +
  geom_line(size = 1.1) +
  geom_hline(yintercept = 1, color = "black", lty = 2) +
  labs(x = "Epidemiological Date",
       y = "Effective reproduction number") +
  scale_x_date(date_breaks       = "1 weeks", 
               date_minor_breaks = "1 weeks") + 
  theme_bw() +
  theme(panel.spacing    = unit(0.4,"cm"),
        axis.text.x      = element_text(angle = 45, hjust = 1),
        axis.title.x     = element_text(size = 14, face = "bold"),
        axis.title.y     = element_text(size = 14, face = "bold"),
        legend.position  = "bottom",
        legend.title     = element_blank(),
        legend.box       = "vertical", 
        legend.margin    = margin())
