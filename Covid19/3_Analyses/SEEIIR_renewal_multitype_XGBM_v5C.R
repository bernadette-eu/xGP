###########################################################################
#
# Paths
#
###########################################################################
Sys.setenv(HOME = "C:/R_Home") # See https://gemini.google.com/app/581c2a4c864865ef about the greek letter path issue

pc_path        <- "C://Users//lbour//OneDrive//Desktop1//"

###########################################################################
#
# Load libraries:
#
###########################################################################
lib <- c("ggplot2",
         "tidyverse",
         "rvest",
         "vroom",
         "rstan",
         #"cmdstanr",
         "bayesplot",
         "gridExtra",
         "Bernadette",
         "readxl",
         "compiler",
         "extraDistr",
         "bridgesampling"
)
lapply(lib, require, character.only = TRUE)

#---- Set system locale to English:
Sys.setlocale("LC_ALL", "English")

###########################################################################
#
# Setup:
#
###########################################################################

sessionInfo()

#--- Paths:
main_path       <- paste0(pc_path,   "xGP_project//GP multitype epidemic England//")
stan_file_path  <- paste0(main_path, "2_Stan_files//Revision1//",   "SEEIIR_renewal_multitype_XGBM_v4C", ".stan")
store_model_out <- paste0(main_path, "4_Stan_Outputs//Revision1//", "SEEIIR_renewal_multitype_MXGBM", ".RData")

load(paste0(main_path, "4_Stan_Outputs//Submission//SEEIIR_renewal_multitype_XGBM_v5C.RData"))

# Update the IFR to have a better-looking validation graph:
#cov_data$ifr <- c(0.07, 0.021, 0.22, 7.2)/100 #3D
cov_data$ifr <- c(0.009, 0.014, 0.19, 7.2)/100 #3E

###########################################################################
#
# Stan model, data and initialisation:
#
###########################################################################

#---- Stan options:
rstan_options(auto_write = TRUE)

#---- Import the Stan model:
compilation_time_start <- Sys.time()
m1 <- rstan::stan_model(stan_file_path)
compilation_time_end <- Sys.time()
duration_compilation <- compilation_time_end - compilation_time_start
# 
# cov_data <- list(A             = length(ifr$Group_mapping),
#                  N             = nrow(dt_mortality_analysis_period),
#                  y_data        = dt_mortality_analysis_period[,-c(1:5)],
#                  n_weeks       = 70,
#                  gbm_no_days   = 3,
#                  n_pop         = aggr_age$PopTotal,
#                  N0            = 6,
#                  EpidemicStart = 7,
#                  # Infection-to-death distribution:
#                  ifr           = ifr[,-1],
#                  # Generation time interval:
#                  I_D           = ditd,
#                  GT            = gen_time2,
#                  contact_matrix= aggr_cm,
#                  # Priors:
#                  #p_tau         = 0.03,
#                  p_phiD        = 5,
#                  p_sigma_x     = rep(2, nrow(ifr)),
#                  p_sigma_mu    = 2,
#                  # Debugging:
#                  inference     = 1,
#                  doprint       = 0)

cov_data <- c(cov_data, dI = 4)

#---- Specify parameters to monitor:
parameters <- c("phiD",
                "sigma_x",
                "sigma_mu",
                "rho_weekly",
                'rho_daily',
                "E_casesAge",
                "E_deathsAge",               
                "E_cases",
                "E_deaths",
                "log_like_age",
                "log_lik",
                "deviance",
                "Susceptibles"
)

#---- Set initial values:
fit_optim <- rstan::optimizing(m1,
                               data    = cov_data,
                               seed    = 50,
                               hessian = TRUE)

sampler_init <- function(){
  list(#tau         = fit_optim$par[grepl("tau",        names(fit_optim$par), fixed = TRUE)],
       E_cases_N0  = fit_optim$par[grepl("E_cases_N0[",names(fit_optim$par), fixed = TRUE)],
       phiD        = fit_optim$par[grepl("phiD",       names(fit_optim$par), fixed = TRUE)],
       eta_noise   = fit_optim$par[grepl("eta_noise[", names(fit_optim$par), fixed = TRUE)],
       sigma_x     = fit_optim$par[grepl("sigma_x[",   names(fit_optim$par), fixed = TRUE)],
       sigma_mu    = fit_optim$par[grepl("sigma_mu",   names(fit_optim$par), fixed = TRUE)]
  )
}# End function

###########################################################################
#
# Fit and sample from the posterior using NUTS:
#
###########################################################################

parallel::detectCores()

#---- MCMC options:
nChains       <- 4
nBurn         <- 1500 ## Number of warm-up samples per chain after thinning
nPost         <- 1500 ## Number of post-warm-up samples per chain after thinning
nThin         <- 1
adapt_delta   <- 0.99
max_treedepth <- 19

nIter   <- (nBurn + nPost) * nThin
nBurnin <- nBurn * nThin

#---- Stan options:
options(mc.cores = nChains)

#---- Test execution:
# test <- rstan::sampling(m1,
#                         data    = cov_data,
#                         init    = sampler_init,
#                         pars    = parameters,
#                         chains  = 1,
#                         warmup  = 50,
#                         iter    = 100,
#                         seed    = 1,
#                         control = list(max_treedepth = max_treedepth,
#                                        adapt_delta   = adapt_delta#,
#                                        #metric        = "dense_e"
#                         ),
#                         show_messages   = TRUE)
# 
# rm(test)

time_start_nuts1 <- Sys.time()
nuts_fit_1 <- rstan::sampling(m1,
                              data    = cov_data,
                              init    = sampler_init,
                              #pars    = parameters,
                              chains  = nChains,
                              warmup  = nBurnin,
                              iter    = nIter,
                              seed    = 1,
                              control = list(max_treedepth = max_treedepth,
                                             adapt_delta   = adapt_delta#,
                                             #metric        = "dense_e"
                              ),
                              show_messages   = TRUE)
time_end_nuts1 <- Sys.time()
duration_nuts1 <- time_end_nuts1 - time_start_nuts1

# Warmup and sampling times for each chain:
# https://cran.r-project.org/web/packages/rstan/vignettes/stanfit-objects.html
print( get_elapsed_time(nuts_fit_1) )

 

###########################################################################
#
# Store HMC output and MCMC summaries
#
###########################################################################
save(cov_data,
     sampler_init,
     nuts_fit_1,
     duration_nuts1,
     file = store_model_out)

load(file = store_model_out)

