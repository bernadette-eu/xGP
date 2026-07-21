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
lib <- c("dplyr",
         "rstan",
         "lubridate"
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
main_path       <- paste0(pc_path, "xGP_project//GP multitype epidemic England//")
stan_file_path  <- paste0(main_path, "2_Stan_files//Revision1//", "mxGP one country and multiple age groups SE GP - Optimized 1 - with day cp", ".stan")
store_model_out <- paste0(main_path, "4_Stan_Outputs//Revision1//", "mxGP", ".RData") # 3D where the ifr$AgrIFR <- c(0.07, 0.021, 0.22, 7.2)/100 are considered

load(paste0(main_path, "4_Stan_Outputs//Submission//SEEIIR_renewal_multitype_IGBM_v3C.RData"))
rm(nuts_fit_1)

#--- Intervals death and generation
cov_data$ifr   <- c(0.009, 0.014, 0.19, 7.2)/100

data_ENG2 <- list(gp_days_cp = 3,
                  N0         = cov_data$N0,
                  T          = dim(cov_data$y_data)[1],
                  N_age_groups = dim(cov_data$y_data)[2],
                  deaths     = t(cov_data$y_data),
                  Cont_mat   = cov_data$contact_matrix,
                  pop        = cov_data$n_pop,
                  s_int_rev1 = rev(cov_data$GT),
                  s_int_rev2 = rev(cov_data$I_D),
                  ifr        = cov_data$ifr,
                  epi_start  = cov_data$EpidemicStart,
                  dI         = 4)

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

#---- Set initial values:
fit_optim <- rstan::optimizing(m1,
                               data    = data_ENG2,
                               seed    = 10,
                               hessian = TRUE)
# 
# sampler_init <- function(){
#   list(
#     cases0           = fit_optim$par[grepl("cases0[",           names(fit_optim$par), fixed = TRUE)],
#     tau              = fit_optim$par[grepl("tau",               names(fit_optim$par), fixed = TRUE)],
#     phi              = fit_optim$par[grepl("phi",               names(fit_optim$par), fixed = TRUE)],
#     eta_country      = fit_optim$par[grepl("eta_country[",      names(fit_optim$par), fixed = TRUE)],
#     eta_age_group    = fit_optim$par[grepl("eta_age_group[",    names(fit_optim$par), fixed = TRUE)],
#     sigma_age_group  = fit_optim$par[grepl("sigma_age_group[",  names(fit_optim$par), fixed = TRUE)],
#     lambda_age_group = fit_optim$par[grepl("lambda_age_group[", names(fit_optim$par), fixed = TRUE)],
#     sigma_country    = fit_optim$par[grepl("sigma_country",     names(fit_optim$par), fixed = TRUE)],
#     lambda_country   = fit_optim$par[grepl("lambda_country",    names(fit_optim$par), fixed = TRUE)]
#   )
# }# End function

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
#                         #init    = sampler_init,
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
                              data    = data_ENG2,
                              init    = 0,
                              chains  = nChains,
                              warmup  = nBurnin,
                              iter    = nIter,
                              seed    = 1,
                              control = list(max_treedepth = max_treedepth,
                                             adapt_delta   = adapt_delta
                              ),
                              show_messages   = TRUE)
time_end_nuts1 <- Sys.time()
duration_nuts1 <- time_end_nuts1 - time_start_nuts1

# Warmup and sampling times for each chain:
# https://cran.r-project.org/web/packages/rstan/vignettes/stanfit-objects.html
print( get_elapsed_time(nuts_fit_1) )

time_run <- get_elapsed_time(nuts_fit_1)
apply(time_run, 1, sum)

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