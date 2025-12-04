lib <- c("rstan",
         "Bernadette",
         "bridgesampling"
         )
lapply(lib, require, character.only = TRUE)

main_path            <- "C://Users//lbour//OneDrive//Desktop1//xgp_project//Zika//"
store_model_out_iGP  <- paste0(main_path, "4_Stan_Outputs//", "CHIKV_FWI_2014_iGP",  ".RData")
store_model_out_xGP  <- paste0(main_path, "4_Stan_Outputs//", "CHIKV_FWI_2014_xGP",  ".RData")
store_model_out_mxGP <- paste0(main_path, "4_Stan_Outputs//", "CHIKV_FWI_2014_mxGP", ".RData")

load(paste0(main_path, "1_Data//zikachik.Rdata"))
source(paste0(main_path, "2_Stan_files//priors.R"))
`%nin%` <- Negate(`%in%`)

###########################################################################
#
# Region: French West Indies, CHIKV disease -----------------
#
###########################################################################

## Select ZIKV in the French polynesian islands:
key_CHIKV_filtered <- key %>% filter(REGION_ID == 1)

zikachik$DATE2 <- lubridate::ymd(zikachik$DATE)

main_dataset_CHIKV <- 
  zikachik %>% 
  filter(REGION_ID == 1,
         VIRUS == "CHIKV",
         DATE2 <= "2015-04-26"
  ) %>% 
  arrange(ISLAND_ID, DATE2)

## Locate the min and max dates per island, so as to keep the common records:

main_dataset_CHIKV %>% 
  dplyr::group_by(ISLAND) %>% 
  dplyr::summarize(
    date_min = min(DATE2, na.rm=T),
    date_max = max(DATE2, na.rm=T)
  )

## Decision: keep records between 2014-11-09 and 2015-02-22.
main_dataset_CHIKV <- 
  zikachik %>% 
  filter(REGION_ID == 1,
         VIRUS == "CHIKV",
         DATE2 >= "2014-01-12" & DATE2 <= "2014-12-07"
  ) %>% 
  arrange(ISLAND_ID, DATE2)

key_CHIKV_filtered <-
  key_CHIKV_filtered %>% 
  filter(ISLAND %nin% "GUYANE")

###########################################################################
#
# Plot the clinical cases per island:
#
###########################################################################
ggplot(main_dataset_CHIKV,
       aes(date       = DATE2,
           New_Deaths = NCASES)) +
  geom_point(aes(x   = DATE2,
                 y   = NCASES,
                 fill= "Reported cases")) +
  facet_wrap(. ~ ISLAND, #ISLAND_ID
             scales = "free_y", 
             ncol   = 3,
             strip.position = "top",
             shrink = T) +
  scale_x_date(date_breaks = "2 weeks"#,
               #date_labels =  "%b %Y"
  ) +
  scale_fill_manual(values = c('Reported cases' = 'black'), guide = "none") +
  labs(x = "Date",
       y = "CHIKV clinical cases") +
  theme_bw() +
  theme(panel.spacing    = unit(0.4,"cm"),
        axis.text.x      = element_text(size = 14, angle = 45, hjust = 1),
        axis.text.y      = element_text(size = 16),
        axis.title.x     = element_text(size = 18, face = "bold"),
        axis.title.y     = element_text(size = 18, face = "bold"),
        legend.text      = element_text(size = 16),
        strip.text.x     = element_text(size = 14),
        legend.title     = element_blank(),
        legend.box       = "vertical", 
        legend.margin    = margin())

###########################################################################
#
# Prepare the datasets for inclusion in Stan:
#
###########################################################################
O_t        <- main_dataset_CHIKV %>% 
  select(ISLAND_ID, DATE2, NCASES) %>% 
  pivot_wider(names_from = ISLAND_ID, values_from = NCASES)

Ostar_t    <- main_dataset_CHIKV %>% 
  select(ISLAND_ID, DATE2, Ostar) %>% 
  pivot_wider(names_from = ISLAND_ID, values_from = Ostar)

sumO_t    <-  main_dataset_CHIKV %>% 
  select(ISLAND_ID, DATE2, CUM_NCASES) %>% 
  pivot_wider(names_from = ISLAND_ID, values_from = CUM_NCASES)

pop       <- key_CHIKV_filtered$POP # This is ordered by ISLAND_ID right from the beginning

weather_prec_names <- colnames(main_dataset_CHIKV)[grepl("PrecCm_", colnames(main_dataset_CHIKV), fixed = TRUE)]

C         <- length(weather_prec_names)

weather   <- vector(mode = "list", length = nrow(key_CHIKV_filtered))
for (i in 1:nrow(key_CHIKV_filtered)) {
  temp         <- subset(main_dataset_CHIKV[,c("DATE2", "ISLAND_ID", weather_prec_names)], ISLAND_ID == i+6)
  weather[[i]] <- temp[,weather_prec_names]
  colnames( weather[[i]] ) <- paste0(i,"_",colnames( weather[[i]] ))
}

weather <- do.call(cbind, weather)

###########################################################################
#
# Stan model, data and initialisation:
#
###########################################################################
compilation_time_start <- Sys.time()
m_igp  <- rstan::stan_model(paste0(main_path, "2_Stan_files//hierarchical_TSIR_IGP.stan"))
m_xgp  <- rstan::stan_model(paste0(main_path, "2_Stan_files//hierarchical_TSIR_XGP.stan"))
m_mxgp <- rstan::stan_model(paste0(main_path, "2_Stan_files//hierarchical_TSIR_mXGP.stan"))
compilation_time_end <- Sys.time()
duration_compilation <- compilation_time_end - compilation_time_start

###########################################################################
#
# HMC initialisation:
#
###########################################################################
inference                <- 1
likelihood_variance_type <- "linear"

if( likelihood_variance_type == "quadratic") l_variance_type <- 0 else if(likelihood_variance_type == "linear") l_variance_type <- 1

standata_preprocessed <-
  nlist(W                         = nrow(O_t),
        K                         = nrow(key_CHIKV_filtered),
        O_t                       = O_t[, -c(1)],
        Ostar_t                   = Ostar_t[, -c(1)],
        sumO_t                    = sumO_t[, -c(1)],
        pop                       = pop, 
        C                         = C,
        weather                   = weather,
        likelihood_variance_type  = l_variance_type,
        inference                 = 1
  )

ok_dists <- nlist("normal",
                  student_t = "t",
                  "cauchy",
                  "gamma",
                  "exponential")

parallel::detectCores()
#########################################################
#
# iGP
#
#########################################################
prior_sigma         <- Bernadette::normal(location = 0, scale = 1)
prior_nb_dispersion <- Bernadette::normal(location = 0, scale = 1)

prior_dist_sigma             <-
  prior_mean_sigma          <- prior_scale_sigma      <- prior_df_sigma  <-
  prior_shape_sigma         <- prior_rate_sigma       <-
  prior_dist_nb_dispersion  <-
  prior_mean_nb_dispersion  <- prior_scale_nb_dispersion <- prior_df_nb_dispersion <-
  prior_shape_nb_dispersion <- prior_rate_nb_dispersion  <- NULL

#---- Prior distribution for the volatilities sigma (handle_prior() from priors.R):
prior_params_sigma        <- handle_prior(prior_sigma, ok_dists = ok_dists)
names(prior_params_sigma) <- paste0(names(prior_params_sigma),  "_sigma")

for (i in names(prior_params_sigma)) assign(i, prior_params_sigma[[i]])

#---- Prior distribution for the dispersion (handle_prior() from priors.R):
prior_params_dispersion        <- handle_prior(prior_nb_dispersion, ok_dists = ok_dists)
names(prior_params_dispersion) <- paste0(names(prior_params_dispersion),  "_nb_dispersion")

for (i in names(prior_params_dispersion)) assign(i, prior_params_dispersion[[i]])

#---- Create entries in the data block of the .stan file:
standata <- c(standata_preprocessed,
              nlist(
                ###
                prior_dist_sigma     = prior_dist_sigma,
                prior_mean_sigma     = prior_mean_sigma,
                prior_scale_sigma    = prior_scale_sigma,
                prior_df_sigma       = prior_df_sigma,
                prior_shape_sigma    = prior_shape_sigma,
                prior_rate_sigma     = prior_rate_sigma,
                ###
                prior_dist_nb_dispersion  = prior_dist_nb_dispersion,
                prior_mean_nb_dispersion  = prior_mean_nb_dispersion,
                prior_scale_nb_dispersion = prior_scale_nb_dispersion,
                prior_df_nb_dispersion    = prior_df_nb_dispersion,
                prior_shape_nb_dispersion = prior_shape_nb_dispersion,
                prior_rate_nb_dispersion  = prior_rate_nb_dispersion,
                timepoints                = 1:standata_preprocessed$W))

#---- Set initial values:
fit_optim    <- rstan::optimizing(m_igp, data = standata, seed = 1, hessian = TRUE)
sampler_init <- function(){
  list(x0           = fit_optim$par[grepl("x_gp_raw[",     names(fit_optim$par), fixed = TRUE)],
       x_noise      = fit_optim$par[grepl("x_gp_common[",  names(fit_optim$par), fixed = TRUE)],
       sigma        = fit_optim$par[grepl("sigma[",        names(fit_optim$par), fixed = TRUE)],
       length_scale = fit_optim$par[grepl("length_scale[", names(fit_optim$par), fixed = TRUE)],
       bW           = fit_optim$par[grepl("bW[",           names(fit_optim$par), fixed = TRUE)],
       rho          = fit_optim$par[names(fit_optim$par) %in% "phi"]
  )
}

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
rstan_options(auto_write = TRUE)
options(mc.cores = nChains)

time_start_nuts1 <- Sys.time()
nuts_fit_iGP <- rstan::sampling(m_igp,
                                data    = standata,
                                init    = sampler_init,
                                chains  = nChains,
                                warmup  = nBurnin,
                                iter    = nIter,
                                seed    = 1,
                                thin    = nThin,
                                control = list(max_treedepth = max_treedepth,
                                               adapt_delta   = adapt_delta
                                ),
                                show_messages   = FALSE)
time_end_nuts1 <- Sys.time()
duration_nuts1 <- time_end_nuts1 - time_start_nuts1

time_run_igp <- get_elapsed_time(nuts_fit_iGP)
apply(time_run_igp, 1, sum)

bridgeiGP <- bridgesampling::bridge_sampler(samples = nuts_fit_iGP)

save(standata,
     nuts_fit_iGP,
     time_run_igp,
     bridgeiGP,
     file = store_model_out_iGP)

nuts_fit_1_summary <- summary(nuts_fit_iGP
                              )$summary

rm(standata, nuts_fit_iGP, sampler_init, fit_optim)
load(store_model_out_iGP)

#########################################################
#
# xGP
#
#########################################################
prior_sigma         <- Bernadette::normal(location = 0, scale = 1)
prior_sigma_common  <- Bernadette::normal(location = 0, scale = 1)
prior_nb_dispersion <- Bernadette::normal(location = 0, scale = 1)

#---- Useless assignments to pass R CMD check
prior_dist_sigma_common         <-
  prior_mean_sigma_common       <- prior_scale_sigma_common      <- prior_df_sigma_common  <-
  prior_shape_sigma_common      <- prior_rate_sigma_common       <-
  prior_dist_sigma          <-
  prior_mean_sigma        <- prior_scale_sigma       <- prior_df_sigma   <-
  prior_shape_sigma       <- prior_rate_sigma        <-
  prior_dist_nb_dispersion    <-
  prior_mean_nb_dispersion  <- prior_scale_nb_dispersion <- prior_df_nb_dispersion <-
  prior_shape_nb_dispersion <- prior_rate_nb_dispersion  <- NULL

ok_dists <- nlist("normal",
                  student_t = "t",
                  "cauchy",
                  "gamma",
                  "exponential")

#---- Prior distribution for the volatilities sigma_mu and sigma_x (handle_prior() from priors.R):
prior_params_sigma_common        <- handle_prior(prior_sigma_common, ok_dists = ok_dists)
names(prior_params_sigma_common) <- paste0(names(prior_params_sigma_common),  "_sigma_common")

for (i in names(prior_params_sigma_common)) assign(i, prior_params_sigma_common[[i]])

prior_params_sigma        <- handle_prior(prior_sigma, ok_dists = ok_dists)
names(prior_params_sigma) <- paste0(names(prior_params_sigma),  "_sigma")

for (i in names(prior_params_sigma)) assign(i, prior_params_sigma[[i]])

#---- Prior distribution for the dispersion (handle_prior() from priors.R):
prior_params_dispersion        <- handle_prior(prior_nb_dispersion, ok_dists = ok_dists)
names(prior_params_dispersion) <- paste0(names(prior_params_dispersion),  "_nb_dispersion")

for (i in names(prior_params_dispersion)) assign(i, prior_params_dispersion[[i]])

#---- Create entries in the data block of the .stan file:
standata <- c(standata_preprocessed,
              nlist(
                ###
                prior_dist_sigma_common     = prior_dist_sigma_common,
                prior_mean_sigma_common     = prior_mean_sigma_common,
                prior_scale_sigma_common    = prior_scale_sigma_common,
                prior_df_sigma_common       = prior_df_sigma_common,
                prior_shape_sigma_common    = prior_shape_sigma_common,
                prior_rate_sigma_common     = prior_rate_sigma_common,
                ###
                prior_dist_sigma     = prior_dist_sigma,
                prior_mean_sigma     = prior_mean_sigma,
                prior_scale_sigma    = prior_scale_sigma,
                prior_df_sigma       = prior_df_sigma,
                prior_shape_sigma    = prior_shape_sigma,
                prior_rate_sigma     = prior_rate_sigma,
                ###
                prior_dist_nb_dispersion  = prior_dist_nb_dispersion,
                prior_mean_nb_dispersion  = prior_mean_nb_dispersion,
                prior_scale_nb_dispersion = prior_scale_nb_dispersion,
                prior_df_nb_dispersion    = prior_df_nb_dispersion,
                prior_shape_nb_dispersion = prior_shape_nb_dispersion,
                prior_rate_nb_dispersion  = prior_rate_nb_dispersion,
                timepoints                = 1:standata_preprocessed$W
              ))

#---- Set initial values:
fit_optim    <- rstan::optimizing(m_xgp, data = standata, seed = 1, hessian = TRUE)
sampler_init <- function(){
  list(x0           = fit_optim$par[grepl("x_gp_raw[",     names(fit_optim$par), fixed = TRUE)],
       x_noise      = fit_optim$par[grepl("x_gp_common[",  names(fit_optim$par), fixed = TRUE)],
       sigma        = fit_optim$par[grepl("sigma_x",       names(fit_optim$par), fixed = TRUE)],
       length_scale = fit_optim$par[grepl("length_scale_x",names(fit_optim$par), fixed = TRUE)],
       bW           = fit_optim$par[grepl("bW[",           names(fit_optim$par), fixed = TRUE)],
       sigma_mu     = fit_optim$par[names(fit_optim$par) %in% "sigma_common"],
       phiD         = fit_optim$par[names(fit_optim$par) %in% "length_scale_common"],
       rho          = fit_optim$par[names(fit_optim$par) %in% "phi"]
  )
}

options(mc.cores = nChains)
time_start_nuts2<- Sys.time()
nuts_fit_xGP <- rstan::sampling(m_xgp,
                                data    = standata,
                                init    = sampler_init,
                                chains  = nChains,
                                warmup  = nBurnin,
                                iter    = nIter,
                                seed    = 1,
                                thin    = nThin,
                                control = list(max_treedepth = max_treedepth,
                                               adapt_delta   = adapt_delta
                                ),
                                show_messages   = FALSE)
time_end_nuts2 <- Sys.time()
duration_nuts2 <- time_end_nuts2 - time_start_nuts2

time_run_xgp <- get_elapsed_time(nuts_fit_xGP)
apply(time_run_xgp, 1, sum)

bridgexGP <- bridgesampling::bridge_sampler(samples = nuts_fit_xGP)

save(standata,
     nuts_fit_xGP,
     time_run_xgp,
     bridgexGP,
     file = store_model_out_xGP)

nuts_fit_1_summary <- summary(nuts_fit_xGP)$summary

rm(standata, nuts_fit_xGP, sampler_init, fit_optim)
load(store_model_out_xGP)

#########################################################
#
# mxGP
#
#########################################################

prior_sigma         <- Bernadette::normal(location = 0, scale = 1)
prior_sigma_common  <- Bernadette::normal(location = 0, scale = 1)
prior_nb_dispersion <- Bernadette::normal(location = 0, scale = 1)

#---- Useless assignments to pass R CMD check
prior_dist_sigma_common         <-
  prior_mean_sigma_common       <- prior_scale_sigma_common      <- prior_df_sigma_common  <-
  prior_shape_sigma_common      <- prior_rate_sigma_common       <-
  prior_dist_sigma          <-
  prior_mean_sigma        <- prior_scale_sigma       <- prior_df_sigma   <-
  prior_shape_sigma       <- prior_rate_sigma        <-
  prior_dist_nb_dispersion    <-
  prior_mean_nb_dispersion  <- prior_scale_nb_dispersion <- prior_df_nb_dispersion <-
  prior_shape_nb_dispersion <- prior_rate_nb_dispersion  <- NULL

ok_dists <- nlist("normal",
                  student_t = "t",
                  "cauchy",
                  "gamma",
                  "exponential")

#---- Prior distribution for the volatilities sigma_mu and sigma_x (handle_prior() from priors.R):
prior_params_sigma_common        <- handle_prior(prior_sigma_common, ok_dists = ok_dists)
names(prior_params_sigma_common) <- paste0(names(prior_params_sigma_common),  "_sigma_common")

for (i in names(prior_params_sigma_common)) assign(i, prior_params_sigma_common[[i]])

prior_params_sigma        <- handle_prior(prior_sigma, ok_dists = ok_dists)
names(prior_params_sigma) <- paste0(names(prior_params_sigma),  "_sigma")

for (i in names(prior_params_sigma)) assign(i, prior_params_sigma[[i]])

#---- Prior distribution for the dispersion (handle_prior() from priors.R):
prior_params_dispersion        <- handle_prior(prior_nb_dispersion, ok_dists = ok_dists)
names(prior_params_dispersion) <- paste0(names(prior_params_dispersion),  "_nb_dispersion")

for (i in names(prior_params_dispersion)) assign(i, prior_params_dispersion[[i]])

#---- Create entries in the data block of the .stan file:
standata <- c(standata_preprocessed,
              nlist(
                ###
                prior_dist_sigma_common     = prior_dist_sigma_common,
                prior_mean_sigma_common     = prior_mean_sigma_common,
                prior_scale_sigma_common    = prior_scale_sigma_common,
                prior_df_sigma_common       = prior_df_sigma_common,
                prior_shape_sigma_common    = prior_shape_sigma_common,
                prior_rate_sigma_common     = prior_rate_sigma_common,
                ###
                prior_dist_sigma     = prior_dist_sigma,
                prior_mean_sigma     = prior_mean_sigma,
                prior_scale_sigma    = prior_scale_sigma,
                prior_df_sigma       = prior_df_sigma,
                prior_shape_sigma    = prior_shape_sigma,
                prior_rate_sigma     = prior_rate_sigma,
                ###
                prior_dist_nb_dispersion  = prior_dist_nb_dispersion,
                prior_mean_nb_dispersion  = prior_mean_nb_dispersion,
                prior_scale_nb_dispersion = prior_scale_nb_dispersion,
                prior_df_nb_dispersion    = prior_df_nb_dispersion,
                prior_shape_nb_dispersion = prior_shape_nb_dispersion,
                prior_rate_nb_dispersion  = prior_rate_nb_dispersion,
                timepoints                = 1:standata_preprocessed$W
              ))

#---- Set initial values:
fit_optim    <- rstan::optimizing(m_mxgp, data = standata, seed = 1, hessian = TRUE)
sampler_init <- function(){
  list(x0           = fit_optim$par[grepl("x_gp_raw[",     names(fit_optim$par), fixed = TRUE)],
       x_noise      = fit_optim$par[grepl("x_gp_common[",  names(fit_optim$par), fixed = TRUE)],
       sigma        = fit_optim$par[grepl("sigma[",        names(fit_optim$par), fixed = TRUE)],
       length_scale = fit_optim$par[grepl("length_scale[", names(fit_optim$par), fixed = TRUE)],
       bW           = fit_optim$par[grepl("bW[",           names(fit_optim$par), fixed = TRUE)],
       sigma_mu     = fit_optim$par[names(fit_optim$par) %in% "sigma_common"],
       phiD         = fit_optim$par[names(fit_optim$par) %in% "length_scale_common"],
       rho          = fit_optim$par[names(fit_optim$par) %in% "phi"]
  )
}

options(mc.cores = nChains)
time_start_nuts3<- Sys.time()
nuts_fit_mxGP <- rstan::sampling(m_mxgp,
                                 data    = standata,
                                 init    = sampler_init,
                                 chains  = nChains,
                                 warmup  = nBurnin,
                                 iter    = nIter,
                                 seed    = 1,
                                 thin    = nThin,
                                 control = list(max_treedepth = max_treedepth,
                                                adapt_delta   = adapt_delta
                                 ),
                                 show_messages   = FALSE)
time_end_nuts3 <- Sys.time()
duration_nuts3 <- time_end_nuts3 - time_start_nuts3

time_run_mxgp <- get_elapsed_time(nuts_fit_mxGP)
apply(time_run_mxgp, 1, sum)

bridgemxGP <- bridgesampling::bridge_sampler(samples = nuts_fit_mxGP, maxiter = 10000)

save(standata,
     nuts_fit_mxGP,
     time_run_mxgp,
     bridgemxGP,
     file = store_model_out_mxGP)

nuts_fit_1_summary <- summary(nuts_fit_mxGP)$summary
load(store_model_out_mxGP)

###########################################################################
#
# Store HMC output and MCMC summaries
#
###########################################################################
print(bridgeiGP)
bridgesampling::error_measures(bridgeiGP)$re2

print(bridgexGP)
bridgesampling::error_measures(bridgexGP)$re2

print(bridgemxGP)
bridgesampling::error_measures(bridgemxGP)$re2

# More stable with warp3, but cannot estimate the error:
bridgemxGP2 <- bridgesampling::bridge_sampler(samples = nuts_fit_mxGP, maxiter = 5000, method = "warp3")
print(bridgemxGP2)
