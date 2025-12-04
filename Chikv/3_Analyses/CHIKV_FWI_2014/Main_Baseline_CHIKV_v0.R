###########################################################################
#
# Libraries
#
###########################################################################
lib <- c("tidyverse",
         "lubridate",
         "rstan",
         "Bernadette",
         "bridgesampling"
)
lapply(lib, require, character.only = TRUE)

###########################################################################
#
# Paths
#
###########################################################################
main_path       <- "...//xgp_project//Zika//"
data_path       <- paste0(main_path, "1_Data//", "zikachik", ".RData")
stan_file_path  <- paste0(main_path, "2_Stan_files//",   "hierarchical_TSIR_v0", ".stan")
store_model_out <- paste0(main_path, "4_Stan_Outputs//", "Main_Model_CHIKV_FWI_v1", ".RData")
path_to_package <- paste0(main_path, "2_Stan_files//")

source(paste0(path_to_package, "priors", ".R"))
`%nin%` <- Negate(`%in%`)

load(data_path)

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
  scale_x_date(date_breaks = "4 weeks"#,
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
        legend.position  = "bottom",
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
m1 <- rstan::stan_model(stan_file_path)
compilation_time_end <- Sys.time()
duration_compilation <- compilation_time_end - compilation_time_start

###########################################################################
#
# HMC initialisation:
#
###########################################################################
inference <- 1

likelihood_variance_type    <- "linear"
prior_sigma_mu              <- Bernadette::normal(location = 0, scale = 1)
prior_sigma_x               <- Bernadette::normal(location = 0, scale = 1)
prior_nb_dispersion         <- Bernadette::normal(location = 0, scale = 1)

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
        inference                 = 1,
        doprint = 0
  )

#---- Useless assignments to pass R CMD check
prior_dist_nb_dispersion    <-
  prior_mean_nb_dispersion  <- prior_scale_nb_dispersion <- prior_df_nb_dispersion <-
  prior_shape_nb_dispersion <- prior_rate_nb_dispersion  <- NULL

ok_dists <- nlist("normal",
                  student_t = "t",
                  "cauchy",
                  "gamma",
                  "exponential")

#---- Prior distribution for the dispersion (handle_prior() from priors.R):
prior_params_dispersion        <- handle_prior(prior_nb_dispersion, ok_dists = ok_dists)
names(prior_params_dispersion) <- paste0(names(prior_params_dispersion),  "_nb_dispersion")

for (i in names(prior_params_dispersion)) assign(i, prior_params_dispersion[[i]])

#---- Create entries in the data block of the .stan file:
standata <- c(standata_preprocessed,
              nlist(
                ###
                prior_dist_nb_dispersion  = prior_dist_nb_dispersion,
                prior_mean_nb_dispersion  = prior_mean_nb_dispersion,
                prior_scale_nb_dispersion = prior_scale_nb_dispersion,
                prior_df_nb_dispersion    = prior_df_nb_dispersion,
                prior_shape_nb_dispersion = prior_shape_nb_dispersion,
                prior_rate_nb_dispersion  = prior_rate_nb_dispersion
              ))

#---- List of parameters that will be monitored:
parameters <- c("mu_b0c",
                "sigma_b0c",
                "phi",
                "b0ic",
                "beta",
                "lp",
                "bW_mat",
                "log_lik",
                "deviance")

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
rstan_options(auto_write = TRUE)
options(mc.cores = nChains)

time_start_nuts1 <- Sys.time()
nuts_fit_1 <- rstan::sampling(m1,
                          data    = standata,
                          chains  = nChains,
                          warmup  = nBurnin,
                          iter    = nIter,
                          seed    = 1,
                          control = list(max_treedepth = max_treedepth,
                                         adapt_delta   = adapt_delta
                          ),
                          show_messages   = FALSE)
time_end_nuts1 <- Sys.time()
duration_nuts1 <- time_end_nuts1 - time_start_nuts1

time_run <- get_elapsed_time(nuts_fit_1)
apply(time_run, 1, sum)

###########################################################################
#
# Store HMC output and MCMC summaries
#
###########################################################################
bridgeBase <- bridgesampling::bridge_sampler(samples = nuts_fit_1)

save(standata,
     nuts_fit_1,
     duration_nuts1,
     key_CHIKV_filtered,
     main_dataset_CHIKV,
     O_t,
     bridgeBase,
     file = store_model_out)

load(file = store_model_out)

nuts_fit_1_summary <- summary(nuts_fit_1)$summary

print(bridgeBase)
bridgeBase$logml
bridgesampling::error_measures(bridgeBase)$percentage
