###########################################################################
#
# Setup
#
###########################################################################

#---- Session info: change locale to English for correct axis
#     label display in the graphs:
# Source: https://stackoverflow.com/questions/15438429/axis-labels-are-not-plotted-in-english
sessionInfo()

#---- Libraries:
lib <- c("ggplot2",
         "tidyverse",
         "boot",
         "readr",
         "dplyr",
         "tidyr",
         "rvest",
         "vroom",
         "rstan",
         "bayesplot",
         "gridExtra",
         #"mgcv",
         "Bernadette",
         #"readxl",
         "compiler",
         "extraDistr",
         "patchwork",
         "latex2exp"
)
lapply(lib, require, character.only = TRUE)

#---- Session info:
sessionInfo()

.libPaths()

###########################################################################
#
# Paths
#
###########################################################################
main_path       <- "C://Users//lbour//OneDrive//Desktop1//xgp_project//GP multitype epidemic England//4_Stan_Outputs//"
store_model_out <- paste0(main_path, "SEEIIR_renewal_multitype_XGBM", ".RData") #SEEIIR_renewal_multitype_XGBM_v5C

#---- Set system locale to English:
Sys.setlocale("LC_ALL", "English")

#---- Data engineering; daily time series of new cases and deaths for a given period:
# source(paste0(path_to_source, "2_Data_Engineering", ".R"))

###########################################################################
#
# Load libraries:
#
###########################################################################
load(file = store_model_out)
# nuts_fit_1  <- fit_ENG_xGP_3cp
posts_mcmc  <- rstan::extract(nuts_fit_1) # Posterior draws 
summaryMCMC <- summary(nuts_fit_1)$summary

summaryMCMC[which(rownames(summaryMCMC) == "Deviance"),]
mean(posts_mcmc$Deviance)

###########################################################################
#
# Estimation of intra-class coefficient for the xGBM model
#
###########################################################################
xgbm_rho <- posts_mcmc[["sigma_x"]]^2/(posts_mcmc[["sigma_x"]]^2 + posts_mcmc[["sigma_mu"]]^2)

round(c(mean(xgbm_rho), quantile(xgbm_rho, probs = c(0.025, 0.975))), 3)

###########################################################################
#
# Estimation of intra-class coefficient for the mxGBM model
#
###########################################################################
dt_store <- data.frame(mean = rep(0, dim(posts_mcmc[["sigma_x"]])[2]),
                       low  = rep(0, dim(posts_mcmc[["sigma_x"]])[2]),
                       high = rep(0, dim(posts_mcmc[["sigma_x"]])[2]))
for (i in 1:dim(posts_mcmc[["sigma_x"]])[2]){
  
  tmp_mcmc <- posts_mcmc[["sigma_x"]][,i]^2/(posts_mcmc[["sigma_x"]][,i]^2 + posts_mcmc[["sigma_mu"]]^2)
  
  dt_store[i,1] <- mean(tmp_mcmc)
  dt_store[i,2] <- quantile(tmp_mcmc, probs = 0.025)
  dt_store[i,3] <- quantile(tmp_mcmc, probs = 0.975)
}

round(dt_store, 3)

###########################################################################
#
# Estimation of intra-class coefficient for the xGP model
#
###########################################################################
xgp_rho <- posts_mcmc[["sigma_age_group"]]^2/(posts_mcmc[["sigma_age_group"]]^2 + posts_mcmc[["sigma_country"]]^2)

round(c(mean(xgp_rho), quantile(xgp_rho, probs = c(0.025, 0.975))), 3)

###########################################################################
#
# Estimation of intra-class coefficient for the mxGP model
#
###########################################################################
dt_store <- data.frame(mean = rep(0, dim(posts_mcmc[["sigma_age_group"]])[2]),
                       low  = rep(0, dim(posts_mcmc[["sigma_age_group"]])[2]),
                       high = rep(0, dim(posts_mcmc[["sigma_age_group"]])[2]))

for (i in 1:dim(posts_mcmc[["sigma_age_group"]])[2]){
  
  tmp_mcmc <- posts_mcmc[["sigma_age_group"]][,i]^2/(posts_mcmc[["sigma_age_group"]][,i]^2 + posts_mcmc[["sigma_country"]]^2)
  
  dt_store[i,1] <- mean(tmp_mcmc)
  dt_store[i,2] <- quantile(tmp_mcmc, probs = 0.025)
  dt_store[i,3] <- quantile(tmp_mcmc, probs = 0.975)
}

round(dt_store, 3)

###########################################################################
#
# Traceplots of hyper-parameters
#
###########################################################################
posterior_1 <- as.array(nuts_fit_1)
colnames(posterior_1[1,,])

color_scheme_set("viridis")

model <- "xGBM"

if (model %in% c("xGBM", "xGP") ){
  
  sigma_mu_raw <- as.data.frame(posterior_1[, , 287]) # sigma_mu
  sigma_x_raw  <- as.data.frame(posterior_1[, , 286]) # sigma_x
  phi_raw      <- as.data.frame(posterior_1[, , 5]) # phi
  
  colnames(sigma_x_raw) <- colnames(sigma_mu_raw) <- colnames(phi_raw) <- as.character(1:ncol(sigma_x_raw))
  sigma_x_raw$Iteration  <- 1:nrow(sigma_x_raw)
  sigma_mu_raw$Iteration <- 1:nrow(sigma_mu_raw)
  phi_raw$Iteration      <- 1:nrow(phi_raw)
  
  sigma_x_raw$Variable  <- rep("sigma_x", nrow(sigma_x_raw))
  sigma_mu_raw$Variable <- rep("sigma_mu", nrow(sigma_mu_raw))
  phi_raw$Variable      <- rep("phiD", nrow(phi_raw))
  
  sigmas_merged <- rbind(sigma_x_raw, sigma_mu_raw, phi_raw)
  
  sigmas_merged_long <- 
    pivot_longer(sigmas_merged, 
                 -c(Variable, Iteration), 
                 values_to = "value", 
                 names_to  = "Chain")
  
  sigmas_merged_long <-
    sigmas_merged_long %>%
    mutate(Parameter = factor(Variable,
                              levels = c("sigma_x", "sigma_mu", "phiD"),
                              labels = c(bquote(sigma[x]),
                                         bquote(sigma[mu]),
                                         bquote(phi[.("D")]) ) ) )
  sigmas_merged_long$Variable<- NULL
  sigmas_merged_long$Chain <- as.numeric(sigmas_merged_long$Chain)
  
} else if (model %in% c("mxGBM", "mxGP") ){
  
  sigma_mu_raw  <- as.data.frame(posterior_1[, , colnames(posterior_1[1,,])[colnames(posterior_1[1,,]) == "sigma_mu"]]) 
  sigma_x_raw_1 <- as.data.frame(posterior_1[, , colnames(posterior_1[1,,])[colnames(posterior_1[1,,]) == "sigma_x[1]"]])
  sigma_x_raw_2 <- as.data.frame(posterior_1[, , colnames(posterior_1[1,,])[colnames(posterior_1[1,,]) == "sigma_x[2]"]])
  sigma_x_raw_3 <- as.data.frame(posterior_1[, , colnames(posterior_1[1,,])[colnames(posterior_1[1,,]) == "sigma_x[3]"]])
  sigma_x_raw_4 <- as.data.frame(posterior_1[, , colnames(posterior_1[1,,])[colnames(posterior_1[1,,]) == "sigma_x[3]"]])
  phi_raw       <- as.data.frame(posterior_1[, , colnames(posterior_1[1,,])[colnames(posterior_1[1,,]) == "phiD"]]) 
  
  colnames(sigma_mu_raw) <- 
    colnames(sigma_x_raw_1) <- 
    colnames(sigma_x_raw_2) <- 
    colnames(sigma_x_raw_3) <- 
    colnames(sigma_x_raw_4) <- 
    colnames(phi_raw) <- as.character(1:ncol(phi_raw))
  
  sigma_mu_raw$Iteration   <- 1:nrow(phi_raw)
  sigma_x_raw_1$Iteration  <- 1:nrow(phi_raw)
  sigma_x_raw_2$Iteration  <- 1:nrow(phi_raw)
  sigma_x_raw_3$Iteration  <- 1:nrow(phi_raw)
  sigma_x_raw_4$Iteration  <- 1:nrow(phi_raw)
  phi_raw$Iteration        <- 1:nrow(phi_raw)
  
  sigma_mu_raw$Variable  <- rep("sigma_mu", nrow(sigma_mu_raw))
  sigma_x_raw_1$Variable <- rep("sigma_x_1", nrow(phi_raw))
  sigma_x_raw_2$Variable <- rep("sigma_x_2", nrow(phi_raw))
  sigma_x_raw_3$Variable <- rep("sigma_x_3", nrow(phi_raw))
  sigma_x_raw_4$Variable <- rep("sigma_x_4", nrow(phi_raw))
  phi_raw$Variable       <- rep("phi", nrow(phi_raw))
  
  sigmas_merged <- rbind(sigma_mu_raw, 
                         sigma_x_raw_1, 
                         sigma_x_raw_2, 
                         sigma_x_raw_3, 
                         sigma_x_raw_4, 
                         phi_raw)
  
  sigmas_merged_long <- 
    pivot_longer(sigmas_merged, 
                 -c(Variable, Iteration), 
                 values_to = "value", 
                 names_to  = "Chain")
  
  sigmas_merged_long <-
    sigmas_merged_long %>%
    mutate(Parameter = factor(Variable,
                              levels = c("sigma_mu", 
                                         "sigma_x_1", 
                                         "sigma_x_2", 
                                         "sigma_x_3",
                                         "sigma_x_4",
                                         "phi"),
                              labels = c(bquote(sigma[mu]),
                                         expression(sigma[1]),
                                         expression(sigma[2]),
                                         expression(sigma[3]),
                                         expression(sigma[4]),
                                         bquote(phi[.("D")]) ) ) )
  sigmas_merged_long$Variable <- NULL
  sigmas_merged_long$Chain    <- as.numeric(sigmas_merged_long$Chain)

}

# Save as pdf: 8x8:
# Figure caption: 
# COVID-19 epidemic in England. Trace plot of the four chains of the xGBM model hyper-parameters.
# Figure name: COVID19_mxGBM_trace
ggplot(sigmas_merged_long, 
       aes(x = Iteration, 
           y = value, 
           colour = as.factor(Chain))) +
  facet_wrap(~Parameter, 
             ncol = 1, 
             scales = "free", 
             labeller = label_parsed)+
  geom_line(alpha = 0.7) + 
  scale_colour_discrete(name = "Chain") +
  labs(x = "Iteration",
       y = "Value") +
  theme_bw() +
  theme(panel.spacing    = unit(0.4,"cm"),
        axis.text.x      = element_text(size = 12, angle = 0, hjust = 1),
        axis.text.y      = element_text(size = 12),
        axis.title.x     = element_text(size = 16, face = "bold"),
        axis.title.y     = element_text(size = 16, face = "bold"),
        legend.text      = element_text(size = 16),
        legend.title     = element_text(size = 16),
        strip.text.x     = element_text(size = 14),
        legend.position  = "bottom",
        legend.box       = "vertical",
        legend.margin    = margin())

###########################################################################
#
# Load the posterior draws:
#
###########################################################################

if( "y_deaths" %in% names(cov_data)) {
  age_bands <- colnames(cov_data$y_deaths)
} else { 
  age_bands <- colnames(cov_data$y_data)
}# End if

start_date <- "2020-03-02"
end_date   <- "2020-09-27"
dates      <- seq(as.Date(start_date), as.Date(end_date), by = "days")

##########################################################################
#
# mxGBM goodness of fit:
#
###########################################################################

#---- Posterior random draws of age-specific mortality counts:
# model_out = nuts_fit_1
# cov_data  = stan_data
# reported_data = data_rt_wide
deaths_random_draws <- function(cov_data,
                                model_out,
                                start_date,
                                end_date,
                                model = c("GBM", "GP")){
  
  set.seed(1)
  
  #---- Checks:
  if(class(model_out)[1] != "stanfit") stop("Provide an object of class 'stanfit' using rstan::sampling()")
  
  #---- Posterior draws:
  posterior_draws <- rstan::extract(model_out)
  
  if (model == "GBM") {
    E_deathsByAge   <- posterior_draws[["E_deathsAge"]]
    E_deaths        <- posterior_draws[["E_deaths"]]
    phi             <- posterior_draws[["phiD"]]
    
  } else {
    E_deathsByAge   <- posterior_draws[["E_deaths"]]
    E_deathsByAge   <- aperm(E_deathsByAge, perm = c(1, 3, 2))
    E_deaths        <- posterior_draws[["total_E_deaths"]]
    phi             <- posterior_draws[["phi"]]
    
  }
    
  mcmc_length     <- dim(E_deathsByAge)[1]
  ts_length       <- dim(E_deathsByAge)[2]
  age_grps        <- dim(E_deathsByAge)[3]
  dates           <- seq(as.Date(start_date), as.Date(end_date), by = "days")
    
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

# nuts_fit = nuts_fit_1
posterior_age_specific_counts <- function(nuts_fit,
                                          cov_data,
                                          start_date,
                                          end_date,
                                          model = c("GBM", "GP")){
  
  #---- Checks:
  if(class(nuts_fit)[1] != "stanfit") stop("Provide an object of class 'stanfit' using rstan::sampling()")
  
  #---- Posterior draws:
  posterior_draws <- rstan::extract(nuts_fit)
  
  #---- Placeholder for output tables:
  data_deaths_cols      <- c("Date", "Group", "New_Deaths", "median", "low", "high")
  data_deaths           <- data.frame(matrix(ncol = length(data_deaths_cols), nrow = 0))
  colnames(data_deaths) <- data_deaths_cols
  
  if( "y_deaths" %in% names(cov_data)) {
    group_id <- colnames(cov_data$y_deaths)
  } else { 
    group_id <- colnames(cov_data$y_data)
  }# End if
  
  # Loop over each age group:
  for (i in 1:cov_data$A){
    
    # word <- c("E_deathsByAge", "E_deaths")
    # 
    # if( any( word %in% names(posterior_draws)) ){
    #   
    #   fit_deaths <- posterior_draws$E_deathsByAge[,,i] #E_deathsByAge[,,i]  #
    #   
    # } else {
    #   fit_deaths <- posterior_draws$E_deathsAge[,,i] #E_deathsByAge[,,i]  #
    # }
    
    if (model == "GBM") {
      fit_deaths <- posterior_draws$E_deathsAge[,,i] #E_deathsByAge[,,i]  #
    } else {
      tmp <- posterior_draws[["E_deaths"]]
      tmp <- aperm(tmp, perm = c(1, 3, 2))
      fit_deaths <- tmp[,,i]
    }
    
    dates <- seq(as.Date(start_date), as.Date(end_date), by = "days")
    
    dt_deaths_age_grp  <- data.frame(Date  = dates,
                                     Group = rep(group_id[i], length(dates) ))

    filter_out_deaths_cols <- c("Index", "Right", "Week_ID", "New_Deaths")
    

    dt_mortality_analysis_period_melt <- 
      cbind(data.frame(Date = dates), cov_data$y_data)
    dt_mortality_analysis_period_melt <- dt_mortality_analysis_period_melt[c("Date", group_id[i])]
    
    dt_mortality_analysis_period_melt <- 
      dt_mortality_analysis_period_melt %>%
      pivot_longer(!Date, names_to = "Group", values_to = "New_Deaths")
    
    dt_deaths_age_grp <- 
      dt_deaths_age_grp %>%
      dplyr::left_join(dt_mortality_analysis_period_melt,
                       by = c("Date"  = "Date",
                              "Group" = "Group"))
    
    # Add quantiles from the model outputs:
   
    dt_deaths_age_grp$median <- apply(fit_deaths, 2, median)
    dt_deaths_age_grp$low    <- apply(fit_deaths, 2, quantile, probs = c(0.025))
    dt_deaths_age_grp$high   <- apply(fit_deaths, 2, quantile, probs = c(0.975))
    dt_deaths_age_grp$low25  <- apply(fit_deaths, 2, quantile, probs = c(0.25))
    dt_deaths_age_grp$high75 <- apply(fit_deaths, 2, quantile, probs = c(0.75))
    
    data_deaths <- rbind(data_deaths, dt_deaths_age_grp)
    
  }# End for
  
  output <- list(Deaths = data_deaths)
  
  return(output)
  
}# End function

posterior_aggregated_counts <- function(cov_data,
                                        nuts_fit,
                                        start_date,
                                        end_date,
                                        model = c("GBM", "GP")){
  
  #---- Checks:
  if(class(nuts_fit)[1] != "stanfit") stop("Provide an object of class 'stanfit' using rstan::sampling()")
  
  dates <- seq(as.Date(start_date), as.Date(end_date), by = "days")
  
  #---- Posterior draws:
  posterior_draws <- rstan::extract(nuts_fit)
  
  #---- Deaths (posterior median, 50% and 95% credible intervals):
  if (model == "GBM") fit_dead <- posterior_draws$E_deaths else fit_dead <- posterior_draws$total_E_deaths
  
  median_fit_dead <- apply(fit_dead, 2, median)
  low_fit_dead    <- apply(fit_dead, 2, quantile, probs = c(0.025)) 
  high_fit_dead   <- apply(fit_dead, 2, quantile, probs = c(0.975))
  low25_fit_dead  <- apply(fit_dead, 2, quantile, probs = c(0.25))
  high75_fit_dead <- apply(fit_dead, 2, quantile, probs = c(0.75))
  
  deaths_output        <- data.frame(Date       = dates,
                                     New_Deaths = apply(cov_data$y_data,1, sum))
  deaths_output$median <- median_fit_dead
  deaths_output$low    <- low_fit_dead
  deaths_output$high   <- high_fit_dead
  deaths_output$low25  <- low25_fit_dead
  deaths_output$high75 <- high75_fit_dead
 
  output <- list(Deaths = deaths_output)
  
  return(output)
  
}# End function

aggregated_fits        <- posterior_aggregated_counts(cov_data,
                                                      nuts_fit = nuts_fit_1,
                                                      start_date,
                                                      end_date,
                                                      model = "GBM")
aggregated_deaths_data <- aggregated_fits$Deaths

#---- Generate random samples of age-stratified deaths from the model:
deaths_rng <- deaths_random_draws(cov_data,
                                  model_out = nuts_fit_1,
                                  start_date,
                                  end_date,
                                  model = "GBM")

aggr_deaths_rng        <- deaths_rng$aggregated_deaths

aggregated_deaths_data <- 
  aggregated_deaths_data %>% 
  left_join(aggr_deaths_rng, by = c("Date"))

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
  # scale_y_continuous(
  #   limits = c(0,   max(aggregated_deaths_data$high)*1.4),    
  #   breaks = seq(0, max(aggregated_deaths_data$high)*1.4, 100) 
  # ) +
  scale_fill_manual(values = c('Reported' = "black",
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

###########################################################################
#
# Goodness of fit to new age-stratified deaths:
#
###########################################################################
age_specific_fits  <- posterior_age_specific_counts(nuts_fit = nuts_fit_1,
                                                    cov_data,
                                                    start_date,
                                                    end_date,
                                                    model = "GBM")
age_specific_deaths_data <- age_specific_fits$Deaths
age_deaths_rng           <- deaths_rng$age_specific_deaths
age_specific_deaths_data <- 
  age_specific_deaths_data %>% 
  left_join(age_deaths_rng,  by = c("Date", "Group"))

fit_deaths_age <- 
  ggplot(age_specific_deaths_data,
         aes(Date       = Date,
             New_Deaths = New_Deaths)) +
  facet_wrap(. ~ Group, 
             scales = "free_y", 
             ncol   = 1,
             strip.position = "top") +
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
  scale_fill_manual(values = c('Reported' = "black",
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

#---- Model fit to observed deaths:
(fit_deaths_age | fit_deaths_all) + plot_annotation(tag_levels =  "A")

# Preserve transparency when saving to .eps
# Source: https://www.sthda.com/english/wiki/saving-high-resolution-ggplots-how-to-preserve-semi-transparency
grDevices::cairo_ps(filename = "C://Users//lbour//OneDrive//Desktop1//xgp_project//GP multitype epidemic England//6_Model_plots//COVID19_xGBM_model_fit_deaths.eps",
                    width     = 12, 
                    height    = 10, 
                    pointsize = 12,
                    fallback_resolution = 300)
print( (fit_deaths_age | fit_deaths_all) + plot_annotation(tag_levels =  "A") )
dev.off()

###########################################################################
#
# Age-specific probability of infection (Rho_t) trajectories vs population
#
###########################################################################
data_rho_cols     <- c("Date", "Group", "median", "low", "high", "low25", "high75")
data_rho           <- data.frame(matrix(ncol = length(data_rho_cols), nrow = 0))
colnames(data_rho) <- data_rho_cols

for (i in 1:cov_data$A){

  fit_rho  <- posts_mcmc$rho_weekly[,,i] #-c(1:rm_dates)

  dt_rho_country  <- data.frame(Date   = dates[seq(1, length(dates), cov_data$gbm_no_days)], # Keep every 3rd element, due to having changes every 3 days
                                Group  = rep(age_bands[i], dim(fit_rho)[2]) )

  # Add quantiles from the model outputs:
  dt_rho_country$median <- apply(fit_rho, 2, median)
  dt_rho_country$low    <- apply(fit_rho, 2, quantile, probs = c(0.025))
  dt_rho_country$high   <- apply(fit_rho, 2, quantile, probs = c(0.975))
  dt_rho_country$low25  <- apply(fit_rho, 2, quantile, probs = c(0.25))
  dt_rho_country$high75 <- apply(fit_rho, 2, quantile, probs = c(0.75))

  data_rho  <- rbind(data_rho, dt_rho_country)

}# End for
rm(i, dt_rho_country)

datalogit_rho_cols     <- c("Date", "Group", "median", "low", "high", "low25", "high75")
datalogit_rho           <- data.frame(matrix(ncol = length(datalogit_rho_cols), nrow = 0))
colnames(datalogit_rho) <- datalogit_rho_cols

for (i in 1:cov_data$A){
  
  fitlogit_rho  <- boot::logit(posts_mcmc$rho_weekly[,,i]) #-c(1:rm_dates)
  
  dtlogit_rho_country  <- data.frame(Date   = dates[seq(1, length(dates), cov_data$gbm_no_days)], # Keep every 3rd element, due to having changes every 3 days
                                Group  = rep(age_bands[i], dim(fitlogit_rho)[2]) )
  
  # Add quantiles from the model outputs:
  dtlogit_rho_country$median <- apply(fitlogit_rho, 2, median)
  dtlogit_rho_country$low    <- apply(fitlogit_rho, 2, quantile, probs = c(0.025))
  dtlogit_rho_country$high   <- apply(fitlogit_rho, 2, quantile, probs = c(0.975))
  dtlogit_rho_country$low25  <- apply(fitlogit_rho, 2, quantile, probs = c(0.25))
  dtlogit_rho_country$high75 <- apply(fitlogit_rho, 2, quantile, probs = c(0.75))
  
  datalogit_rho  <- rbind(datalogit_rho, dtlogit_rho_country)
  
}# End for

#---- Check the transformed model outputs:
table(data_rho$Group)

#---- Remove redundant objects:
rm(i, dt_rho_country)

#---- Posterior median infection probability per group:
rho_plot <- 
ggplot(data_rho,
       aes(x      = Date,
           y      = median,
           colour = Group)) +
  facet_wrap(. ~ Group, 
             scales = "free_y", 
             ncol   = 1,
             strip.position = "top") +
  geom_line(aes(x     = Date,
                y     = median,
                color = "Median"),
            size = 1.3) +
  geom_ribbon(aes(x    = Date,
                  ymin = low25,
                  ymax = high75,
                  fill = "50% CrI")) +
  geom_ribbon(aes(x    = Date,
                  ymin = low,
                  ymax = high,
                  fill = "95% CrI"),
              alpha = 0.5) +
  labs(x = "Epidemiological Date",
       y = "Transmissibility") +
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

logit_rho_plot <- 
  ggplot(datalogit_rho,
         aes(x      = Date,
             y      = median,
             colour = Group)) +
  facet_wrap(. ~ Group, 
             scales = "free_y", 
             ncol   = 1,
             strip.position = "top") +
  geom_line(aes(x     = Date,
                y     = median,
                color = "Median"),
            size = 1.3) +
  geom_ribbon(aes(x    = Date,
                  ymin = low25,
                  ymax = high75,
                  fill = "50% CrI")) +
  geom_ribbon(aes(x    = Date,
                  ymin = low,
                  ymax = high,
                  fill = "95% CrI"),
              alpha = 0.5) +
  labs(x = "Epidemiological Date",
       y = "logit(Transmissibility)") +
  scale_x_date(date_breaks       = "1 month", 
               date_minor_breaks = "1 month") + 
  scale_y_continuous(
    limits = c(round(min(datalogit_rho$low),0),   round(max(datalogit_rho$high),0) )#,    
    #breaks = seq(round(min(data_rho$low),0), round(max(data_rho$high),0), 4) 
  ) +
  scale_fill_manual(values = c("50% CrI" = "gray40",
                               "95% CrI" = "gray70"
  )) +
  scale_colour_manual(name   = '',
                      values = c('Median' = "black")) +
  geom_hline(yintercept = 0.0,
             colour   = "black", 
             linetype = "dashed") +
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

# Source: https://github.com/tidyverse/ggplot2/issues/5619
# PDF = 12 x 10
library(patchwork)
combined <-
  (logit_rho_plot +
     rho_plot ) + 
  plot_layout(guides = 'collect') &
  theme(legend.position = 'bottom', legend.direction = 'horizontal')

# Preserve transparency when saving to .eps
# Source: https://www.sthda.com/english/wiki/saving-high-resolution-ggplots-how-to-preserve-semi-transparency
grDevices::cairo_ps(filename = "C://Users//lbour//OneDrive//Desktop1//xgp_project//GP multitype epidemic England//6_Model_plots//COVID19_xGBM_transmissibility.eps",
                    width     = 12, 
                    height    = 10, 
                    pointsize = 12,
                    fallback_resolution = 300)
print(combined + plot_annotation( tag_levels = list(c("A", "B", "C") ))  )
dev.off()

###########################################################################
#
# Effective reproduction number - Next generation matrix approach 
# SE2I2R transmission model
#
###########################################################################

# # NOTE: K is not a symmetric matrix:
# # https://github.com/CBDRH/covoid/blob/master/R/1-sir-cm-model.R
# R <- max( Re( eigen(K)$values ) )
start_date <- "2020-03-02"
end_date   <- "2020-09-27"
dates     <- seq(as.Date(start_date), as.Date(end_date), by = "1 days")

effective_rt_seeiir <- function(cov_data,
                                model_out,
                                dates,
                                dI = 6.5, # = 6.5 same as the mean of the generation-time interval distribution
                                progress_bar = FALSE){
  
  # model_out <- nuts_fit_1
  # rm(aggr_age, model_out, nuts_fit_1)
  
  `%nin%` <- Negate(`%in%`)
  
  #---- Checks:
  if(class(model_out)[1] != "stanfit") stop("Provide an object of class 'stanfit' using rstan::sampling()")
  
  #---- Posterior draws:
  posterior_draws <- rstan::extract(model_out)
  pop_diag        <- 1/(cov_data$n_pop)
  
  #---- Functions:
  # Source: https://stackoverflow.com/questions/28153822/how-to-apply-a-function-to-a-list-of-matrices-elementwise
  eigen_mat  <- function(mat) max( Re( eigen(mat)$values ) )
  
  #---- Assistant matrices:
  age_grps             <- cov_data$A
  ones_mat             <- matrix(1L, nrow = age_grps, ncol = age_grps)  
  reciprocal_age_distr <- matrix(rep(pop_diag, age_grps), ncol  = age_grps, nrow  = age_grps, byrow = TRUE)
  
  #---- Output storage:
  # The R matrix must have the same dimensions as the matrix posterior_draws$beta_N:
  
  beta_draws   <- posterior_draws$rho_daily
  chain_length <- nrow(beta_draws)
  ts_length    <- dim(beta_draws)[2]
  R_mat        <- matrix(0L, nrow = chain_length, ncol = ts_length)
  
  # Step 1 - Calculate Q^{-1} at iteration i:  
  Q_inverse <- solve (diag(rep(1/dI, age_grps)) )
  #Q_inverse <- dI*ones_mat
  
  if(progress_bar == TRUE) pb = txtProgressBar(min = 0, max = chain_length, initial = 0, style = 3) 
  
  #---- Calculation of the effective at iteration i and time point j:
  for (i in 1:chain_length) {
    
    for (j in 1:ts_length){
      
      # Step 2 - Create the non-zero part of the B_t matrix: 
      # NOTE: Remember to perform element-wise multiplication for this part:
      
      # Common GBM across groups:
      if( is.na(dim(beta_draws)[3]) ){
        
        B_inside_tmp <- 
          #dI*
          beta_draws[i,j] *
          cov_data$contact_matrix *
          matrix(rep( posterior_draws$Susceptibles[i,,][j,], age_grps),
                 ncol  = age_grps, 
                 nrow  = age_grps, 
                 byrow = FALSE) *
          reciprocal_age_distr
        
        # Multiple GBMs: 
      } else if ( !is.na(dim(beta_draws)[3]) ){
        
        B_inside_tmp <- 
          #dI*
          matrix( rep(beta_draws[i,,][j,], age_grps),
                  ncol  = age_grps, 
                  nrow  = age_grps, 
                  byrow = FALSE) *
          cov_data$contact_matrix *
          matrix(rep( posterior_draws$Susceptibles[i,,][j,], age_grps),
                 ncol  = age_grps, 
                 nrow  = age_grps, 
                 byrow = FALSE) *
          reciprocal_age_distr
        
      }# End if  
      
      # Step 3 - Create the B_t matrix:
      # B_tmp <- rbind(cbind(zero_mat, B_inside_tmp),
      #                cbind(zero_mat, zero_mat) ) 
      
      B_tmp <- B_inside_tmp
      
      # Step 4 calculate the B_t %*% Q^{-1} matrix:
      BQinv_tmp <- B_tmp %*% Q_inverse
      
      #BQinv_tmp <- B_inside_tmp
      
      # Step 5 - calculate the spectral radius of B_t %*% Q^{-1} at iteration i and time point j:
      R_mat[i,j] <- eigen_mat(BQinv_tmp)
      
    }# End for
    
    #---- Cleanup:
    #rm(B_inside_tmp, B_tmp, BQinv_tmp)
    
    if(progress_bar == TRUE) setTxtProgressBar(pb,i)
    
  }# End for
  
  data_repnumber                <- data.frame(Date  = dates)
  data_repnumber$eff_rt_median  <- apply(R_mat, 2, median)
  data_repnumber$eff_rt_low0025 <- apply(R_mat, 2, quantile, probs = c(0.025)) # c(0.025)
  data_repnumber$eff_rt_low25   <- apply(R_mat, 2, quantile, probs = c(0.25))  # c(0.025)
  data_repnumber$eff_rt_high75  <- apply(R_mat, 2, quantile, probs = c(0.75))  # c(0.975)
  data_repnumber$eff_rt_high975 <- apply(R_mat, 2, quantile, probs = c(0.975)) # c(0.975)
  
  #---- Export:
  return(data_repnumber)
  
}# End function

#---- Compilation of the function will not be required when I include the function in the package:
effective_rt_seeiir <- compiler::cmpfun(effective_rt_seeiir)
effRt <- effective_rt_seeiir(cov_data,
                             model_out = nuts_fit_1,
                             dates,
                             dI = 6.5,                                
                             progress_bar = TRUE)

save(effRt,
     file = "C://Users//lbour//OneDrive//Desktop1//xgp_project//GP multitype epidemic England//6_Model_plots//effRt_xBM_NGM.RData")

load(file = "C://Users//lbour//OneDrive//Desktop1//xgp_project//GP multitype epidemic England//6_Model_plots//effRt_xBM_NGM.RData")

# Keep records from every 3 days:
effRt <- effRt %>% as.data.frame %>% slice(which(row_number() %% 3 == 1))

plotRtNGM <- 
ggplot(effRt,
       aes(Date   = Date,
           Median = eff_rt_median)) +
  geom_line(aes(x     = Date,
                y     = eff_rt_median,
                color = "Median"),
            size = 1.3) +
  geom_ribbon(aes(x    = Date,
                  ymin = eff_rt_low25,
                  ymax = eff_rt_high75,
                  fill = "50% CrI"),
              alpha = 0.5) +
  geom_ribbon(aes(x    = Date,
                  ymin = eff_rt_low0025,
                  ymax = eff_rt_high975,
                  fill = "95% CrI"),
              alpha = 0.5) +
  labs(x = "Epidemiological Date",
       y = "Effective reproduction number (NGM)") +
  scale_x_date(date_breaks       = "2 weeks",
               date_minor_breaks = "2 week") +
  #geom_hline(yintercept = 1, color = "black") +
  scale_y_continuous(
  #   limits = c( min(rt_ngmplot_data$low25)*0.8,   max(rt_ngmplot_data$high75)*1.1),    # max(death_plot_data$high)*1.2
    breaks = seq( 0.0, 45, 1) #max(death_plot_data$high)*1.2
  ) +
  scale_fill_manual(values = c("50% CrI" = "gray50",
                               "95% CrI" = "gray80"
  )) +
  scale_colour_manual(name   = '',
                      values = c('Median' = "black")) +
  geom_hline(yintercept = 1.0,
             colour   = "black", 
             linetype = "dashed") +
  theme_bw() +
  theme(panel.spacing    = unit(0.2,"cm"),
        axis.text.x      = element_text(size = 12, angle = 45, hjust = 1),
        axis.text.y      = element_text(size = 12),
        axis.title.x     = element_text(size = 16, face = "bold"),
        axis.title.y     = element_text(size = 16, face = "bold"),
        legend.position  = "bottom",
        legend.title     = element_blank(),
        legend.text      = element_text(size = 14),
        legend.margin    = margin(),
        strip.text.x     = element_text(size = 16)
  )

#---- Overall R_t from renewal equation:
datesRt                      <- seq(as.Date(start_date), as.Date(end_date), by = "1 days")
dataRtRenewal                <- data.frame(Date  = datesRt)
dataRtRenewal$eff_rt_median  <- apply(posts_mcmc$Rt1, 2, median)
dataRtRenewal$eff_rt_low0025 <- apply(posts_mcmc$Rt1, 2, quantile, probs = c(0.025)) # c(0.025)
dataRtRenewal$eff_rt_low25   <- apply(posts_mcmc$Rt1, 2, quantile, probs = c(0.25))  # c(0.025)
dataRtRenewal$eff_rt_high75  <- apply(posts_mcmc$Rt1, 2, quantile, probs = c(0.75))  # c(0.975)
dataRtRenewal$eff_rt_high975 <- apply(posts_mcmc$Rt1, 2, quantile, probs = c(0.975)) # c(0.975)

# Keep records from every 3 days:
dataRtRenewal <- dataRtRenewal %>% as.data.frame %>% slice(which(row_number() %% 3 == 1))

plotRtRenewal <- 
ggplot(dataRtRenewal,
       aes(Date   = Date,
           Median = eff_rt_median)) +
  geom_line(aes(x     = Date,
                y     = eff_rt_median,
                color = "Median"),
            size = 1.3) +
  geom_ribbon(aes(x    = Date,
                  ymin = eff_rt_low25,
                  ymax = eff_rt_high75,
                  fill = "50% CrI"),
              alpha = 0.5) +
  geom_ribbon(aes(x    = Date,
                  ymin = eff_rt_low0025,
                  ymax = eff_rt_high975,
                  fill = "95% CrI"),
              alpha = 0.5) +
  labs(x = "Epidemiological Date",
       y = "Effective reproduction number (renewal equation)") +
  scale_x_date(date_breaks       = "2 weeks",
               date_minor_breaks = "2 week") +
  #geom_hline(yintercept = 1, color = "black") +
  scale_y_continuous(
    #   limits = c( min(rt_ngmplot_data$low25)*0.8,   max(rt_ngmplot_data$high75)*1.1),    # max(death_plot_data$high)*1.2
    breaks = seq( 0.0, 8, 1) #max(death_plot_data$high)*1.2
  ) +
  scale_fill_manual(values = c("50% CrI" = "gray50",
                               "95% CrI" = "gray80"
  )) +
  scale_colour_manual(name   = '',
                      values = c('Median' = "black")) +
  geom_hline(yintercept = 1.0,
             colour   = "black", 
             linetype = "dashed") +
  theme_bw() +
  theme(panel.spacing    = unit(0.2,"cm"),
        axis.text.x      = element_text(size = 12, angle = 45, hjust = 1),
        axis.text.y      = element_text(size = 12),
        axis.title.x     = element_text(size = 16, face = "bold"),
        axis.title.y     = element_text(size = 16, face = "bold"),
        legend.position  = "bottom",
        legend.title     = element_blank(),
        legend.text      = element_text(size = 14),
        legend.margin    = margin(),
        strip.text.x     = element_text(size = 16)
  )

combinedRt <-
  (plotRtNGM +
     plotRtRenewal ) + 
  plot_layout(guides = 'collect') &
  theme(legend.position = 'bottom', legend.direction = 'horizontal')
print(combinedRt + plot_annotation( tag_levels = list(c("A", "B") ))  )

# Preserve transparency when saving to .eps
# Source: https://www.sthda.com/english/wiki/saving-high-resolution-ggplots-how-to-preserve-semi-transparency
grDevices::cairo_ps(filename = "C://Users//lbour//OneDrive//Desktop1//xgp_project//GP multitype epidemic England//6_Model_plots//COVID19_xGBM_model_fit_deaths.eps",
                    width     = 12, 
                    height    = 10, 
                    pointsize = 12,
                    fallback_resolution = 300)
print( (fit_deaths_age | fit_deaths_all) + plot_annotation(tag_levels =  "A") )
dev.off()
