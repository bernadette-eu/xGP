###########################################################################
#
# Setup
#
###########################################################################

#---- Session info: change locale to English for correct axis
#     label display in the graphs:
sessionInfo()

#---- Libraries:
lib <- c("ggplot2",
         "tidyverse",
         "readr",
         "dplyr",
         "tidyr",
         "rvest",
         "vroom",
         "rstan",
         "bayesplot",
         "gridExtra",
         "Bernadette",
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
main_path       <- "...//Chikv//4_Stan_Outputs//"
store_model_out <- paste0(main_path, "CHIKV_FWI_2014_mxGP", ".RData")

#---- Set system locale to English:
Sys.setlocale("LC_ALL", "English")

###########################################################################
#
# Load libraries:
#
###########################################################################
load(file = store_model_out)

posts_mcmc         <- rstan::extract(nuts_fit_1) 
nuts_fit_1_summary <- summary(nuts_fit_1)$summary

options(scipen = 999)
round(nuts_fit_1_summary, 3)

###########################################################################
#
# Estimation of intra-class coefficient for the xGBM model
#
###########################################################################
xgbm_rho <- posts_mcmc[["sigma_mu"]]^2/(posts_mcmc[["sigma_x"]]^2 + posts_mcmc[["sigma_mu"]]^2)

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

  tmp_mcmc <- posts_mcmc[["sigma_mu"]]^2/(posts_mcmc[["sigma_x"]][,i]^2 + posts_mcmc[["sigma_mu"]]^2)

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
xgp_rho <- posts_mcmc[["sigma_common"]]^2/(posts_mcmc[["sigma_x"]]^2 + posts_mcmc[["sigma_common"]]^2)

round(c(mean(xgp_rho), quantile(xgp_rho, probs = c(0.025, 0.975))), 3)

###########################################################################
#
# Estimation of intra-class coefficient for the mxGP model
#
###########################################################################
dt_store <- data.frame(mean = rep(0, dim(posts_mcmc[["sigma"]])[2]),
                       low  = rep(0, dim(posts_mcmc[["sigma"]])[2]),
                       high = rep(0, dim(posts_mcmc[["sigma"]])[2]))
for (i in 1:dim(posts_mcmc[["sigma"]])[2]){

  tmp_mcmc <- posts_mcmc[["sigma_common"]]^2/(posts_mcmc[["sigma"]][,i]^2 + posts_mcmc[["sigma_common"]]^2)

  dt_store[i,1] <- mean(tmp_mcmc)
  dt_store[i,2] <- quantile(tmp_mcmc, probs = 0.025)
  dt_store[i,3] <- quantile(tmp_mcmc, probs = 0.975)
}

round(dt_store, 3)

#########################################################
# 
# Island-specific transmission coefficient
# 
#########################################################

beta_mat_cols      <- c("Date", "Group", "median", "low", "low25", "high75", "high")
beta_mat           <- data.frame(matrix(ncol = length(beta_mat_cols), nrow = 0))
colnames(beta_mat) <- beta_mat_cols

for (i in 1:standata$K){
  
  fit_beta          <- posts_mcmc$beta[,,i] 
  
  dt_effcr_group  <- data.frame(Date     = unique(main_dataset_CHIKV$DATE2),
                                Group    = rep(key_CHIKV_filtered$ISLAND[i], standata$W)
  )
  
  # Add quantiles from the model outputs:
  dt_effcr_group$median <- apply(fit_beta, 2, median)
  dt_effcr_group$low    <- apply(fit_beta, 2, quantile, probs = c(0.025))
  dt_effcr_group$low25  <- apply(fit_beta, 2, quantile, probs = c(0.25))
  dt_effcr_group$high75 <- apply(fit_beta, 2, quantile, probs = c(0.75))
  dt_effcr_group$high   <- apply(fit_beta, 2, quantile, probs = c(0.975))
  
  beta_mat  <- rbind(beta_mat, dt_effcr_group)
  
}# End for

beta_mat$Group_f <- factor(beta_mat$Group, levels = key_CHIKV_filtered$ISLAND)

#---- Plot of posterior estimates:
posterior_beta <- 
  ggplot(beta_mat) +
  facet_wrap(. ~ Group, 
             scales = "free_y", 
             ncol   = 1,
             strip.position = "top") +
  geom_line(aes(x     = Date,
                y     = median,
                color = "Median"),
            size  = 1.0) +
  geom_ribbon(aes(x    = Date,
                  ymin = low25,
                  ymax = high75,
                  fill = "50% CrI"),
              alpha = 0.5) +
  geom_ribbon(aes(x    = Date,
                  ymin = low,
                  ymax = high,
                  fill = "95% CrI"),
              alpha = 0.5) +
  labs(x = "Date",
       y = "Transmission") +
  scale_x_date(date_breaks       = "2 weeks",
               date_minor_breaks = "2 weeks") + 
  scale_fill_manual(values = c("50% CrI" = "gray30",
                               "95% CrI" = "gray70"
  )) +
  scale_colour_manual(name   = '',
                      values = c('Median' = "black")) +
  theme_bw() +
  theme(strip.placement  = "outside",
        axis.text.x      = element_text(size = 12, angle = 45, hjust = 1),
        axis.text.y      = element_text(size = 12),
        axis.title.x     = element_text(size = 16, face = "bold"),
        axis.title.y     = element_text(size = 16, face = "bold"),
        legend.text      = element_text(size = 14),
        strip.text.x     = element_text(size = 16),
        legend.position  = "bottom",
        legend.title     = element_blank() )
  
#########################################################
# 
# Expected cases - goodness of fit
# 
#########################################################

cases_random_draws <- function(dispersion_type = "linear",
                                cov_data,
                                model_out,
                                main_dataset,
                                key_mat
                                ){
  
  set.seed(1)
  
  #---- Checks:
  if(class(model_out)[1] != "stanfit") stop("Provide an object of class 'stanfit' using rstan::sampling()")
  
  #---- Posterior draws:
  posterior_draws <- rstan::extract(model_out)
  lp              <- posterior_draws[["lp"]]
  phi             <- posterior_draws[["phi"]]
  
  mcmc_length     <- dim(lp)[1]
  ts_length       <- dim(lp)[2]
  dates           <- main_dataset$DATE2
  
  death_draws           <- array(NA, c(mcmc_length, ts_length, cov_data$K))
  data_deaths_cols      <- c("Date",
                             "Group",
                             'median',
                             "rng_low",
                             "rng_low25",
                             "rng_high75",
                             "rng_high")
  age_rng_draws           <- data.frame(matrix(ncol = length(data_deaths_cols),
                                               nrow = 0))
  colnames(age_rng_draws) <- data_deaths_cols
  aggregated_rng_draws    <- age_rng_draws
  
  message(" > Estimation for island-spefic infections")
  
  for (k in 1:cov_data$K) {
    message(paste0(" > Estimation in island ", k))
    
    for (i in 1:mcmc_length) {
      for (j in 1:ts_length) {
        
        if(dispersion_type == "quadratic"){
          death_draws[i,j,k] <- rnbinom(1,
                                        mu   = lp[i,j,k],
                                        size = phi[i])
          
        } else {
          death_draws[i,j,k] <- rnbinom(1,
                                        mu   = lp[i,j,k],
                                        size = (lp[i,j,k] / phi[i]))
          
        }# End if
      }# End for
    }# End for
  }# End for
  
  for (k in 1:cov_data$K) {
    rng_draws_age_grp  <- data.frame(Date  = dates,
                                     Group = rep(key_mat$ISLAND[k], cov_data$W)
                                     )
    rng_draws_age_grp$median     <- apply(death_draws[,,k], 2, quantile, probs = c(0.50))
    rng_draws_age_grp$rng_low    <- apply(death_draws[,,k], 2, quantile, probs = c(0.025))
    rng_draws_age_grp$rng_low25  <- apply(death_draws[,,k], 2, quantile, probs = c(0.25))
    rng_draws_age_grp$rng_high75 <- apply(death_draws[,,k], 2, quantile, probs = c(0.75))
    rng_draws_age_grp$rng_high   <- apply(death_draws[,,k], 2, quantile, probs = c(0.975))
    
    age_rng_draws <- rbind(age_rng_draws, rng_draws_age_grp)
  }# End for
  
  ###################################
  age_rng_draws <- age_rng_draws %>% 
                   distinct(Date, Group, .keep_all = TRUE)
  
  return(age_rng_draws)
  
}# End function
cases_random_draws <- compiler::cmpfun(cases_random_draws)

cases_rng <- cases_random_draws(dispersion_type = "linear",
                                cov_data        = standata,
                                model_out       = nuts_fit_1,
                                main_dataset    = main_dataset_CHIKV,
                                key_mat         = key_CHIKV_filtered)

colnames(O_t)[-1] <- key_CHIKV_filtered$ISLAND

age_specific_clinicalcases_simdata_long <- tidyr::pivot_longer(O_t, 
                                                              cols      = -c("DATE2"), 
                                                              values_to = "Cases", 
                                                              names_to  = "Group") 
colnames(age_specific_clinicalcases_simdata_long)[1] <- "Date"

age_specific_clinicalcases_data <- age_specific_clinicalcases_simdata_long %>% 
                                   left_join(cases_rng, by = c("Date", "Group"))

posterior_clinicalcases <-
  ggplot(age_specific_clinicalcases_data,
         aes(Date  = Date,
             Cases = Cases)) +
  facet_wrap(. ~ Group, 
             scales = "free_y", 
             ncol   = 1,
             strip.position = "top") +
  geom_point(aes(x   = Date, 
                 y   = Cases)
  ) +
  geom_line(aes(x     = Date,
                y     = median,
                color = 'Median'),
            size = 1.0) +
  geom_ribbon(aes(x    = Date,
                  ymin = rng_low25,
                  ymax = rng_high75,
                  fill = "50% CrI"),
              alpha = 0.5) +
  geom_ribbon(aes(x    = Date,
                  ymin = rng_low,
                  ymax = rng_high,
                  fill = "95% CrI"),
              alpha = 0.5) +
  labs(x = "Date",
       y = "Clinical cases",
       fill = "") +
  scale_x_date(date_breaks       = "2 weeks", 
               date_minor_breaks = "2 weeks") +
  scale_fill_manual(values = c("50% CrI" = "gray30",
                               "95% CrI" = "gray70"
  )) +
  scale_colour_manual(name   = '',
                      values = c('Median' = "black")) +
  theme_bw() +
  theme(strip.placement  = "outside",
        axis.text.x      = element_text(size = 12, angle = 45, hjust = 1),
        axis.text.y      = element_text(size = 12),
        axis.title.x     = element_text(size = 16, face = "bold"),
        axis.title.y     = element_text(size = 16, face = "bold"),
        legend.text      = element_text(size = 14),
        strip.text.x     = element_text(size = 16),
        legend.position  = "bottom",
        legend.title     = element_blank() )
  
#########################################################
# 
# Island-specific effective reproduction number
# 
#########################################################

rt_mat_cols      <- c("Date", "Group", "median", "low", "low25", "high75", "high")
rt_mat           <- data.frame(matrix(ncol = length(rt_mat_cols), nrow = 0))
colnames(rt_mat) <- rt_mat_cols

for (i in 1:standata$K){
  
  fit_rt          <- posts_mcmc$rt_eff_island2[,,i] 
  
  dt_effcr_group  <- data.frame(Date     = unique(main_dataset_CHIKV$DATE2),
                                Group    = rep(key_CHIKV_filtered$ISLAND[i], standata$W)
  )
  
  # Add quantiles from the model outputs:
  dt_effcr_group$median <- apply(fit_rt, 2, median)
  dt_effcr_group$low    <- apply(fit_rt, 2, quantile, probs = c(0.025))
  dt_effcr_group$low25  <- apply(fit_rt, 2, quantile, probs = c(0.25))
  dt_effcr_group$high75 <- apply(fit_rt, 2, quantile, probs = c(0.75))
  dt_effcr_group$high   <- apply(fit_rt, 2, quantile, probs = c(0.975))
  
  rt_mat  <- rbind(rt_mat, dt_effcr_group)
  
}# End for

rt_mat$Group_f <- factor(rt_mat$Group, levels = key_CHIKV_filtered$ISLAND)

#---- Plot of posterior estimates:
posterior_rt <- 
  ggplot(rt_mat) +
  facet_wrap(. ~ Group, 
             scales = "free_y", 
             ncol   = 1,
             strip.position = "top") +
  geom_line(aes(x     = Date,
                y     = median,
                color = "Median"),
            size  = 1.0) +
  geom_ribbon(aes(x    = Date,
                  ymin = low25,
                  ymax = high75,
                  fill = "50% CrI"),
              alpha = 0.5) +
  geom_ribbon(aes(x    = Date,
                  ymin = low,
                  ymax = high,
                  fill = "95% CrI"),
              alpha = 0.5) +
  labs(x = "Date",
       y = "Effective reproduction number") +
  geom_hline(yintercept = 1, linetype = 2 ) + 
  scale_x_date(date_breaks       = "2 weeks", 
               date_minor_breaks = "2 weeks") + 
  scale_fill_manual(values = c("50% CrI" = "gray30",
                               "95% CrI" = "gray70"
  )) +
  scale_colour_manual(name   = '',
                      values = c('Median' = "black")) +
  theme_bw() +
  theme(strip.placement  = "outside",
        axis.text.x      = element_text(size = 12, angle = 45, hjust = 1),
        axis.text.y      = element_text(size = 12),
        axis.title.x     = element_text(size = 16, face = "bold"),
        axis.title.y     = element_text(size = 16, face = "bold"),
        legend.text      = element_text(size = 14),
        strip.text.x     = element_text(size = 16),
        legend.position  = "bottom",
        legend.title     = element_blank() )

#########################################################
#
# NegBin overdispersion parameter: 
#
#########################################################

# Create the appropriate dataset in a format acceptable by ggplot:
dt_phiD_post  <- data.frame(Posterior = posts_mcmc[["phi"]])
dt_phiD_post2 <- reshape2::melt(dt_phiD_post)

posterior_phi <- 
ggplot(dt_phiD_post2, 
       aes(x    = value)) +
  geom_density(alpha = 0.8) +
  scale_x_continuous(
    limits = c(min(dt_phiD_post2$value)*0.95, 
               max(dt_phiD_post2$value)*1.05)
  ) +
  labs(x = expression(phi),
       y = "Density") +
  theme_bw() +
  theme(panel.spacing    = unit(0.4,"cm"),
        axis.text.x      = element_text(angle = 0, hjust = 1),
        axis.title.x     = element_text(size = 14, face = "bold"),
        axis.title.y     = element_text(size = 14, face = "bold"),
        legend.position  = "bottom",
        legend.title     = element_blank(),
        legend.box       = "vertical", 
        legend.margin    = margin())


###########################################################################
#
# Traceplots of hyper-parameters
#
###########################################################################
posterior_1 <- as.array(nuts_fit_1)
colnames(posterior_1[1,,])

color_scheme_set("viridis")

if (model == "xGP"){
  
  sigma_mu_raw <- as.data.frame(posterior_1[, , 71]) # sigma_mu
  sigma_x_raw  <- as.data.frame(posterior_1[, , 72]) # sigma_x
  phi_raw      <- as.data.frame(posterior_1[, , 118]) # phi
  
  colnames(sigma_x_raw) <- colnames(sigma_mu_raw) <- colnames(phi_raw) <- as.character(1:ncol(sigma_x_raw))
  sigma_x_raw$Iteration  <- 1:nrow(sigma_x_raw)
  sigma_mu_raw$Iteration <- 1:nrow(sigma_mu_raw)
  phi_raw$Iteration      <- 1:nrow(phi_raw)
  
  sigma_x_raw$Variable  <- rep("sigma_x", nrow(sigma_x_raw))
  sigma_mu_raw$Variable <- rep("sigma_mu", nrow(sigma_mu_raw))
  phi_raw$Variable      <- rep("phi", nrow(phi_raw))
  
  sigmas_merged <- rbind(sigma_x_raw, sigma_mu_raw, phi_raw)
  
  sigmas_merged_long <- 
    pivot_longer(sigmas_merged, 
                 -c(Variable, Iteration), 
                 values_to = "value", 
                 names_to  = "Chain")
  
  sigmas_merged_long <-
    sigmas_merged_long %>%
    mutate(Parameter = factor(Variable,
                              levels = c("sigma_x", "sigma_mu", "phi"),
                              labels = c(bquote(sigma[x]),
                                         bquote(sigma[mu]),
                                         bquote(phi[.("C")]) ) ) )
  sigmas_merged_long$Variable<- NULL
  sigmas_merged_long$Chain <- as.numeric(sigmas_merged_long$Chain)
  
} else if (model == "mxGP"){
  
  sigma_mu_raw              <- as.data.frame(posterior_1[, , colnames(posterior_1[1,,])[colnames(posterior_1[1,,]) == "sigma_common"]]) 
  sigma_length_scale_mu_raw <- as.data.frame(posterior_1[, , colnames(posterior_1[1,,])[colnames(posterior_1[1,,]) == "length_scale_common"]])
  
  sigma_x_raw_1  <- as.data.frame(posterior_1[, , colnames(posterior_1[1,,])[colnames(posterior_1[1,,]) == "sigma[1]"]])
  sigma_x_raw_2  <- as.data.frame(posterior_1[, , colnames(posterior_1[1,,])[colnames(posterior_1[1,,]) == "sigma[2]"]])
  sigma_x_raw_3  <- as.data.frame(posterior_1[, , colnames(posterior_1[1,,])[colnames(posterior_1[1,,]) == "sigma[3]"]])
  
  sigma_length_scale_1_raw <- as.data.frame(posterior_1[, , colnames(posterior_1[1,,])[colnames(posterior_1[1,,]) == "length_scale[1]"]])
  sigma_length_scale_2_raw <- as.data.frame(posterior_1[, , colnames(posterior_1[1,,])[colnames(posterior_1[1,,]) == "length_scale[2]"]])
  sigma_length_scale_3_raw <- as.data.frame(posterior_1[, , colnames(posterior_1[1,,])[colnames(posterior_1[1,,]) == "length_scale[3]"]])

  phi_raw      <- as.data.frame(posterior_1[, , colnames(posterior_1[1,,])[colnames(posterior_1[1,,]) == "phi"]]) 
  
  colnames(sigma_mu_raw) <- 
    colnames(sigma_length_scale_mu_raw) <-
    colnames(sigma_x_raw_1) <- 
    colnames(sigma_x_raw_2) <- 
    colnames(sigma_x_raw_3) <- 
    colnames(sigma_length_scale_1_raw) <- 
    colnames(sigma_length_scale_2_raw) <- 
    colnames(sigma_length_scale_3_raw) <- 
    colnames(phi_raw) <- as.character(1:ncol(phi_raw))
  
  sigma_mu_raw$Iteration  <- 1:nrow(phi_raw)
  sigma_length_scale_mu_raw$Iteration  <- 1:nrow(phi_raw)
  sigma_x_raw_1$Iteration  <- 1:nrow(phi_raw)
  sigma_x_raw_2$Iteration  <- 1:nrow(phi_raw)
  sigma_x_raw_3$Iteration  <- 1:nrow(phi_raw)
  sigma_length_scale_1_raw$Iteration  <- 1:nrow(phi_raw)
  sigma_length_scale_2_raw$Iteration  <- 1:nrow(phi_raw)
  sigma_length_scale_3_raw$Iteration  <- 1:nrow(phi_raw)
  phi_raw$Iteration  <-1:nrow(phi_raw)
  
  
  sigma_x_raw$Variable  <- rep("sigma_x", nrow(sigma_x_raw))
  sigma_mu_raw$Variable <- rep("sigma_mu", nrow(sigma_mu_raw))

  sigma_mu_raw$Variable              <- rep("sigma_mu", nrow(phi_raw))
  sigma_length_scale_mu_raw$Variable <- rep("length_scale_mu", nrow(phi_raw))
  sigma_x_raw_1$Variable             <- rep("sigma_x_1", nrow(phi_raw))
  sigma_x_raw_2$Variable             <- rep("sigma_x_2", nrow(phi_raw))
  sigma_x_raw_3$Variable             <- rep("sigma_x_3", nrow(phi_raw))
  sigma_length_scale_1_raw$Variable  <- rep("length_scale_1", nrow(phi_raw))
  sigma_length_scale_2_raw$Variable  <- rep("length_scale_2", nrow(phi_raw))
  sigma_length_scale_3_raw$Variable  <- rep("length_scale_3", nrow(phi_raw))
  phi_raw$Variable                   <- rep("phi", nrow(phi_raw))
  
  sigmas_merged <- rbind(sigma_mu_raw, 
                         sigma_length_scale_mu_raw, 
                         sigma_x_raw_1, 
                         sigma_x_raw_2, 
                         sigma_x_raw_3, 
                         sigma_length_scale_1_raw,
                         sigma_length_scale_2_raw,
                         sigma_length_scale_3_raw,
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
                                         "length_scale_mu", 
                                         "sigma_x_1", 
                                         "sigma_x_2", 
                                         "sigma_x_3",
                                         "length_scale_1",
                                         "length_scale_2",
                                         "length_scale_3",
                                         "phi"),
                              labels = c(bquote(sigma[mu]),
                                         expression("\u2113"[mu]),
                                         expression(sigma[paste(x, ",", 1)]),
                                         expression(sigma[paste(x, ",", 2)]),
                                         expression(sigma[paste(x, ",", 3)]),
                                         expression(l[paste(x, ",", 1)]),
                                         expression(l[paste(x, ",", 2)]),
                                         expression(l[paste(x, ",", 3)]),
                                         bquote(phi[.("C")]) ) ) )
  sigmas_merged_long$Variable<- NULL
  sigmas_merged_long$Chain <- as.numeric(sigmas_merged_long$Chain)

} else if (model == "iBM"){
  
  sigma_x_raw_1  <- as.data.frame(posterior_1[, , colnames(posterior_1[1,,])[colnames(posterior_1[1,,]) == "sigma[1]"]])
  sigma_x_raw_2  <- as.data.frame(posterior_1[, , colnames(posterior_1[1,,])[colnames(posterior_1[1,,]) == "sigma[2]"]])
  sigma_x_raw_3  <- as.data.frame(posterior_1[, , colnames(posterior_1[1,,])[colnames(posterior_1[1,,]) == "sigma[3]"]])
  
  phi_raw      <- as.data.frame(posterior_1[, , colnames(posterior_1[1,,])[colnames(posterior_1[1,,]) == "phi"]]) 
  
  colnames(sigma_x_raw_1) <- 
  colnames(sigma_x_raw_2) <- 
  colnames(sigma_x_raw_3) <- 
  colnames(phi_raw) <- as.character(1:ncol(phi_raw))
  
  sigma_x_raw_1$Iteration  <- 1:nrow(phi_raw)
  sigma_x_raw_2$Iteration  <- 1:nrow(phi_raw)
  sigma_x_raw_3$Iteration  <- 1:nrow(phi_raw)
  phi_raw$Iteration        <-1:nrow(phi_raw)
  
  sigma_x_raw_1$Variable             <- rep("sigma_x_1", nrow(phi_raw))
  sigma_x_raw_2$Variable             <- rep("sigma_x_2", nrow(phi_raw))
  sigma_x_raw_3$Variable             <- rep("sigma_x_3", nrow(phi_raw))
  phi_raw$Variable                   <- rep("phi", nrow(phi_raw))
  
  sigmas_merged <- rbind(
                         sigma_x_raw_1, 
                         sigma_x_raw_2, 
                         sigma_x_raw_3, 
                         phi_raw)
  
  sigmas_merged_long <- 
    pivot_longer(sigmas_merged, 
                 -c(Variable, Iteration), 
                 values_to = "value", 
                 names_to  = "Chain")
  
  sigmas_merged_long <-
    sigmas_merged_long %>%
    mutate(Parameter = factor(Variable,
                              levels = c("sigma_x_1", 
                                         "sigma_x_2", 
                                         "sigma_x_3",
                                         "phi"),
                              labels = c(expression(sigma[paste(x, ",", 1)]),
                                         expression(sigma[paste(x, ",", 2)]),
                                         expression(sigma[paste(x, ",", 3)]),
                                         bquote(phi[.("C")]) ) ) )
  sigmas_merged_long$Variable<- NULL
  sigmas_merged_long$Chain <- as.numeric(sigmas_merged_long$Chain)
  
}

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

#########################################################
#
# Print the combined plot
#
#########################################################
library(patchwork)
combined <-
  (posterior_rt +
   posterior_beta + # / 
   posterior_clinicalcases) + 
  plot_layout(guides = 'collect') &
  theme(legend.position = 'bottom', legend.direction = 'horizontal')

grDevices::cairo_ps(filename = "....//Chikv//6_Model_plots//CHIKV_FP_2015_xBM_posteriorpredictive_rt.eps",
         width     = 10, 
         height    = 10, 
         pointsize = 12,
         fallback_resolution = 300)
print(combined + plot_annotation( tag_levels = list(c("A", "B", "C") ))  )
dev.off()
