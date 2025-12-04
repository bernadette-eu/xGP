mod1_diagnostics <- rstan::get_sampler_params(nuts_fit_1)
posts_1          <- rstan::extract(nuts_fit_1)
lp_cp            <- bayesplot::log_posterior(nuts_fit_1)
np_cp            <- bayesplot::nuts_params(nuts_fit_1)
posterior_1      <- as.array(nuts_fit_1)

##########################################################################
#
# MCMC summaries
#
###########################################################################
nuts_fit_1_summary <- summary(nuts_fit_1)$summary

options(scipen = 999)
round(nuts_fit_1_summary, 3)

############################################################
# Look at treedepth__
############################################################
unique(np_cp$Parameter)

check_treedepth <- subset(np_cp, Parameter == "treedepth__")

ggplot(check_treedepth,
       aes(Iteration = Iteration,
           Value     = Value)) +
  facet_wrap(. ~ Chain, 
             scales = "free_y", 
             ncol   = 1,
             strip.position = "right") +
  geom_line(aes(x     = Iteration,
                y     = Value),
            size = 1.3) +
  labs(x = "Iteration",
       y = "treedepth__") +
  theme_bw() +
  theme(panel.spacing    = unit(0.2,"cm"),
        axis.text.x      = element_text(angle = 0, hjust = 1),
        axis.title.x     = element_text(size = 14, face = "bold"),
        axis.title.y     = element_text(size = 14, face = "bold"),
        legend.position  = "bottom",
        legend.title     = element_blank(),
        legend.box       = "vertical", 
        legend.margin    = margin() )

############################################################
# Look at n_leapfrog__
############################################################
check_leapfrog <- subset(np_cp, Parameter == "n_leapfrog__")

ggplot(check_leapfrog,
       aes(Iteration = Iteration,
           Value     = Value)) +
  facet_wrap(. ~ Chain, 
             scales = "free_y", 
             ncol   = 1,
             strip.position = "right") +
  geom_line(aes(x     = Iteration,
                y     = Value),
            size = 1.3) +
  labs(x = "Iteration",
       y = "n_leapfrog__") +
  theme_bw() +
  theme(panel.spacing    = unit(0.2,"cm"),
        axis.text.x      = element_text(angle = 0, hjust = 1),
        axis.title.x     = element_text(size = 14, face = "bold"),
        axis.title.y     = element_text(size = 14, face = "bold"),
        legend.position  = "bottom",
        legend.title     = element_blank(),
        legend.box       = "vertical", 
        legend.margin    = margin() )


############################################################
# Look at accept_stat__
############################################################
check_acceptstat<- subset(np_cp, Parameter == "accept_stat__")

ggplot(check_acceptstat,
       aes(Iteration = Iteration,
           Value     = Value)) +
  facet_wrap(. ~ Chain, 
             scales = "free_y", 
             ncol   = 1,
             strip.position = "right") +
  geom_line(aes(x     = Iteration,
                y     = Value),
            size = 1.3) +
  labs(x = "Iteration",
       y = "accept_stat__") +
  theme_bw() +
  theme(panel.spacing    = unit(0.2,"cm"),
        axis.text.x      = element_text(angle = 0, hjust = 1),
        axis.title.x     = element_text(size = 14, face = "bold"),
        axis.title.y     = element_text(size = 14, face = "bold"),
        legend.position  = "bottom",
        legend.title     = element_blank(),
        legend.box       = "vertical", 
        legend.margin    = margin() )

############################################################
# Pairs plots
############################################################
nuts_fit_1@model_pars

# Transmissibility vs volatilities:
pairs(nuts_fit_1, pars = c("beta_N[40,1]",
                           "beta_N[40,2]",
                           #"beta_N[40,3]",
                           "sigmaBM[1]",
                           "sigmaBM[2]"#,
                           #"sigmaBM[3]"
                           ))

pairs(nuts_fit_1, pars = c("beta_N[50,1]",
                           "beta_N[50,2]",
                           "beta_N[50,3]",
                           "sigmaBM[1]",
                           "sigmaBM[2]",
                           "sigmaBM[3]"))

# Contact matrix:
pairs(nuts_fit_1, pars = c("cm_sample[1,1]",
                           "cm_sample[1,2]",
                           "cm_sample[1,3]",
                           "cm_sample[2,1]",
                           "cm_sample[2,2]",
                           "cm_sample[2,3]",
                           "cm_sample[3,1]",
                           "cm_sample[3,2]",
                           "cm_sample[3,3]"))

pairs(nuts_fit_1, pars = c("beta_N[150,1]",
                           "beta_N[150,2]",
                           "beta_N[150,3]",
                           "cm_sample[1,1]",
                           "cm_sample[1,2]",
                           "cm_sample[1,3]",
                           "cm_sample[2,1]"))

# Transmissibility vs volatilities vs overdispersion:
pairs(nuts_fit_1, pars = c("beta_N[40,1]",
                           "beta_N[40,2]",
                           #"beta_N[40,3]",
                           "sigmaBM[1]",
                           "sigmaBM[2]",
                           #"sigmaBM[3]",
                           "phiD"))

# Volatilities vs overdispersion:
pairs(nuts_fit_1, pars = c("sigmaBM[1]",
                           "sigmaBM[2]",
                           #"sigmaBM[3]",
                           "phiD"))
