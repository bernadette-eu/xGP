lib <- c("loo",
         "bridgesampling"
)
lapply(lib, require, character.only = TRUE)

#--- Paths:
main_path <- "...//Chikv//4_Stan_Outputs//"

#---- Paths to model outputs:
model_ids <- c("Main_Model_CHIKV_FWI_v1",       # Baseline model
               "CHIKV_FWI_2014_IGBM_v5_chol",   
               "CHIKV_FWI_2014_XGBM_v3_chol",   
               "CHIKV_FWI_2014_MXGBM_v3_chol", 
               "CHIKV_FWI_2014_iGP",
               "CHIKV_FWI_2014_xGP",
               "CHIKV_FWI_2014_mxGP"
               )

store_model_out_model1 <- paste0(main_path, model_ids[1], ".RData")
store_model_out_model2 <- paste0(main_path, model_ids[2], ".RData")
store_model_out_model3 <- paste0(main_path, model_ids[3], ".RData")
store_model_out_model4 <- paste0(main_path, model_ids[4], ".RData")
nChains                <- 6

#---- Loo, model 1:
load(file = store_model_out_model1)

log_lik_1 <- loo::extract_log_lik(nuts_fit_1, merge_chains = FALSE)
r_eff_1   <- loo::relative_eff(exp(log_lik_1),    cores = nChains)
loo_1     <- loo::loo(log_lik_1, r_eff = r_eff_1, cores = nChains)
print(loo_1)

# Export to latex format:
model1_summary <- c(paste0(round(loo_1$estimates[1,1],1),  " (",  round(loo_1$estimates[1,2],1), ")"),
                    round(loo_1$estimates[3,1] - 2*loo_1$estimates[2,1], 1),
                    paste0(round(loo_1$estimates[2,1],1),  " (",  round(loo_1$estimates[2,2],1), ")"),
                    paste0(round(loo_1$estimates[3,1],1),  " (",  round(loo_1$estimates[3,2],1), ")"),
                    paste0(round(bridgeBase$logml, 1), " (",  bridgesampling::error_measures(bridgeBase)$percentage, ")")
)
names(model1_summary)<- c("ELPD", 
                          "Adequacy", 
                          "Complexity", 
                          "LooIC", 
                          "MarginalLikelihood")

xtable(t(data.frame(model1_summary)))

model1_summary2 <- c(round(loo_1$estimates[3,1] - 2*loo_1$estimates[2,1], 1),
                     paste0(round(loo_1$estimates[2,1],1),  " (",  round(loo_1$estimates[2,2],1), ")"),
                     paste0(round(loo_1$estimates[3,1],1),  " (",  round(loo_1$estimates[3,2],1), ")"),
                     paste0(round(bridgeBase$logml, 1), " (",  bridgesampling::error_measures(bridgeBase)$percentage, ")")
)
names(model1_summary2)<- c("Adequacy", "Complexity", "LooIC", "MarginalLikelihood")

xtable(t(data.frame(model1_summary2)))

rm(nuts_fit_1)

#---- Loo, model 2:
load(file = store_model_out_model2)

log_lik_2 <- loo::extract_log_lik(nuts_fit_1, merge_chains = FALSE)
r_eff_2   <- loo::relative_eff(exp(log_lik_2),    cores = nChains)
loo_2     <- loo::loo(log_lik_2, r_eff = r_eff_2, cores = nChains)
print(loo_2)

# Export to latex format:
model2_summary <- c(paste0(round(waic_2$estimates[3,1],1), " (",  round(waic_2$estimates[3,2],1), ")"),
                    paste0(round(loo_2$estimates[3,1],1),  " (",  round(loo_2$estimates[3,2],1), ")"),
                    paste0(round(loo_2$estimates[1,1],1),  " (",  round(loo_2$estimates[1,2],1), ")"),
                    round(loo_2$estimates[3,1] - 2*loo_2$estimates[2,1], 1),
                    paste0(round(loo_2$estimates[2,1],1),  " (",  round(loo_2$estimates[2,2],1), ")"),
                    paste0(round(bridgeIgbm$logml, 1), " (",  bridgesampling::error_measures(bridgeIgbm)$percentage, ")")
)
names(model2_summary)<- c("WAIC", "LooIC", "ELPD", "Adequacy", "Complexity", "MarginalLikelihood")

xtable(t(data.frame(model2_summary)))

model2_summary2 <- c(round(loo_2$estimates[3,1] - 2*loo_2$estimates[2,1], 1),
                     paste0(round(loo_2$estimates[2,1],1),  " (",  round(loo_2$estimates[2,2],1), ")"),
                     paste0(round(loo_2$estimates[3,1],1),  " (",  round(loo_2$estimates[3,2],1), ")"),
                     paste0(round(bridgeIgbm$logml, 1), " (",  bridgesampling::error_measures(bridgeIgbm)$percentage, ")")
)
names(model2_summary2)<- c("Adequacy", 
                           "Complexity", 
                           "LooIC", 
                           "MarginalLikelihood")

xtable(t(data.frame(model2_summary2)))

rm(nuts_fit_1)

#---- Loo, model 3:
load(file = store_model_out_model3)

log_lik_3 <- loo::extract_log_lik(nuts_fit_1, merge_chains = FALSE)
r_eff_3   <- loo::relative_eff(exp(log_lik_3),    cores = nChains)
loo_3     <- loo::loo(log_lik_3, r_eff = r_eff_3, cores = nChains)
print(loo_3)

# Export to latex format:
model3_summary <- c(paste0(round(waic_3$estimates[3,1],1), " (",  round(waic_3$estimates[3,2],1), ")"),
                    paste0(round(loo_3$estimates[3,1],1),  " (",  round(loo_3$estimates[3,2],1), ")"),
                    paste0(round(loo_3$estimates[1,1],1),  " (",  round(loo_3$estimates[1,2],1), ")"),
                    round(loo_3$estimates[3,1] - 2*loo_3$estimates[2,1], 1),
                    paste0(round(loo_3$estimates[2,1],1),  " (",  round(loo_3$estimates[2,2],1), ")"),
                    paste0(round(bridgeXgbm$logml, 1), " (",  bridgesampling::error_measures(bridgeXgbm)$percentage, ")")
)
names(model3_summary)<- c("WAIC", "LooIC", "ELPD", "Adequacy", "Complexity", "MarginalLikelihood")

xtable(t(data.frame(model3_summary)))

model3_summary2 <- c(round(loo_3$estimates[3,1] - 2*loo_3$estimates[2,1], 1),
                     paste0(round(loo_3$estimates[2,1],1),  " (",  round(loo_3$estimates[2,2],1), ")"),
                     paste0(round(loo_3$estimates[3,1],1),  " (",  round(loo_3$estimates[3,2],1), ")"),
                     paste0(round(bridgeXgbm$logml, 1), " (",  bridgesampling::error_measures(bridgeXgbm)$percentage, ")")
)
names(model3_summary2)<- c("Adequacy", 
                           "Complexity", 
                           "LooIC", 
                           "MarginalLikelihood")

xtable(t(data.frame(model3_summary2)))

rm(nuts_fit_1)

#---- Loo, model 4:
load(file = store_model_out_model4)

log_lik_4 <- loo::extract_log_lik(nuts_fit_1, merge_chains = FALSE)
r_eff_4   <- loo::relative_eff(exp(log_lik_4),    cores = nChains)
loo_4     <- loo::loo(log_lik_4, r_eff = r_eff_4, cores = nChains)

print(loo_4)

# Export to latex format:
model4_summary <- c(paste0(round(waic_4$estimates[3,1],1), " (",  round(waic_4$estimates[3,2],1), ")"),
                    paste0(round(loo_4$estimates[3,1],1),  " (",  round(loo_4$estimates[3,2],1), ")"),
                    paste0(round(loo_4$estimates[1,1],1),  " (",  round(loo_4$estimates[1,2],1), ")"),
                    round(loo_4$estimates[3,1] - 2*loo_4$estimates[2,1], 1),
                    paste0(round(loo_4$estimates[2,1],1),  " (",  round(loo_4$estimates[2,2],1), ")"),
                    paste0(round(bridgemxgbm$logml,1), " (",  bridgesampling::error_measures(bridgemxgbm)$percentage, ")")
)
names(model3_summary)<- c("WAIC", "LooIC", "ELPD", "Adequacy", "Complexity", "MarginalLikelihood")

xtable(t(data.frame(model4_summary)))

model4_summary2 <- c(round(loo_4$estimates[3,1] - 2*loo_4$estimates[2,1], 1),
                     paste0(round(loo_4$estimates[2,1],1),  " (",  round(loo_4$estimates[2,2],1), ")"),
                     paste0(round(loo_4$estimates[3,1],1),  " (",  round(loo_4$estimates[3,2],1), ")"),
                     paste0(round(bridgemxgbm$logml,1), " (",  bridgesampling::error_measures(bridgemxgbm)$percentage, ")")
)
names(model4_summary2)<- c("Adequacy", 
                           "Complexity", 
                           "LooIC", 
                           "MarginalLikelihood")

xtable(t(data.frame(model4_summary2)))

rm(nuts_fit_1)

#---- Loo, model 5:
load(paste0(main_path, model_ids[5], ".RData"))

log_lik_5 <- loo::extract_log_lik(nuts_fit_iGP, merge_chains = FALSE)
r_eff_5   <- loo::relative_eff(exp(log_lik_5),    cores = nChains)
loo_5     <- loo::loo(log_lik_5, r_eff = r_eff_5, cores = nChains)

# Export to latex format:
model5_summary <- c(paste0(round(waic_5$estimates[3,1],1), " (",  round(waic_5$estimates[3,2],1), ")"),
                    paste0(round(loo_5$estimates[3,1],1),  " (",  round(loo_5$estimates[3,2],1), ")"),
                    paste0(round(loo_5$estimates[1,1],1),  " (",  round(loo_5$estimates[1,2],1), ")"),
                    round(loo_5$estimates[3,1] - 2*loo_5$estimates[2,1], 1),
                    paste0(round(loo_5$estimates[2,1],1),  " (",  round(loo_5$estimates[2,2],1), ")"),
                    paste0(round(bridgeiGP$logml,1), " (",  bridgesampling::error_measures(bridgeiGP)$percentage, ")")
)
names(model5_summary)<- c("WAIC", "LooIC", "ELPD", "Adequacy", "Complexity", "MarginalLikelihood")

xtable(t(data.frame(model5_summary)))

model5_summary2 <- c(round(loo_5$estimates[3,1] - 2*loo_5$estimates[2,1], 1),
                     paste0(round(loo_5$estimates[2,1],1),  " (",  round(loo_5$estimates[2,2],1), ")"),
                     paste0(round(loo_5$estimates[3,1],1),  " (",  round(loo_5$estimates[3,2],1), ")"),
                     paste0(round(bridgeiGP$logml,1), " (",  bridgesampling::error_measures(bridgeiGP)$percentage, ")")
)
names(model5_summary2)<- c("Adequacy", 
                           "Complexity", 
                           "LooIC", 
                           "MarginalLikelihood")

xtable(t(data.frame(model5_summary2)))

rm(nuts_fit_iGP)

#---- Loo, model 6:
load(paste0(main_path, model_ids[6], ".RData"))

log_lik_6 <- loo::extract_log_lik(nuts_fit_xGP, merge_chains = FALSE)
r_eff_6   <- loo::relative_eff(exp(log_lik_6),    cores = nChains)
loo_6     <- loo::loo(log_lik_6, r_eff = r_eff_6, cores = nChains)
print(loo_6)

# Export to latex format:
model6_summary <- c(paste0(round(waic_6$estimates[3,1],1), " (",  round(waic_6$estimates[3,2],1), ")"),
                    paste0(round(loo_6$estimates[3,1],1),  " (",  round(loo_6$estimates[3,2],1), ")"),
                    paste0(round(loo_6$estimates[1,1],1),  " (",  round(loo_6$estimates[1,2],1), ")"),
                    round(loo_6$estimates[3,1] - 2*loo_6$estimates[2,1], 1),
                    paste0(round(loo_6$estimates[2,1],1),  " (",  round(loo_6$estimates[2,2],1), ")"),
                    paste0(round(bridgexGP$logml,1), " (",  bridgesampling::error_measures(bridgexGP)$percentage, ")")
)
names(model6_summary)<- c("WAIC", "LooIC", "ELPD", "Adequacy", "Complexity", "MarginalLikelihood")

xtable(t(data.frame(model6_summary)))

model6_summary2 <- c(round(loo_6$estimates[3,1] - 2*loo_6$estimates[2,1], 1),
                     paste0(round(loo_6$estimates[2,1],1),  " (",  round(loo_6$estimates[2,2],1), ")"),
                     paste0(round(loo_6$estimates[3,1],1),  " (",  round(loo_6$estimates[3,2],1), ")"),
                     paste0(round(bridgexGP$logml,1), " (",  bridgesampling::error_measures(bridgexGP)$percentage, ")")
)
names(model6_summary2)<- c("Adequacy", 
                           "Complexity", 
                           "LooIC", 
                           "MarginalLikelihood")

xtable(t(data.frame(model6_summary2)))

rm(nuts_fit_xGP)

#---- Loo, model 7:
load(paste0(main_path, model_ids[7], ".RData"))

log_lik_7 <- loo::extract_log_lik(nuts_fit_mxGP, merge_chains = FALSE)
r_eff_7   <- loo::relative_eff(exp(log_lik_7),    cores = nChains)
loo_7     <- loo::loo(log_lik_7, r_eff = r_eff_7, cores = nChains)

print(loo_7)
#print(waic_7)

# Export to latex format:
model7_summary <- c(paste0(round(waic_7$estimates[3,1],1), " (",  round(waic_7$estimates[3,2],1), ")"),
                    paste0(round(loo_7$estimates[3,1],1),  " (",  round(loo_7$estimates[3,2],1), ")"),
                    paste0(round(loo_7$estimates[1,1],1),  " (",  round(loo_7$estimates[1,2],1), ")"),
                    round(loo_7$estimates[3,1] - 2*loo_7$estimates[2,1], 1),
                    paste0(round(loo_7$estimates[2,1],1),  " (",  round(loo_7$estimates[2,2],1), ")"),
                    paste0(round(bridgemxGP$logml,1), " (",  bridgesampling::error_measures(bridgemxGP)$percentage, ")")
)
names(model7_summary)<- c("WAIC", "LooIC", "ELPD", "Adequacy", "Complexity", "MarginalLikelihood")

xtable(t(data.frame(model7_summary)))

model7_summary2 <- c(round(loo_7$estimates[3,1] - 2*loo_7$estimates[2,1], 1),
                     paste0(round(loo_7$estimates[2,1],1),  " (",  round(loo_7$estimates[2,2],1), ")"),
                     paste0(round(loo_7$estimates[3,1],1),  " (",  round(loo_7$estimates[3,2],1), ")"),
                     paste0(round(bridgemxGP$logml,1), " (",  bridgesampling::error_measures(bridgemxGP)$percentage, ")")
)
names(model7_summary2)<- c("Adequacy", 
                           "Complexity", 
                           "LooIC", 
                           "MarginalLikelihood")

xtable(t(data.frame(model7_summary2)))

rm(nuts_fit_mxGP)

#---- Model comparison:
comp2 <- loo::loo_compare(loo_1, loo_2, loo_3, loo_4, loo_5, loo_6, loo_7)
print(comp2)

#---- Report LooIC, Adequacy, Complexity, Delta LooIC and compare models by their information criteria:
round(cbind(comp2[,1] * -2, comp2[,2] * 2), 1)

# Approximate 95% CI:
round(comp2[2,1] + c(-1,1) * 1.96 * comp2[2,2] , 2)

looic_comparison <- round(cbind(comp2[,1] * -2, comp2[,2] * 2), 1)
looic_comparison <- looic_comparison[order(rownames(looic_comparison)),]

#---- Posterior model probabilities:
post1 <- post_prob(bridgeBase, 
                   bridgeIgbm, 
                   bridgeXgbm, 
                   bridgemxgbm,
                   bridgeiGP,
                   bridgexGP,
                   bridgemxGP)
round( print(post1), 3)

#---- Pseudo-BMA+ weights for stacking:
pseudobma_weights(cbind(loo_1$pointwise[,1], 
                        loo_2$pointwise[,1], 
                        loo_3$pointwise[,1],
                        loo_4$pointwise[,1]))