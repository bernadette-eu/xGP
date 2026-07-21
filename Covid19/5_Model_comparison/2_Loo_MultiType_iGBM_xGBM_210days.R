library(loo)

#--- Paths:
# pc_path        <- "C://Users//bouranis//"
# project_path   <- paste0(pc_path,"OneDrive - aueb.gr/BERNADETTE/")
# lib_path       <- paste0(project_path,
#                          "8_Package/bernadette/R/")
# experiments_path <- "10_Stan_Project_code/13_Exchange_BM"

#---- Paths to model outputs:
main_path       <- "C://Users//lbour//OneDrive//Desktop1//xgp_project//GP multitype epidemic England//4_Stan_Outputs//"

# Set 1 (Analysis period = 210 days)
# 1. SEEIIR_renewal_multitype_CGBM_v1.R    [cGBM]
# 2. SEEIIR_renewal_multitype_XGBM_v4E2.R  [xGBM]
# 3. SEEIIR_renewal_multitype_XGBM_v5C.R   [mxGBM]
# 4. SEEIIR_renewal_multitype_IGBM_v3C.R   [iGBM]

model_ids <- c(#"SEEIIR_renewal_multitype_CGBM_v1",
              "SEEIIR_renewal_multitype_IGBM_v3E",
              "SEEIIR_renewal_multitype_XGBM",
              "SEEIIR_renewal_multitype_MXGBM",
              "iGP_3cp_loglik", #"fit_ENG_LABROS_DATA_iGP_3cp_2000iter",
              "mGP_3cp_loglik", #"fit_ENG_LABROS_DATA_xGP_3cp_2000iter",
              "mxGP_3cp_loglik" #"fit_ENG_LABROS_DATA_mxGP_3cp_2000iter"
              # "SEEIIR_logLik_iGBM",  # "SEEIIR_renewal_multitype_IGBM_v3C",
              # "SEEIIR_logLik_xGBM",  # "SEEIIR_renewal_multitype_XGBM_v4E2",
              # "SEEIIR_logLik_mxGBM", # "SEEIIR_renewal_multitype_XGBM_v5C",
              # "fit_ENG_LABROS_DATA_iGP_2000iter",
              # "fit_ENG_LABROS_DATA_xGP_2000iter",
              # "fit_ENG_LABROS_DATA_mxGP_2000iter"
              )

store_model_out_model1 <- paste0(main_path, 
                                 model_ids[1],
                                 ".RData")
store_model_out_model2 <- paste0(main_path,
                                 model_ids[2],
                                 ".RData")
store_model_out_model3 <- paste0(main_path, 
                                 model_ids[3],
                                 ".RData")
store_model_out_model4 <- paste0(main_path, 
                                 model_ids[4],
                                 ".RData")
store_model_out_model5 <- paste0(main_path, 
                                 model_ids[5],
                                 ".RData")
store_model_out_model6 <- paste0(main_path,
                                 model_ids[6],
                                 ".RData")
# store_model_out_model7 <- paste0(main_path, 
#                                  model_ids[7],
#                                  ".RData")

nChains   <- 6

#---- Loo, model 1:
load(file = store_model_out_model1)

log_lik_1 <- loo::extract_log_lik(nuts_fit_1, merge_chains = FALSE, parameter_name = "log_like_age") #iGBM_gqs
r_eff_1   <- loo::relative_eff(exp(log_lik_1),    cores = nChains)
loo_1     <- loo::loo(log_lik_1, r_eff = r_eff_1, cores = nChains)
waic_1    <- loo::waic(log_lik_1)

print(loo_1)
print(waic_1)

# Export to latex format:
model1_summary <- c(paste0(round(waic_1$estimates[3,1],1), " (",  round(waic_1$estimates[3,2],1), ")"),
                    paste0(round(loo_1$estimates[3,1],1),  " (",  round(loo_1$estimates[3,2],1), ")"),
                    paste0(round(loo_1$estimates[1,1],1),  " (",  round(loo_1$estimates[1,2],1), ")"),
                    round(loo_1$estimates[3,1] - 2*loo_1$estimates[2,1], 1),
                    paste0(round(loo_1$estimates[2,1],1),  " (",  round(loo_1$estimates[2,2],1), ")")
)
names(model1_summary)<- c("WAIC", "LooIC", "ELPD", "Adequacy", "Complexity")

xtable(t(data.frame(model1_summary)))

model1_summary2 <- c(round(loo_1$estimates[3,1] - 2*loo_1$estimates[2,1], 1),
                     paste0(round(loo_1$estimates[2,1],1),  " (",  round(loo_1$estimates[2,2],1), ")"),
                     paste0(round(loo_1$estimates[3,1],1),  " (",  round(loo_1$estimates[3,2],1), ")")
)
names(model1_summary2)<- c("Adequacy", "Complexity", "LooIC")

xtable(t(data.frame(model1_summary2)))

rm(nuts_fit_1)

#---- Loo, model 2:
load(file = store_model_out_model2)

log_lik_2 <- loo::extract_log_lik(nuts_fit_1, merge_chains = FALSE, parameter_name = "log_like_age")
r_eff_2   <- loo::relative_eff(exp(log_lik_2),    cores = nChains)
loo_2     <- loo::loo(log_lik_2, r_eff = r_eff_2, cores = nChains)
waic_2    <- loo::waic(log_lik_2)

print(loo_2)
print(waic_2)

# Export to latex format:
model2_summary <- c(paste0(round(waic_2$estimates[3,1],1), " (",  round(waic_2$estimates[3,2],1), ")"),
                    paste0(round(loo_2$estimates[3,1],1),  " (",  round(loo_2$estimates[3,2],1), ")"),
                    paste0(round(loo_2$estimates[1,1],1),  " (",  round(loo_2$estimates[1,2],1), ")"),
                    round(loo_2$estimates[3,1] - 2*loo_2$estimates[2,1], 1),
                    paste0(round(loo_2$estimates[2,1],1),  " (",  round(loo_2$estimates[2,2],1), ")")
)
names(model2_summary)<- c("WAIC", "LooIC", "ELPD", "Adequacy", "Complexity")

xtable(t(data.frame(model2_summary)))

model2_summary2 <- c(round(loo_2$estimates[3,1] - 2*loo_2$estimates[2,1], 1),
                     paste0(round(loo_2$estimates[2,1],1),  " (",  round(loo_2$estimates[2,2],1), ")"),
                     paste0(round(loo_2$estimates[3,1],1),  " (",  round(loo_2$estimates[3,2],1), ")")
)
names(model2_summary2)<- c("Adequacy", "Complexity", "LooIC")

xtable(t(data.frame(model2_summary2)))

rm(nuts_fit_1)

#---- Loo, model 3:
load(file = store_model_out_model3)

log_lik_3 <- loo::extract_log_lik(nuts_fit_1, merge_chains = FALSE, parameter_name = "log_like_age")
r_eff_3   <- loo::relative_eff(exp(log_lik_3),    cores = nChains)
loo_3     <- loo::loo(log_lik_3, r_eff = r_eff_3, cores = nChains)
waic_3    <- loo::waic(log_lik_3)

print(loo_3)
print(waic_3)

# Export to latex format:
model3_summary <- c(paste0(round(waic_3$estimates[3,1],1), " (",  round(waic_3$estimates[3,2],1), ")"),
                    paste0(round(loo_3$estimates[3,1],1),  " (",  round(loo_3$estimates[3,2],1), ")"),
                    paste0(round(loo_3$estimates[1,1],1),  " (",  round(loo_3$estimates[1,2],1), ")"),
                    round(loo_3$estimates[3,1] - 2*loo_3$estimates[2,1], 1),
                    paste0(round(loo_3$estimates[2,1],1),  " (",  round(loo_3$estimates[2,2],1), ")")
)
names(model3_summary)<- c("WAIC", "LooIC", "ELPD", "Adequacy", "Complexity")

xtable(t(data.frame(model3_summary)))

model3_summary2 <- c(round(loo_3$estimates[3,1] - 2*loo_3$estimates[2,1], 1),
                     paste0(round(loo_3$estimates[2,1],1),  " (",  round(loo_3$estimates[2,2],1), ")"),
                     paste0(round(loo_3$estimates[3,1],1),  " (",  round(loo_3$estimates[3,2],1), ")")
)
names(model3_summary2)<- c("Adequacy", "Complexity", "LooIC")

xtable(t(data.frame(model3_summary2)))

rm(nuts_fit_1)

#---- Loo, model 4:
load(file = store_model_out_model4)

log_lik_4 <- loo::extract_log_lik(iGP_gqs, merge_chains = FALSE, parameter_name = "log_like_age")
r_eff_4   <- loo::relative_eff(exp(log_lik_4),    cores = nChains)
loo_4     <- loo::loo(log_lik_4, r_eff = r_eff_4, cores = nChains)
waic_4    <- loo::waic(log_lik_4)

print(loo_4)
print(waic_4)

# Export to latex format:
model4_summary <- c(paste0(round(waic_4$estimates[3,1],1), " (",  round(waic_4$estimates[3,2],1), ")"),
                    paste0(round(loo_4$estimates[3,1],1),  " (",  round(loo_4$estimates[3,2],1), ")"),
                    paste0(round(loo_4$estimates[1,1],1),  " (",  round(loo_4$estimates[1,2],1), ")"),
                    round(loo_4$estimates[3,1] - 2*loo_4$estimates[2,1], 1),
                    paste0(round(loo_4$estimates[2,1],1),  " (",  round(loo_4$estimates[2,2],1), ")")
)
names(model4_summary)<- c("WAIC", "LooIC", "ELPD", "Adequacy", "Complexity")

xtable(t(data.frame(model4_summary)))

model4_summary2 <- c(round(loo_4$estimates[3,1] - 2*loo_4$estimates[2,1], 1),
                     paste0(round(loo_4$estimates[2,1],1),  " (",  round(loo_4$estimates[2,2],1), ")"),
                     paste0(round(loo_4$estimates[3,1],1),  " (",  round(loo_4$estimates[3,2],1), ")")
)
names(model4_summary2)<- c("Adequacy", "Complexity", "LooIC")

xtable(t(data.frame(model4_summary2)))

rm(iGP_gqs)

#---- Loo, model 5:
load(file = store_model_out_model5)

log_lik_5 <- loo::extract_log_lik(mGP_gqs, merge_chains = FALSE, parameter_name = "log_like_age")
r_eff_5   <- loo::relative_eff(exp(log_lik_5),    cores = nChains)
loo_5     <- loo::loo(log_lik_5, r_eff = r_eff_5, cores = nChains)
waic_5    <- loo::waic(log_lik_5)

print(loo_5)
print(waic_5)

# Export to latex format:
model5_summary <- c(paste0(round(waic_5$estimates[3,1],1), " (",  round(waic_5$estimates[3,2],1), ")"),
                    paste0(round(loo_5$estimates[3,1],1),  " (",  round(loo_5$estimates[3,2],1), ")"),
                    paste0(round(loo_5$estimates[1,1],1),  " (",  round(loo_5$estimates[1,2],1), ")"),
                    round(loo_5$estimates[3,1] - 2*loo_5$estimates[2,1], 1),
                    paste0(round(loo_5$estimates[2,1],1),  " (",  round(loo_5$estimates[2,2],1), ")")
)
names(model5_summary)<- c("WAIC", "LooIC", "ELPD", "Adequacy", "Complexity")

xtable(t(data.frame(model5_summary)))

model5_summary2 <- c(round(loo_5$estimates[3,1] - 2*loo_5$estimates[2,1], 1),
                     paste0(round(loo_5$estimates[2,1],1),  " (",  round(loo_5$estimates[2,2],1), ")"),
                     paste0(round(loo_5$estimates[3,1],1),  " (",  round(loo_5$estimates[3,2],1), ")")
)
names(model5_summary2)<- c("Adequacy", "Complexity", "LooIC")

xtable(t(data.frame(model5_summary2)))

rm(mGP_gqs)

#---- Loo, model 6:
load(file = store_model_out_model6)

log_lik_6 <- loo::extract_log_lik(mxGP_gqs, merge_chains = FALSE, parameter_name = "log_like_age")
r_eff_6   <- loo::relative_eff(exp(log_lik_6),    cores = nChains)
loo_6     <- loo::loo(log_lik_6, r_eff = r_eff_6, cores = nChains)
waic_6    <- loo::waic(log_lik_6)

print(loo_6)
print(waic_6)

# Export to latex format:
model6_summary <- c(paste0(round(waic_6$estimates[3,1],1), " (",  round(waic_6$estimates[3,2],1), ")"),
                    paste0(round(loo_6$estimates[3,1],1),  " (",  round(loo_6$estimates[3,2],1), ")"),
                    paste0(round(loo_6$estimates[1,1],1),  " (",  round(loo_6$estimates[1,2],1), ")"),
                    round(loo_6$estimates[3,1] - 2*loo_6$estimates[2,1], 1),
                    paste0(round(loo_6$estimates[2,1],1),  " (",  round(loo_6$estimates[2,2],1), ")")
)
names(model6_summary)<- c("WAIC", "LooIC", "ELPD", "Adequacy", "Complexity")

xtable(t(data.frame(model6_summary)))

model6_summary2 <- c(round(loo_6$estimates[3,1] - 2*loo_6$estimates[2,1], 1),
                     paste0(round(loo_6$estimates[2,1],1),  " (",  round(loo_6$estimates[2,2],1), ")"),
                     paste0(round(loo_6$estimates[3,1],1),  " (",  round(loo_6$estimates[3,2],1), ")")
)
names(model6_summary2)<- c("Adequacy", "Complexity", "LooIC")

xtable(t(data.frame(model6_summary2)))

rm(mxGP_gqs)

# #---- Loo, model 7:
# load(file = store_model_out_model7)
# 
# log_lik_7 <- loo::extract_log_lik(fit_ENG_mxGP, merge_chains = FALSE, parameter_name = "log_likelihood")
# r_eff_7   <- loo::relative_eff(exp(log_lik_7),    cores = nChains)
# loo_7     <- loo::loo(log_lik_7, r_eff = r_eff_7, cores = nChains)
# waic_7    <- loo::waic(log_lik_7)
# 
# print(loo_7)
# print(waic_7)
# 
# # Export to latex format:
# model7_summary <- c(paste0(round(waic_7$estimates[3,1],1), " (",  round(waic_7$estimates[3,2],1), ")"),
#                     paste0(round(loo_7$estimates[3,1],1),  " (",  round(loo_7$estimates[3,2],1), ")"),
#                     paste0(round(loo_7$estimates[1,1],1),  " (",  round(loo_7$estimates[1,2],1), ")"),
#                     round(loo_7$estimates[3,1] - 2*loo_7$estimates[2,1], 1),
#                     paste0(round(loo_7$estimates[2,1],1),  " (",  round(loo_7$estimates[2,2],1), ")")
# )
# names(model7_summary)<- c("WAIC", "LooIC", "ELPD", "Adequacy", "Complexity")
# 
# xtable(t(data.frame(model7_summary)))
# 
# model7_summary2 <- c(round(loo_7$estimates[3,1] - 2*loo_7$estimates[2,1], 1),
#                      paste0(round(loo_7$estimates[2,1],1),  " (",  round(loo_7$estimates[2,2],1), ")"),
#                      paste0(round(loo_7$estimates[3,1],1),  " (",  round(loo_7$estimates[3,2],1), ")")
# )
# names(model7_summary2)<- c("Adequacy", "Complexity", "LooIC")
# 
# xtable(t(data.frame(model7_summary2)))
# 
# rm(fit_ENG_mxGP)

#---- Model comparison:
comp2 <- loo::loo_compare(loo_1, loo_2, loo_3, loo_4, loo_5, loo_6)
looic_comparison <- round(cbind(comp2[,1] * -2, comp2[,2] * 2), 1)
looic_comparison <- looic_comparison[order(rownames(looic_comparison)),]

looic_comparison[,2] * 1.96

#---- Report LooIC, Adequacy, Complexity, Delta LooIC and compare models by their information criteria:
# https://discourse.mc-stan.org/t/where-did-loo-difference-estimates-go/8178/18

cbind(comp2[,1] * -2, comp2[,2] * 2)

# Approximate 95% CI:
round(-10.1 + c(-1,1) * 1.96 * 1.1 , 2)

round(-2.3  + c(-1,1) * 1.96 * 0.9 , 2)

comp3 <- loo::loo_compare(loo_1, loo_3)
print(comp3)

# Approximate 95% CI:
round(-1.6 + c(-1,1) * 1.96 * 1.7 , 2)

#---- Compare the three models with WAIC:
# comp_waic <- loo::loo_compare(waic_1, waic_2)
# print(comp_waic)

# Interpretation: When that difference, elpd_diff, is positive then 
# the expected predictive accuracy for the second model is higher. 
# A negative elpd_diff favors the first model.

# Pairwise comparisons between each model and the model with the largest ELPD 
# (the model in the first row). For this reason the elpd_diff column will always 
# have the value 0 in the first row (i.e., the difference between the preferred model 
# and itself) and negative values in subsequent rows for the remaining models.