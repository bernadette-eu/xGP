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

main_path       <- "C://Users//lbour//OneDrive//Desktop1//xgp_project//GP multitype epidemic England//4_Stan_Outputs//"
store_model_out <- paste0(main_path, "SEEIIR_renewal_multitype_XGBM", ".RData") #SEEIIR_renewal_multitype_XGBM_v5C
#covData_out     <- paste0(main_path, "SEEIIR_renewal_multitype_XGBM_v5C", ".RData") #SEEIIR_renewal_multitype_XGBM_v5C

#---- Set system locale to English:
Sys.setlocale("LC_ALL", "English")
`%nin%` <- Negate(`%in%`)

load(file = store_model_out)
#load(file = covData_out)

posts_mcmc <- rstan::extract(nuts_fit_1) # Posterior draws   # fit_ENG_xGP

###########################################################################
#
# Age distribution - adjusted figures of estimated infections 
# (Table 2, Ward 2021):
#
###########################################################################
# x = age_distr_ENG
# user_AgeGrp = lookup_table_ENG
aggregate_react_cumcases <- function(x, 
                                     user_AgeGrp) 
{
  options(dplyr.summarise.inform = FALSE)
  
  if (nrow(user_AgeGrp) != nrow(x)) 
    stop("The mapped age group labels do not correspond to the age group labels of the aggregated age distribution matrix.\n")
  
  REACT2_cumul_cases <- data.frame(AgeGrp = c("0-14", "15-44", "45-64", "65-74", "74-100"),
                                   Mean   = c(0, 1536 * 1e+03, 895 * 1e+03 ,181 * 1e+03, 166 * 1e+03),
                                   Lower  = c(0, 1437 * 1e+03, 837 * 1e+03 ,153 * 1e+03, 131 * 1e+03),
                                   Upper  = c(0, 1635 * 1e+03, 953 * 1e+03 ,209 * 1e+03, 201 * 1e+03) )
  
  REACT2_cumul_cases$AgeGrpStart <- sapply(1:nrow(REACT2_cumul_cases), function(x) {
    min(as.numeric(strsplit(REACT2_cumul_cases$AgeGrp, "-")[[x]]))
  })
  
  REACT2_cumul_cases$AgeGrpEnd <- sapply(1:nrow(REACT2_cumul_cases), function(x) {
    max(as.numeric(strsplit(REACT2_cumul_cases$AgeGrp, "-")[[x]]))
  })
  
  temp_x <- x
  temp_x$AgeGrp <- gsub("\\+", "-100", temp_x$AgeGrp)
  temp_x$AgeGrpEnd <- sapply(1:nrow(temp_x), function(x) {
    max(as.numeric(strsplit(temp_x$AgeGrp, "-")[[x]]))
  })
  
  temp_x$Group_mapping <- user_AgeGrp$Mapping
  
  # Calculating temp_x$Mean, temp_x$Lower, temp_x$Upper within the same loop was not working, replacing the
  # Lower and Upper fields with the values of the oldest age group. Workaround = 3 separate loops.
  if ("AgeGrpStart" %nin% colnames(temp_x)) 
    temp_x$AgeGrpStart <- sapply(1:nrow(temp_x), function(x) {
      min(as.numeric(strsplit(temp_x$AgeGrp, "-")[[x]]))
    })
  
  for (i in 1:nrow(temp_x)) {
    for (j in 1:nrow(REACT2_cumul_cases)) {
      if ( (temp_x$AgeGrpStart[i] >= REACT2_cumul_cases$AgeGrpStart[j]) & 
           (temp_x$AgeGrpEnd[i]   <= REACT2_cumul_cases$AgeGrpEnd[j]) ) 
        
        temp_x$Mean[i]  <- REACT2_cumul_cases$Mean[j]
      
    }
  }
  
  for (i in 1:nrow(temp_x)) {
    for (j in 1:nrow(REACT2_cumul_cases)) {
      if ( (temp_x$AgeGrpStart[i] >= REACT2_cumul_cases$AgeGrpStart[j]) & 
           (temp_x$AgeGrpEnd[i]   <= REACT2_cumul_cases$AgeGrpEnd[j]) ) 
        
        temp_x$Lower[i] <- REACT2_cumul_cases$Lower[j]
      
    }
  }
  
  for (i in 1:nrow(temp_x)) {
    for (j in 1:nrow(REACT2_cumul_cases)) {
      if ( (temp_x$AgeGrpStart[i] >= REACT2_cumul_cases$AgeGrpStart[j]) & 
           (temp_x$AgeGrpEnd[i]   <= REACT2_cumul_cases$AgeGrpEnd[j]) ) 
        
        temp_x$Upper[i] <- REACT2_cumul_cases$Upper[j]
      
    }
  }
  
  output <- temp_x %>% 
    as.data.frame() %>%
    dplyr::group_by(Group_mapping) %>% 
    dplyr::mutate(PopPerc = prop.table(PopTotal), 
                  Mean    = sum(Mean  * PopPerc),
                  Lower   = sum(Lower * PopPerc),
                  Upper   = sum(Upper * PopPerc)
    ) %>% 
    dplyr::select(dplyr::one_of(c("Group_mapping", 
                                  "Mean", "Lower", "Upper"))) %>% 
    dplyr::group_by(Group_mapping) %>% 
    dplyr::slice(1) %>% 
    as.data.frame()
  
  return(output)
}

age_distr_ENG         <- age_distribution(country = "United Kingdom", year = 2020)
age_distr_ENG$PopTotal<- age_distr_ENG$PopTotal*(sum(cov_data$n_pop)/sum(age_distr_ENG$PopTotal))

# Mapping of the age distribution:
lookup_table_ENG <- data.frame(Initial = age_distr_ENG$AgeGrp,
                               Mapping = c(rep("0-19",  4),
                                           rep("20-39", 4), 
                                           rep("40-59", 4),
                                           rep("60+",   4)))

REACT2_adjusted_cumul_cases <- aggregate_react_cumcases(x           = age_distr_ENG, 
                                                        user_AgeGrp = lookup_table_ENG) 
colnames(REACT2_adjusted_cumul_cases)[1] <- "Group"

###########################################################################
#
# Data engineering
#
###########################################################################

#---- Country and analysis period:
country    <- "England"
start_date <- "2020-03-02"
end_date   <- "2020-09-27"
dates     <- seq(as.Date(start_date), as.Date(end_date), by = "days")

###########################################
# Estimated infections (per group)
###########################################

#---- Transform model outputs to the appropriate format:
data_cases_cols      <- c("Date", "Group", "New_Cases", "median", "low", "high")
data_cases           <- data.frame(matrix(ncol = length(data_cases_cols), nrow = 0))
colnames(data_cases) <- data_cases_cols

for (i in 1:cov_data$A){

  if( "y_deaths" %in% names(cov_data)) {
    group_id <- colnames(cov_data$y_deaths)
  } else { 
    group_id <- colnames(cov_data$y_data)
  }# End if
  
  dt_cases_age_grp  <- data.frame(Date  = dates,
                                  Group = rep(group_id[i], length(dates)))
  
  word <- c("E_casesAge", "E_casesdByAge", "cases")
  
  if ( any( word[1] %in% names(posts_mcmc)) ){
    fit_cases <- posts_mcmc$E_casesAge[,,i] 
  } else if ( any( word[2] %in% names(posts_mcmc)) ){
    fit_cases <- posts_mcmc$E_casesdByAge[,,i] 
  } else {
    fit_cases <- posts_mcmc$cases[,i,] 
  }
  
  # Add quantiles from the model outputs:
  dt_cases_age_grp$median <- apply(fit_cases, 2, median)
  dt_cases_age_grp$low    <- apply(fit_cases, 2, quantile, probs = c(0.025))
  dt_cases_age_grp$high   <- apply(fit_cases, 2, quantile, probs = c(0.975))
  dt_cases_age_grp$low25  <- apply(fit_cases, 2, quantile, probs = c(0.25))
  dt_cases_age_grp$high75 <- apply(fit_cases, 2, quantile, probs = c(0.75))

  data_cases  <- rbind(data_cases, dt_cases_age_grp)

}# End for

# #---- Check the transformed model outputs:
# table(data_cases$Group)

#---- Remove redundant objects:
rm(i,
   dt_cases_age_grp)

###########################################################################
#
# Cumulative group-specific cases
#
###########################################################################

# NOTE: Please edit the following accordingly, in order to not plot the observed cumulative infections if they are not available.

# Locate duplicates, after rounding to 0 dps:
data_cases %>% 
  mutate_if(is.numeric, round, 0) %>% 
  dplyr::group_by(Group) %>% 
  filter(duplicated(median) | duplicated(median, fromLast = TRUE))

cumul_data <- 
  data_cases %>% 
  mutate_if(is.numeric, round, 0) %>% 
  dplyr::select(-dplyr::one_of(c("Date"))) %>% 
  dplyr::group_by(Group) %>% 
  dplyr::mutate(Date      = dates, 
                median    = cumsum(median * !duplicated(median)),
                low       = cumsum(low    * !duplicated(low)),
                high      = cumsum(high   * !duplicated(high)),
                low25     = cumsum(low25  * !duplicated(low25)),
                high75    = cumsum(high75 * !duplicated(high75))
  )

#---- Store the plot:
# Setting: 7 times 12 inches, pdf format
# # Name: "SEEIIR_ENG_Age_model4_Model_building_Experiment3_CumulCases_REACT"
cumul_data_sub                  <- subset(cumul_data, Group != "0-19")
REACT2_adjusted_cumul_cases_sub <- subset(REACT2_adjusted_cumul_cases, Group != "0-19")

# Save as pdf 13x10, then convert to .eps:
validationPlot <- 
  ggplot(cumul_data_sub,
       aes(Date      = Date,
           New_Cases = median)) +
  facet_wrap(. ~ Group, 
             scales = "free_y", 
             ncol   = 1,
             strip.position = "top",
             shrink = T) +
  geom_line(aes(x     = Date,
                y     = median,
                color = "Median"),
            size = 1.3) +
  geom_ribbon(aes(x    = Date,
                  ymin = low25,
                  ymax = high75,
                  fill = "50% CrI"),
              alpha = 0.5) +
  geom_ribbon(aes(x    = Date,
                  ymin = low,
                  ymax = high,
                  fill = "95% CrI"),
              alpha = 0.5)  +
  labs(x = "Epidemiological Date",
       y = "Cumulative infections") +
  scale_x_date(date_breaks       = "1 month", 
               date_minor_breaks = "1 month") + 
  scale_y_continuous(labels = scales::scientific, n.breaks = 3) +
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
        )+
  geom_vline(xintercept = as.Date("2020-07-17"), colour = "black" ) +
  geom_hline(data = REACT2_adjusted_cumul_cases_sub,
             aes(yintercept = Lower),
             colour   = "black", 
             linetype = "dashed") +
  geom_hline(data = REACT2_adjusted_cumul_cases_sub,
             aes(yintercept = Upper),
             colour   = "black", 
             linetype = "dashed") +
  geom_hline(data = REACT2_adjusted_cumul_cases_sub,
             aes(yintercept = Mean),
             colour   = "black", 
             linetype = "solid") 




###########################################################################
#
# Cumulative overall cases
#
###########################################################################
fit_cases        <- posts_mcmc$E_cases
median_fit_cases <- apply(fit_cases, 2, median)
low_fit_cases    <- apply(fit_cases, 2, quantile, probs = c(0.025)) 
high_fit_cases   <- apply(fit_cases, 2, quantile, probs = c(0.975)) 
low25_fit_cases  <- apply(fit_cases, 2, quantile, probs = c(0.25)) 
high75_fit_cases <- apply(fit_cases, 2, quantile, probs = c(0.75))

cases_plot_data        <- data.frame(Date = dates)
cases_plot_data$median <- median_fit_cases
cases_plot_data$low    <- low_fit_cases
cases_plot_data$high   <- high_fit_cases
cases_plot_data$low25  <- low25_fit_cases
cases_plot_data$high75 <- high75_fit_cases

cumulCasesOverall <- 
  cases_plot_data %>% 
  mutate_if(is.numeric, round, 0) %>% 
  dplyr::mutate(median    = cumsum(median * !duplicated(median)),
                low       = cumsum(low    * !duplicated(low)),
                high      = cumsum(high   * !duplicated(high)),
                low25     = cumsum(low25  * !duplicated(low25)),
                high75    = cumsum(high75 * !duplicated(high75))
  )

REACT2_cumul_cases <- data.frame(AgeGrp = c("0-14", "15-44", "45-64", "65-74", "74-100"),
                                 Mean   = c(0, 1536 * 1e+03, 895 * 1e+03 ,181 * 1e+03, 166 * 1e+03),
                                 Lower  = c(0, 1437 * 1e+03, 837 * 1e+03 ,153 * 1e+03, 131 * 1e+03),
                                 Upper  = c(0, 1635 * 1e+03, 953 * 1e+03 ,209 * 1e+03, 201 * 1e+03) )
react2AggrCases <- apply(REACT2_cumul_cases[,-1], 2, sum)
react2AggrCases <- data.frame(Mean  = react2AggrCases[1],
                              Lower = react2AggrCases[2],
                              Upper = react2AggrCases[3])

#---- Store the plot:
# Setting: 7 times 12 inches, pdf format
# # Name: "SEEIIR_ENG_Age_model4_Model_building_Experiment3_CumulCases_REACT"

ggplot(cumulCasesOverall,
       aes(Date      = Date,
           New_Cases = median)) +
  geom_line(aes(x     = Date,
                y     = median,
                color = "Median"),
            size = 1.3) +
  geom_ribbon(aes(x    = Date,
                  ymin = low25,
                  ymax = high75,
                  fill = "50% CrI"),
              alpha = 0.5) +
  geom_ribbon(aes(x    = Date,
                  ymin = low,
                  ymax = high,
                  fill = "95% CrI"),
              alpha = 0.5)  +
  labs(x = "Epidemiological Date",
       y = "Cumulative infections") +
  scale_x_date(date_breaks       = "1 month", 
               date_minor_breaks = "1 month") + 
  scale_y_continuous(labels = scales::scientific, n.breaks = 5) +
  scale_fill_manual(values = c("50% CrI" = "gray20",
                               "95% CrI"  = "gray40")) +
  scale_colour_manual(name   = '',
                      values = c('Median' = "black")) +          
  theme_bw() +
  theme(panel.spacing    = unit(0.2,"cm"),
        axis.text.x      = element_text(angle = 45, hjust = 1),
        axis.title.x     = element_text(size = 14, face = "bold"),
        axis.title.y     = element_text(size = 14, face = "bold"),
        legend.position  = "bottom",
        legend.title     = element_blank(),
        legend.box       = "vertical", 
        legend.margin    = margin())+
  geom_vline(xintercept = as.Date("2020-07-17"), colour = "black" ) +
  geom_hline(data = react2AggrCases,
             aes(yintercept = Lower),
             colour   = "black", 
             linetype = "dashed") +
  geom_hline(data = react2AggrCases,
             aes(yintercept = Upper),
             colour   = "black", 
             linetype = "dashed") +
  geom_hline(data = react2AggrCases,
             aes(yintercept = Mean),
             colour   = "black", 
             linetype = "solid") 




