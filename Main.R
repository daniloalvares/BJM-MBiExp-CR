rm(list = ls(all = TRUE))

# Required sources
library(DBI)
library(dplyr)
library(tidyr)
library(stringr)
library(ggplot2)
library(gridExtra)
library(ggh4x)
library(rstan)
library(statmod)
library(loo)
library(MASS)
library(emdbook)
library(adaptMCMC)
library(pracma)
library(parallel)
library(survival)
library(survminer)
library(pec)
library(timeROC)
library(scales)

# ====================================================================== #
# Loading data and applying inclusion/exclusion criteria
# - Combine Kappa/Lambda-FLC into a single biomarker
# - Cut-off: 99th percentile for biomarkers
# - Eliminate jumps between LoTs
source("functions/1 - Load and Clear Data.R")
# ====================================================================== #
# Splitting the data into training and test sets
prop <- 0.2 # Proportion for test
source("functions/2 - Splitting Data.R")
# ====================================================================== #
# Standardising continuous variables (in log scale) and imputing
# their missing data via mean imputation (zero)
source("functions/3 - Standardisation and Imputation.R")
# ====================================================================== #
# Creating longitudinal and survival data
source("functions/4 - Long and Short Formats.R")
# ====================================================================== #
# Building joint models in Stan considering joint specification and
# corrected two-stage approaches
source("functions/5 - Stan Models.R")
# ====================================================================== #
# Auxiliary functions
source("functions/6 - Auxiliary Functions.R")
# ====================================================================== #


# ====================================================================== #
#       JOINT MODELLING FOR ALL LOTS USING THE JOINT SPECIFICATION       #
# ====================================================================== #
# LoT 1
fit_JM1 <- fit_jm(data = trainlot1, LoT = 1)
# LoT 2
fit_JM2 <- fit_jm(data = trainlot2, LoT = 2)
# LoT 3
fit_JM3 <- fit_jm(data = trainlot3, LoT = 3)
# LoT 4
fit_JM4 <- fit_jm(data = trainlot4, LoT = 4)

# Efficiency (effective sample size): Ideally, > 100
# Convergence (R-hat): Ideally, < 1.05
check_jm_ess_rhat(fit_JM1$fit, fit_JM2$fit, fit_JM3$fit, fit_JM4$fit)

# POSTERIOR SUMMARIES
# LoT 1
round(posterior_summary_jm(fit_JM1$fit, type = "long"), 3)
round(posterior_summary_jm(fit_JM1$fit, type = "surv", LoT = 1), 3)
# LoT 2
round(posterior_summary_jm(fit_JM2$fit, type = "long"), 3)
round(posterior_summary_jm(fit_JM2$fit, type = "surv", LoT = 2), 3)
# LoT 3
round(posterior_summary_jm(fit_JM3$fit, type = "long"), 3)
round(posterior_summary_jm(fit_JM3$fit, type = "surv", LoT = 3), 3)
# LoT 4
round(posterior_summary_jm(fit_JM4$fit, type = "long"), 3)
round(posterior_summary_jm(fit_JM4$fit, type = "surv", LoT = 4), 3)
# ====================================================================== #


# ====================================================================== #
#       JOINT MODELLING FOR ALL LOTS USING THE CORRECTED TWO-STAGE       #
# ====================================================================== #
# LoT 1
fit_LONG1_M <- fit_biexp(data = trainlot1, biomarker = "M-spike")
fit_LONG1_F <- fit_biexp(data = trainlot1, biomarker = "FLC")
fit_CR1 <- fit_surv(data = trainlot1, fit_LONG_M = fit_LONG1_M$fit, 
                    fit_LONG_F = fit_LONG1_F$fit, LoT = 1)
# LoT 2
fit_LONG2_M <- fit_biexp(data = trainlot2, biomarker = "M-spike")
fit_LONG2_F <- fit_biexp(data = trainlot2, biomarker = "FLC")
fit_CR2 <- fit_surv(data = trainlot2, fit_LONG_M = fit_LONG2_M$fit, 
                    fit_LONG_F = fit_LONG2_F$fit, LoT = 2)
# LoT 3
fit_LONG3_M <- fit_biexp(data = trainlot3, biomarker = "M-spike")
fit_LONG3_F <- fit_biexp(data = trainlot3, biomarker = "FLC")
fit_CR3 <- fit_surv(data = trainlot3, fit_LONG_M = fit_LONG3_M$fit, 
                    fit_LONG_F = fit_LONG3_F$fit, LoT = 3)
# LoT 4
fit_LONG4_M <- fit_biexp(data = trainlot4, biomarker = "M-spike")
fit_LONG4_F <- fit_biexp(data = trainlot4, biomarker = "FLC")
fit_Surv4 <- fit_surv(data = trainlot4, fit_LONG_M = fit_LONG4_M$fit, 
                      fit_LONG_F = fit_LONG4_F$fit, LoT = 4)

# Efficiency (effective sample size): Ideally, > 100
# Convergence (R-hat): Ideally, < 1.05
check_long_ess_rhat(fit_LONG1_M$fit, fit_LONG1_F$fit, 
                    fit_LONG2_M$fit, fit_LONG2_F$fit,
                    fit_LONG3_M$fit, fit_LONG3_F$fit, 
                    fit_LONG4_M$fit, fit_LONG4_F$fit)
check_surv_ess_rhat(fit_CR1$fit, fit_CR2$fit, fit_CR3$fit, fit_Surv4$fit)

# POSTERIOR SUMMARIES
# LoT 1
round(posterior_summary_ts(fit_LONG1_M$fit, type = "long"), 3)
round(posterior_summary_ts(fit_LONG1_F$fit, type = "long"), 3)
round(posterior_summary_ts(fit_CR1$fit, type = "surv", LoT = 1), 3)
# LoT 2
round(posterior_summary_ts(fit_LONG2_M$fit, type = "long"), 3)
round(posterior_summary_ts(fit_LONG2_F$fit, type = "long"), 3)
round(posterior_summary_ts(fit_CR2$fit, type = "surv", LoT = 2), 3)
# LoT 3
round(posterior_summary_ts(fit_LONG3_M$fit, type = "long"), 3)
round(posterior_summary_ts(fit_LONG3_F$fit, type = "long"), 3)
round(posterior_summary_ts(fit_CR3$fit, type = "surv", LoT = 3), 3)
# LoT 4
round(posterior_summary_ts(fit_LONG4_M$fit, type = "long"), 3)
round(posterior_summary_ts(fit_LONG4_F$fit, type = "long"), 3)
round(posterior_summary_ts(fit_Surv4$fit, type = "surv", LoT = 4), 3)
# ====================================================================== #


# ====================================================================== #
#                 FITTED AVERAGE LONGITUDINAL TRAJECTORY                 #
# ====================================================================== #
plot_avg_trajectory(fit_JM1, fit_JM2, fit_JM3, fit_JM4,
                    fit_LONG1_M, fit_LONG2_M, fit_LONG3_M, fit_LONG4_M,
                    fit_LONG1_F, fit_LONG2_F, fit_LONG3_F, fit_LONG4_F)
# ====================================================================== #
