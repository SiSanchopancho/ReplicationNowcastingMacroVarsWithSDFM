# SPDX-License-Identifier: GPL-3.0-or-later #
#
# Copyright (C) 2026 Domenic Franjic
#
# This file is part of ReplicationNowcastingMacroVarsWithSDFM.
#
# ReplicationNowcastingMacroVarsWithSDFM is free software: you can redistribute
# it and/or modify it under the terms of the GNU General Public License as
# published by the Free Software Foundation, either version 3 of the License,
# or (at your option) any later version.
#
# ReplicationNowcastingMacroVarsWithSDFM is distributed in the hope that it
# will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
# of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with ReplicationNowcastingMacroVarsWithSDFM. If not, see <https://www.gnu.org/licenses/>.
#

# Empirical application #

library(rstudioapi)
library(zoo)
library(lubridate)
library(BVAR)
library(readxl)
library(alfred)
library(TwoStepSDFM) # Download from Github https://github.com/SiSanchopancho/TwoStepSDFM
library(sandwich)
library(lmtest)
library(tidyr)
library(dplyr)
library(car)
library(stargazer)
library(murphydiagram)
library(dynlm)
library(modelsummary)
library(ggplot2)
setwd(dirname(getActiveDocumentContext()$path))
rm(list = ls())

# # Start data download and clean up loop #
# 
# # Note: Run only once to download and pre-process the data
#
# # loop over FRED-MD vintage files in correspoding directory
# quarterly_series_trans <- read.csv("./fred-qd_2024m12.csv")[2, -1] # This is only used to get the correct transformation code for the GDP series
# monthly_file_list <- list.files(path="./all_fred_md_vintage_dir", pattern="*.csv", full.names = TRUE, recursive = FALSE)
# pb = txtProgressBar(min = 0, max = length(monthly_file_list), initial = 0, style = 3)
# step <- 0
# Realisation <- c()
# for(file_name in monthly_file_list){
# 
#   setTxtProgressBar(pb, step)
# 
#   # load and clean monthly data set
#   fred_md_raw <- read.csv(file_name)
#   trans <- as.integer(fred_md_raw[1, -1])
#   fred_md_clean <- fred_md_raw[-1, -1]
#   dates <- mdy(fred_md_raw[-1, 1])
#   rownames(fred_md_clean) <- dates
#   fred_md_clean <- na.locf(fred_md_clean, fromLast = TRUE, na.rm = FALSE) # impute nas at beginning of sample but keep at the end
# 
#   # transform the data using the bvar package function
#   fred_md <- fred_transform(fred_md_clean, "fred_md", codes = trans, na.rm = FALSE)
#   fred_md <- na.locf(fred_md, fromLast = TRUE, na.rm = FALSE)
#   fred_md_ts <- ts(fred_md, start = c(year(as.Date(rownames(fred_md)[1])), month(as.Date(rownames(fred_md)[1]))), frequency = 12)
#   fred_md_zoo <- as.zoo(fred_md_ts)
# 
#   # Add quarterly GDP vintages extracted from AL-FRED
#   current_vintage <- quantdates::LastDayOfMonth(date = rownames(fred_md)[dim(fred_md)[1]])
# 
#   # Try downloading GDP until it works as it sometimes does not pull the data from the webiste
#   current_gdp_series <- NULL
#   while(is.null(current_gdp_series)){
#     current_gdp_series <- get_alfred_series("GDPC1", "GDP",
#                                             realtime_start = current_vintage,
#                                             realtime_end = current_vintage,
#                                             api_key = # SET YOUR API KEY HERE
# )[, c(1, 3)]
#   }
# 
#   # Clean up the GDP series and make it stationary
#   quarterly_data_clean <- na.locf(current_gdp_series, fromLast = TRUE, na.rm = FALSE)
#   quarterly_data <- fred_transform(quarterly_data_clean[, -1, drop = FALSE], "fred_qd",
#                                    codes = quarterly_series_trans[which(colnames(quarterly_series_trans) %in% "GDPC1")],
#                                    na.rm = FALSE)
#   Realisation[step + 1] <- quarterly_data[dim(quarterly_data)[1], ]
#   fred_qd_ts <- ts(quarterly_data, start = c(year(as.Date(current_gdp_series[1, 1])), quarter(as.Date(current_gdp_series[1, 1]))),
#                    frequency = 4)
#   monthly_time_index <- seq(from = as.yearmon(time(fred_qd_ts)[1]),
#                             to   = as.yearmon(time(fred_qd_ts)[length(fred_qd_ts)]) + 2/12,
#                             by   = 1/12)
#   fred_qd_zoo <- zoo(rep(coredata(fred_qd_ts), each = 3), monthly_time_index)
#   fred_md_zoo <- merge.zoo(fred_md_zoo, fred_qd_zoo)
#   colnames(fred_md_zoo)[dim(fred_md_zoo)[2]] <- "GDP"
# 
#   # re-order and cut data for later usage
#   fred_md_zoo <- fred_md_zoo[, c("GDP",
#                                  colnames(fred_md))]
#   vintage <- na.locf(fred_md_zoo, fromLast = TRUE, na.rm = FALSE)
#   vintage <- window(vintage, start = as.yearmon("1979-01"), end = current_vintage)
#   vintage_df <- fortify.zoo(vintage)[, -1]
# 
#   # create a delay and frequency vector
#   del_m <- colSums(is.na(vintage))
#   freq_m <- c(4, rep(12, length(trans)))
# 
#   # save the data
#   directory_name <- paste0("./clean_data/fredgdp_real_time/vintage_", rownames(fred_md)[dim(fred_md)[1]], "/")
#   dir.create(directory_name)
#   write.table(round(vintage_df, 15), paste0(directory_name, "data.csv"),
#               row.names = FALSE, col.names = FALSE, sep = ",")
#   t_out <- dim(vintage_df)[1]
#   n_out <- dim(vintage_df)[2]
#   write.table(as.Date(time(vintage)), paste0(directory_name, "dates.csv"), row.names = FALSE, col.names = FALSE, sep = ",")
#   write.table(del_m, paste0(directory_name, "delay.csv"), row.names = FALSE, col.names = FALSE, sep = ",")
#   write.table(freq_m, paste0(directory_name, "frequency.csv"), row.names = FALSE, col.names = FALSE, sep = ",")
#   write.table(colnames(vintage_df), paste0(directory_name, "names.csv"), row.names = FALSE, col.names = FALSE, sep = ",")
#   write.table(month(as.Date(rownames(fred_md)[dim(fred_md)[1]])) %% 3, paste0(directory_name, "month_of_qtr_indicator.csv"), row.names = FALSE, col.names = FALSE, sep = ",")
#   fc_window_start <- which(as.Date(time(vintage)) %in% as.Date(rownames(fred_md)[dim(fred_md)[1]])) - 1
#   fc_window_end <- which(as.Date(time(vintage)) %in% as.Date(rownames(fred_md)[dim(fred_md)[1]])) - 1
#   write.table(seq(fc_window_start, fc_window_end, 3), paste0(directory_name, "fc_dates.csv"), row.names = FALSE, col.names = FALSE, sep = ",")
# 
#   step <- step + 1
# }
# close(pb)
# write.table(Realisation, "./clean_data/gdp_realisations.csv", row.names = FALSE, col.names = FALSE, sep = ",")

# End data download and clean up loop #

# prediction block #
rm(list = ls())

# Set up result containers
vintage_file_list <- list.files(path = "./clean_data/fredgdp_real_time", full.names = TRUE, recursive = FALSE)
first_dates <- as.Date(unlist(read.table(paste0(vintage_file_list[1], "/dates.csv"), sep = ",", header = FALSE)))
last_dates <- as.Date(unlist(read.table(paste0(vintage_file_list[length(vintage_file_list)], "/dates.csv"), sep = ",", header = FALSE)))
start_date <- ymd(as.Date(first_dates[length(first_dates)]))
end_date <- ymd(as.Date(last_dates[length(last_dates)]))
no_of_predictions <- time_length(interval(start_date, end_date), "month")
Realisation <- read.table("./clean_data/gdp_realisations.csv", sep = ",", header = FALSE)
nowcasts <- as.zoo(ts(matrix(NA, no_of_predictions, 5),
                      start = c(year(start_date), month(start_date)),
                      end = c(year(end_date), month(end_date)),
                      frequency = 12))
colnames(nowcasts) <- c("SDFM(CV)", "SDFM(BIC)", "DFM", "Realisation", "Target Date")
nowcasts$Realisation <- Realisation[, 1]
month_three_ind <- which(month(time(nowcasts)) %% 3 == 0)
month_two_ind <- which(month(time(nowcasts)) %% 3 == 2)
month_one_ind <- which(month(time(nowcasts)) %% 3 == 1)
nowcasts$`Target Date` <- as.yearqtr(time(nowcasts))

# start prediction loop over all vintages #

vintage_ind <- 1
no_of_factors_old <- 0
for(vintage_name in vintage_file_list){ # Loop only over the last months of the quarters

  # Load and prepare data
  # data <- read.table(paste0(vintage_name, "/data.csv"), sep = ",", header = FALSE)
  data <- read.table(paste0(vintage_name, "/data.csv"), sep = ",", header = FALSE)
  print(dim(data))
  data_scaled <- scale(data)
  scaling <- attributes(data_scaled)$`scaled:scale`
  centering <- attributes(data_scaled)$`scaled:center`
  delay <- unlist(read.table(paste0(vintage_name, "/delay.csv"), sep = ",", header = FALSE))
  frequency <- unlist(read.table(paste0(vintage_name, "/frequency.csv"), sep = ",", header = FALSE))
  dates <- as.Date(unlist(read.table(paste0(vintage_name, "/dates.csv"), sep = ",", header = FALSE)))
  
  # Correct the quarterly delay to always be equal to three
  qtrly_delay_correction_index <- which(delay %% 3 != 0 & frequency == 4)
  if(length(qtrly_delay_correction_index) != 0){
    delay[qtrly_delay_correction_index] <- delay[qtrly_delay_correction_index] + 3 - (delay[qtrly_delay_correction_index] %% 3)
  }
  names <- unlist(read.table(paste0(vintage_name, "/names.csv"), sep = ",", header = FALSE))
  dates <- as.Date(unlist(read.table(paste0(vintage_name, "/dates.csv"), sep = ",", header = FALSE)))
  no_of_obs <- dim(data)[1]
  data_zoo <- as.zoo(ts(data_scaled, end = c(year(dates[no_of_obs]), month(dates[no_of_obs])),
                        frequency = 12))
  colnames(data_zoo) <- names
  variables_of_interest <- which(colnames(data_zoo) == "GDP")
  
  # Estimate the number of factors:
  # We choose a range between 1 and 7 factor for a reasonable interval and a good
  #   testing power. The confidence level for rejecting the null is based on the
  #   original paper of Onatski (2009)
  min_no_of_factors <- 2
  max_no_of_factors <- 7
  no_of_factors <- noOfFactors(data_zoo[, which(frequency == 12)], min_no_of_factors, max_no_of_factors, 0.01)$no_of_factors
  
  # Re-estimate the number of factors if no_of_factors == max_no_of_factors -1
  while(no_of_factors == max_no_of_factors - 1 && max_no_of_factors < 21){
    min_no_of_factors <- min_no_of_factors + 1
    max_no_of_factors <- max_no_of_factors + 1
    no_of_factors <- noOfFactors(data_zoo[, which(frequency == 12)], min_no_of_factors, max_no_of_factors, 0.01)$no_of_factors
  }
  
  # Fit the dense model
  dense_fit <- TwoStepSDFM::nowcast(data = data_zoo, variables_of_interest = variables_of_interest,
                                    max_fcast_horizon = 4, delay = delay,
                                    selected = NULL, sparse = FALSE,
                                    frequency = frequency, no_of_factors = no_of_factors,
                                    conv_crit = 1e-8, decorr_errors = FALSE,
                                    max_ar_lag_order = 1, max_predictor_lag_order = 1
  )
  predictions_DFM <- na.omit(dense_fit$Forecasts$`Fcast GDP`)
  nowcasts$DFM[vintage_ind] <- centering[variables_of_interest] + predictions_DFM[1] * scaling[variables_of_interest]
  
  # SDFM model tuning
  # The model is only tuned if the number of factors changes
  if(no_of_factors_old != no_of_factors){
    no_of_factors_old <- no_of_factors
    cv_results <- crossVal(data = data_zoo, variable_of_interest = variables_of_interest,
                           fcast_horizon = 0, delay = delay, frequency = frequency,
                           no_of_factors = no_of_factors, seed = 16102025, min_ridge_penalty = 0.01,
                           max_ridge_penalty = 10, cv_repititions = 3, cv_size = 100 * no_of_factors,
                           lasso_penalty_type = "selected", min_max_penalty = c(10, sum(frequency == 12)),
                           parallel = TRUE, no_of_cores = floor(parallel::detectCores()/2), max_factor_lag_order = 10,
                           comp_null = 1e-10, conv_crit = 1e-8,   max_ar_lag_order = 1,
                           max_predictor_lag_order = 1)
  }
  
  # SDFM nowcasting
  nowcast_cv <- TwoStepSDFM::nowcast(data = data_zoo, variables_of_interest = variables_of_interest,
                                     max_fcast_horizon = 4, delay = delay,
                                     selected = cv_results$CV$`Min. CV`[3:(3 + no_of_factors - 1)],
                                     frequency = frequency, no_of_factors = no_of_factors,
                                     conv_crit = 1e-8, max_ar_lag_order = 1,
                                     max_predictor_lag_order = 1
  )
  predictions_cv <- na.omit(nowcast_cv$Forecasts$`Fcast GDP`)
  nowcasts$`SDFM(CV)`[vintage_ind] <- centering[variables_of_interest] + predictions_cv[1] * scaling[variables_of_interest]
  
  nowcast_bic <- TwoStepSDFM::nowcast(data = data_zoo, variables_of_interest = variables_of_interest,
                                      max_fcast_horizon = 4, delay = delay,
                                      selected = cv_results$BIC$`Min. BIC`[3:(3 + no_of_factors - 1)],
                                      frequency = frequency, no_of_factors = no_of_factors,
                                      conv_crit = 1e-8, max_ar_lag_order = 1,
                                      max_predictor_lag_order = 1
  )
  predictions_bic <- na.omit(nowcast_bic$Forecasts$`Fcast GDP`)
  nowcasts$`SDFM(BIC)`[vintage_ind] <- centering[variables_of_interest] + predictions_bic[1] * scaling[variables_of_interest]
  
  # Print some stuff for tracking
  cat(paste0("\n#################################################\n",
             "Executed vintage No.: ", vintage_ind, "; ", time(nowcasts)[vintage_ind], ".\n",
             "\n_________________________________________________\n",
             "_________________________________________________\n",
             "Model       SDFM(CV)       SDFM(BIC)       DFM\n",
             "_________________________________________________\n",
             "MSNE    \n",
             " Month 1    ", sprintf("%07.4f", mean((nowcasts$`SDFM(CV)` - nowcasts$Realisation)[month_one_ind]^2, na.rm = TRUE)),
             "        ", sprintf("%07.4f", mean((nowcasts$`SDFM(BIC)` - nowcasts$Realisation)[month_one_ind]^2, na.rm = TRUE)),
             "        ", sprintf("%07.4f", mean((nowcasts$`DFM` - nowcasts$Realisation)[month_one_ind]^2, na.rm = TRUE)), "\n",
             " Month 2    ", sprintf("%07.4f", mean((nowcasts$`SDFM(CV)` - nowcasts$Realisation)[month_two_ind]^2, na.rm = TRUE)),
             "        ", sprintf("%07.4f", mean((nowcasts$`SDFM(BIC)` - nowcasts$Realisation)[month_two_ind]^2, na.rm = TRUE)),
             "        ", sprintf("%07.4f", mean((nowcasts$`DFM` - nowcasts$Realisation)[month_two_ind]^2, na.rm = TRUE)), "\n",
             " Month 3    ", sprintf("%07.4f", mean((nowcasts$`SDFM(CV)` - nowcasts$Realisation)[month_three_ind]^2, na.rm = TRUE)),
             "        ", sprintf("%07.4f", mean((nowcasts$`SDFM(BIC)` - nowcasts$Realisation)[month_three_ind]^2, na.rm = TRUE)),
             "        ", sprintf("%07.4f", mean((nowcasts$`DFM` - nowcasts$Realisation)[month_three_ind]^2, na.rm = TRUE)), "\n",
             "_________________________________________________\n",
             "##################################################\n"
  )
  )
  
  vintage_ind <- 1 + vintage_ind
}

# End prediction loop over all vintages #

# Statistical test block #

# Extract loss function values and period indicators
sq_fcst_error_cv <- as.vector(nowcasts$Realisation - nowcasts$`SDFM(CV)`)[month_three_ind]^2
sq_fcst_error_dense <- as.vector(nowcasts$Realisation - nowcasts$DFM)[month_three_ind]^2
sq_fcst_diff <- as.vector(nowcasts$DFM - nowcasts$`SDFM(CV)`)[month_three_ind]^2
pre_corona_ind <- (round(as.numeric(time(nowcasts)), 3) < round(2020 + (3 - 1)/12, 3))[month_three_ind]
corona_ind <- (round(as.numeric(time(nowcasts)), 3) %in% round(seq(2020 + (3 - 1)/12, 2021 + (3 - 1)/12, by = 1/12), 3))[month_three_ind]
post_corona_ind <- (round(as.numeric(time(nowcasts)), 3) > round(2021 + (3 - 1)/12, 3))[month_three_ind]
dot_com_ind <- (round(as.numeric(time(nowcasts)), 3) %in% round(seq(2001 + (3 - 1)/12, 2001 + (12 - 1)/12, by = 1/12), 3))[month_three_ind]
financial_ind <- (round(as.numeric(time(nowcasts)), 3) %in% round(seq(2007 + (12 - 1)/12, 2009 + (6 - 1)/12, by = 1/12), 3))[month_three_ind]
corona_and_after <- corona_ind | post_corona_ind

# Clark-West approximate normal test for equal predicitive ability of nested models (Clark, T. E., & West, K. D. (2007). Approximately normal tests for equal predictive accuracy in nested models. Journal of econometrics, 138(1), 291-311.)
adJ_loss_differential <- sq_fcst_error_dense - sq_fcst_error_cv + sq_fcst_diff
clark_west_fit <- lm(adJ_loss_differential ~ 1)
model_summary <- summary(clark_west_fit)
model_summary
regular_t_stat <- model_summary$coefficients[3]
regular_t_stat
regular_p_val <- pt(regular_t_stat, df = clark_west_fit$df.residual, lower.tail = FALSE)
regular_p_val

nw_var_cov <- sandwich::vcovHAC(clark_west_fit, lag = 4)
nw_t_stat <- model_summary$coefficients[1] / sqrt(nw_var_cov)
nw_t_stat
nw_p_val <- pt(nw_t_stat, df = clark_west_fit$df.residual, lower.tail = FALSE)
nw_p_val

# Compute MSNE during different periods #

full_msne_cv <- mean(sq_fcst_error_cv)
round(full_msne_cv, 3)
full_msne_dense <- mean(sq_fcst_error_dense)
round(full_msne_dense, 3)
full_reduction <-  1 - full_msne_cv / full_msne_dense
round(100 * full_reduction, 3)

pre_corona_msne_cv <- mean(sq_fcst_error_cv[pre_corona_ind])
round(pre_corona_msne_cv, 3)
pre_corona_msne_dense <- mean(sq_fcst_error_dense[pre_corona_ind])
round(pre_corona_msne_dense, 3)
pre_corona_reduction <-  1 - pre_corona_msne_cv / pre_corona_msne_dense
round(100 * pre_corona_reduction, 3)

corona_msne_cv <- mean(sq_fcst_error_cv[corona_ind])
round(corona_msne_cv, 3)
corona_msne_dense <- mean(sq_fcst_error_dense[corona_ind])
round(corona_msne_dense, 3)
corona_reduction <-  1 - corona_msne_cv / corona_msne_dense
round(100 * corona_reduction, 3)

post_corona_msne_cv <- mean(sq_fcst_error_cv[post_corona_ind])
round(post_corona_msne_cv, 3)
post_corona_msne_dense <- mean(sq_fcst_error_dense[post_corona_ind])
round(post_corona_msne_dense, 3)
post_corona_reduction <-  1 - post_corona_msne_cv / post_corona_msne_dense
round(100 * post_corona_reduction, 3)

dotcom_msne_reduction <-  1 - mean(sq_fcst_error_cv[dot_com_ind]) / mean(sq_fcst_error_dense[dot_com_ind])
round(100 * dotcom_msne_reduction, 3)
financial_msne_reduction <-  1 - mean(sq_fcst_error_cv[financial_ind]) / mean(sq_fcst_error_dense[financial_ind])
round(100 * financial_msne_reduction, 3)

# loss_differential <- ts(sq_fcst_error_dense - sq_fcst_error_cv)
loss_differential <- ts(sq_fcst_error_dense - sq_fcst_error_cv)
corona_dummy <- as.numeric(corona_ind)
post_corona_dummy <- as.numeric(post_corona_ind)
pre_corona_dummy <- as.numeric(pre_corona_ind)

# Giacomini-White-Style regressions (Giacomini, R., & White, H. (2006). Tests of conditional predictive ability. Econometrica, 74(6), 1545-1578.) #
# Note: Due to the expanding/recursive nowcasting scheme the p-values and t-statistics are asymptotically
#   not valid. We still employ the regression type analysis here to conduct a more sophisticated analysis
#   of our regimes in the given regimes. Also, we think that the the t-statistics and p-value are
#   somewhat insightful, as we do not believe that an increase in the sample size from about 83 to 105
#   quarterly observations is enough to trigger the asymptotics.

# Giacomini white tests with period dummies only
gw_test_reg <- dynlm(loss_differential ~ -1
                     + pre_corona_dummy
                     + corona_dummy
                     + post_corona_dummy
)
var_cov <- vcovHAC(gw_test_reg)
coeftest(gw_test_reg, df = gw_test_reg$df.residual, vcov. = var_cov)
Wald_stat <- matrix(gw_test_reg$coefficients, nrow = 1) %*% solve(var_cov) %*% matrix(gw_test_reg$coefficients, ncol = 1)
Wald_stat
1 - pchisq(Wald_stat, length(gw_test_reg$coefficients))
linearHypothesis(gw_test_reg, c("pre_corona_dummy = 0", "corona_dummy = 0", 
                                "post_corona_dummy = 0"), vcov. = var_cov)

# Giacomini white tests with period dummies and lags
forecast::auto.arima(loss_differential, ic = "aic") # Use 2 lags according to aic
gw_test_reg_lags <- dynlm(loss_differential ~ -1
                          + pre_corona_dummy
                          + corona_dummy
                          + post_corona_dummy
                          + L(loss_differential, 1)
                          + L(loss_differential, 2)
)
var_cov_lag <- vcovHAC(gw_test_reg_lags)
coeftest(gw_test_reg_lags, df = gw_test_reg_lags$df.residual, vcov. = var_cov_lag)
Wald_stat_lag <- matrix(gw_test_reg_lags$coefficients, nrow = 1) %*% solve(var_cov_lag) %*% matrix(gw_test_reg_lags$coefficients, ncol = 1)
Wald_stat_lag
1 - pchisq(Wald_stat_lag, length(gw_test_reg_lags$coefficients))
linearHypothesis(gw_test_reg_lags, c("pre_corona_dummy = 0", "corona_dummy = 0", 
                                     "post_corona_dummy = 0", "L(loss_differential, 1) = 0",
                                     "L(loss_differential, 2) = 0"), vcov. = var_cov_lag)


# Giacomini white tests with period dummies and lags
no_crisis_ind <- !corona_ind & !post_corona_ind & !dot_com_ind& !financial_ind
no_crisis_dummy <- as.numeric(no_crisis_ind)
forecast::auto.arima(loss_differential, ic = "aic") # Use 2 lags according to aic
dot_com_dummy <- as.numeric(dot_com_ind)
financial_dummy <- as.numeric(financial_ind)
gw_test_reg_lags_all_crisis <- dynlm(loss_differential ~ -1
                                     + no_crisis_dummy
                                     + corona_dummy
                                     + post_corona_dummy
                                     + dot_com_dummy
                                     + financial_dummy
                                     + L(loss_differential, 1)
                                     + L(loss_differential, 2)
)
var_cov_lag_all_crisis <- vcovHAC(gw_test_reg_lags_all_crisis)
coeftest(gw_test_reg_lags_all_crisis, df = gw_test_reg_lags_all_crisis$fitted.values, vcov. = var_cov_lag_all_crisis)
Wald_stat_lag_all_crisis <- matrix(gw_test_reg_lags_all_crisis$coefficients, nrow = 1) %*% solve(var_cov_lag_all_crisis) %*% matrix(gw_test_reg_lags_all_crisis$coefficients, ncol = 1)
Wald_stat_lag_all_crisis
1 - pchisq(Wald_stat_lag_all_crisis, length(gw_test_reg_lags_all_crisis$coefficients))
linearHypothesis(gw_test_reg_lags_all_crisis, c("no_crisis_dummy = 0", "corona_dummy = 0", 
                                                "post_corona_dummy = 0", "dot_com_dummy = 0",
                                                "financial_dummy = 0", "L(loss_differential, 1) = 0",
                                                "L(loss_differential, 2) = 0"), vcov. = var_cov_lag_all_crisis)
# Giacomini white tests with covid and no-covid
no_covid_dummy <- post_corona_dummy + pre_corona_dummy
gw_test_reg_no_covid <- dynlm(loss_differential ~ -1
                              + no_covid_dummy
                              + corona_dummy
)
var_cov_no_covid <- vcovHAC(gw_test_reg_no_covid)
coeftest(gw_test_reg_no_covid, df = gw_test_reg_no_covid$df.residual, vcov. = var_cov_no_covid)
Wald_stat_no_covid <- matrix(gw_test_reg_no_covid$coefficients, nrow = 1) %*% solve(var_cov_no_covid) %*% matrix(gw_test_reg_no_covid$coefficients, ncol = 1)
Wald_stat_no_covid
1 - pchisq(Wald_stat_no_covid, length(gw_test_reg_no_covid$coefficients))
linearHypothesis(gw_test_reg_no_covid, c("no_covid_dummy = 0", "corona_dummy = 0"), vcov. = var_cov_no_covid)

# Giacomini white tests with covid and no-covid with lags
no_covid_dummy <- post_corona_dummy + pre_corona_dummy
gw_test_reg_no_covid_lags <- dynlm(loss_differential ~ -1
                                   + no_covid_dummy
                                   + corona_dummy
                                   + L(loss_differential, 1)
                                   + L(loss_differential, 2)
)
var_cov_no_covid_lags <- vcovHAC(gw_test_reg_no_covid_lags)
coeftest(gw_test_reg_no_covid_lags, df = gw_test_reg_no_covid_lags$df.residual, vcov. = var_cov_no_covid_lags)
Wald_stat_no_covid_lags <- matrix(gw_test_reg_no_covid_lags$coefficients, nrow = 1) %*% solve(var_cov_no_covid_lags) %*% matrix(gw_test_reg_no_covid_lags$coefficients, ncol = 1)
Wald_stat_no_covid_lags
1 - pchisq(Wald_stat_no_covid, length(gw_test_reg_no_covid$coefficients))
linearHypothesis(gw_test_reg_no_covid_lags, c("no_covid_dummy = 0", "corona_dummy = 0"), vcov. = var_cov_no_covid_lags)


# Plotting
results_combined_zoo <- merge(nowcasts$`SDFM(CV)`, nowcasts$`DFM`, nowcasts$Realisation)
results_combined <- cbind(as.data.frame(as.yearqtr(time(results_combined_zoo)[month_three_ind])), 
                          as.data.frame(coredata(results_combined_zoo[month_three_ind, ])))
colnames(results_combined) <- c("Dates", "SDFM", "DFM", "Realisation")
results_long <- pivot_longer(results_combined[, 1:3], cols = c(SDFM, DFM), names_to = "Model", values_to = "Value")

# Plot skelleton
nowcast_plots <- ggplot() +
  geom_hline(yintercept = 0) +
  geom_bar(data = results_combined, 
           aes(x = Dates, y = Realisation, fill = "Realisation"), 
           stat = "identity", alpha = 0.5, width = 0.1) +
  geom_line(data = results_long, 
            aes(x = Dates, y = Value, color = Model, group = Model, linetype = Model), 
            size = 1) +
  labs(x = "Dates", y = expression(Delta ~ "log gdp")) +
  theme_minimal() 

# Full evaluation period plot
full_nowcast_plots <- nowcast_plots +
  scale_x_yearqtr(format = "%Y-Q%q", expand = c(0, 0)) +
  scale_fill_manual(values = c("Realisation" = "#85C0F9")) +
  scale_color_manual(values = c("SDFM" = "#000000", "DFM" = "#F5793A")) +
  scale_linetype_manual(values = c("DFM" = "twodash", "SDFM" = "solid")) +
  theme(
    axis.text.x = element_text(angle = 30, vjust = 1, hjust = 1),
    text = element_text(size = 30),
    axis.line = element_line(color = "black"), 
    legend.title = element_blank()
  ) +
  theme(legend.position = "none")
full_nowcast_plots

# Pre-Covid plot
nowcast_plot_pre_covid <- nowcast_plots +
  scale_x_yearqtr(format = "%Y-Q%q", expand = c(0, 0), 
                  limits = as.yearqtr(c("1999 Q3", "2019 Q4"))) +
  scale_y_continuous(limits = c(-2, 2)) +
  scale_fill_manual(values = c("Realisation" = "#85C0F9")) +
  scale_color_manual(values = c("SDFM" = "#000000", "DFM" = "#F5793A")) +
  scale_linetype_manual(values = c("DFM" = "twodash", "SDFM" = "solid")) +
  theme(
    axis.text.x = element_text(angle = 30, vjust = 1, hjust = 1),
    text = element_text(size = 50),
    axis.line = element_line(color = "black"), 
    legend.title = element_blank()
  ) +
  theme(legend.position = "none")
nowcast_plot_pre_covid

# Post-Covid plot
nowcast_plot_post_covid <- nowcast_plots +
  scale_x_yearqtr(format = "%Y-Q%q", expand = c(0, 0), 
                  limits = as.yearqtr(c("2020 Q4", "2025 Q4"))) +
  scale_y_continuous(limits = c(-2, 2)) +
  scale_fill_manual(values = c("Realisation" = "#85C0F9")) +
  scale_color_manual(values = c("SDFM" = "#000000", "DFM" = "#F5793A")) +
  scale_linetype_manual(values = c("DFM" = "twodash", "SDFM" = "solid")) +
  theme(
    axis.text.x = element_text(angle = 30, vjust = 1, hjust = 1),
    text = element_text(size = 50),
    axis.line = element_line(color = "black"), 
    legend.title = element_blank()
  ) +
  theme(legend.position = "none")
nowcast_plot_post_covid

