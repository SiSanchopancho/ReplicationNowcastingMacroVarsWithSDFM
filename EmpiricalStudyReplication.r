# SPDX-License-Identifier: GPL-3.0-or-later #
#
# Copyright (C) 2024-2026 Domenic Franjic
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
install.packages("./TwoStepSDFM_0.3.0.3.tar.gz", source = TRUE)
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
library(xtable)
library(readr)
library(foreach)
library(doParallel)
library(doSNOW)
library(parallel)
setwd(dirname(getActiveDocumentContext()$path))
rm(list = ls())


# # # Start data download and clean up loop #
# 
# # Note: Run only once to download and pre-process the data
# 
# # loop over FRED-MD vintage files in correspoding directory
# quarterly_series_trans <- read.csv("./raw_data/historical vintages of fred-qd 2018-05 to 2024-12/fred-qd_2024m12.csv")[2, -1]
# monthly_file_list <- list.files(path="./raw_data/historical-vintages-of-fred-md-2015-01-to-2024-12", pattern="*.csv", full.names = TRUE, recursive = FALSE)
# pb = txtProgressBar(min = 0, max = length(monthly_file_list), initial = 0, style = 3)
# step <- 0
# Realisation <- as.data.frame(matrix(NA, length(monthly_file_list), 3))
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
#     )[, c(1, 3)]
#   }
#   
#   # Clean up the GDP series and make it stationary
#   quarterly_data_clean <- na.locf(current_gdp_series, fromLast = TRUE, na.rm = FALSE)
#   quarterly_data <- fred_transform(quarterly_data_clean[, -1, drop = FALSE], "fred_qd",
#                                    codes = quarterly_series_trans[which(colnames(quarterly_series_trans) %in% "GDPC1")],
#                                    na.rm = FALSE)
#   fred_qd_ts <- as.zoo(ts(quarterly_data, start = c(year(as.Date(current_gdp_series[1, 1])), quarter(as.Date(current_gdp_series[1, 1]))),
#                           frequency = 4))
#
#   Realisation[step + 1, 1] <- fred_qd_ts[dim(fred_qd_ts)[1], ]
#   Realisation[step + 1, 2] <- as.character(current_gdp_series$date[dim(current_gdp_series)[1]])
#   Realisation[step + 1, 3] <- as.character(current_vintage)
#   
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
Rcpp::sourceCpp("RHelper/ARDL.cpp")
source("./RHelper/MidasForecaster.R")

# Set up result containers
vintage_file_list <- list.files(path = "./clean_data/fredgdp_real_time", full.names = TRUE, recursive = FALSE)
first_dates <- as.Date(unlist(read.table(paste0(vintage_file_list[1], "/dates.csv"), sep = ",", header = FALSE)))
last_dates <- as.Date(unlist(read.table(paste0(vintage_file_list[length(vintage_file_list)], "/dates.csv"), sep = ",", header = FALSE)))
start_date <- ymd(as.Date(first_dates[length(first_dates)]))
end_date <- ymd(as.Date(last_dates[length(last_dates)]))
no_of_predictions <- time_length(interval(start_date, end_date), "month")
Realisation <- read.table("./clean_data/gdp_realisations.csv", sep = ",", header = FALSE)

nowcasts <- as.zoo(ts(matrix(NA, no_of_predictions, 7),
                      start = c(year(start_date), month(start_date)),
                      end = c(year(end_date), month(end_date)),
                      frequency = 12))
colnames(nowcasts) <- c("SDFM(CV)", "SDFM(BIC)", "DFM", "DFM EM", "SDFM EM", "Realisation", "Target Date")
nowcasts$`Target Date` <- as.yearqtr(time(nowcasts))
nowcast_ind <- which(nowcasts$`Target Date` %in%  as.numeric(as.yearqtr(as.Date(Realisation$V2))))
reali_nowcast_ind <- which(as.numeric(as.yearqtr(as.Date(Realisation$V2))) %in% nowcasts$`Target Date`)
nowcasts$Realisation[nowcast_ind] <- Realisation[reali_nowcast_ind, 1]
month_three_ind <- which(month(time(nowcasts)) %% 3 == 0)
month_two_ind <- which(month(time(nowcasts)) %% 3 == 2)
month_one_ind <- which(month(time(nowcasts)) %% 3 == 1)
nowcasts$Realisation[c(month_one_ind, month_two_ind)] <- NA
nowcasts$Realisation <- na.locf(nowcasts$Realisation, fromLast = TRUE, na.rm = FALSE)

one_step <- as.zoo(ts(matrix(NA, no_of_predictions, 7),
                      start = c(year(start_date), month(start_date)),
                      end = c(year(end_date), month(end_date)),
                      frequency = 12))
colnames(one_step) <- c("SDFM(CV)", "SDFM(BIC)", "DFM", "DFM EM", "SDFM EM", "Realisation", "Target Date")
one_step$`Target Date` <- nowcasts$`Target Date` + 0.25
osh_ind <- which(one_step$`Target Date` %in%  as.numeric(as.yearqtr(as.Date(Realisation$V2))))
reali_osh_ind <- which(as.numeric(as.yearqtr(as.Date(Realisation$V2))) %in% one_step$`Target Date`)
one_step$Realisation[osh_ind] <- Realisation[reali_osh_ind, 1]
one_step$Realisation[c(month_one_ind, month_two_ind)] <- NA
one_step$Realisation <- na.locf(one_step$Realisation, fromLast = TRUE, na.rm = FALSE)

# start prediction loop over all vintages #
vintage_ind <- 3 # Start at the end of the first quarter
vintage_name <- vintage_file_list[[vintage_ind]]
step_size <- 1
no_of_factors_old <- 0
cv_grid_size <- 100
density_reps <- 100
marginal_density_nowcasts <- matrix(NA, length(vintage_file_list), cv_grid_size * density_reps)
set.seed(07072026)
for(vintage_name in vintage_file_list[seq(vintage_ind, length(vintage_file_list), step_size)]){
  
  # Load and prepare data
  dates <- as.Date(unlist(read.table(paste0(vintage_name, "/dates.csv"), sep = ",", header = FALSE)))
  covid_time <- (dates[length(dates)] >= as.Date("2020-03-01"))
  data <- read.table(paste0(vintage_name, "/data.csv"), sep = ",", header = FALSE)
  delay <- unlist(read.table(paste0(vintage_name, "/delay.csv"), sep = ",", header = FALSE))
  frequency <- unlist(read.table(paste0(vintage_name, "/frequency.csv"), sep = ",", header = FALSE))
  names <- unlist(read.table(paste0(vintage_name, "/names.csv"), sep = ",", header = FALSE))
  data_scaled <- scale(data)
  scaling <- attributes(data_scaled)$`scaled:scale`
  centering <- attributes(data_scaled)$`scaled:center`
  no_of_obs <- dim(data)[1]
  data_zoo <- as.zoo(ts(data_scaled, end = c(year(dates[no_of_obs]), month(dates[no_of_obs])),
                        frequency = 12))
  colnames(data_zoo) <- names
  variables_of_interest <- which(colnames(data_zoo) == "GDP")
  end_of_qtr_ind <- month(dates[length(dates)]) %% 3 == 0
  
  if(end_of_qtr_ind && !covid_time){
    # if(end_of_qtr_ind){
    # Estimate the number of factors:
    # We choose a range between 1 and 7 factor for a reasonable interval and a good
    #   testing power. The confidence level for rejecting the null is based on the
    #   original paper of Onatski (2009)
    min_no_of_factors <- 2
    max_no_of_factors <- 8
    test_res <- noOfFactors(data_zoo[, which(frequency == 12)], min_no_of_factors, max_no_of_factors, 0.01)
    test_res_plot <- plot(test_res)
    test_res_plot$`Eigen Value Plot Test Procedure`
    
    # Re-estimate the number of factors if no_of_factors == max_no_of_factors -1
    while(test_res$test$no_of_factors == max_no_of_factors - 1 && max_no_of_factors < 21){
      min_no_of_factors <- min_no_of_factors + 1
      max_no_of_factors <- max_no_of_factors + 1
      no_of_factors <- noOfFactors(data_zoo[, which(frequency == 12)], min_no_of_factors, max_no_of_factors, 0.01)
    }
    test_res_plot <- plot(test_res)
    test_res_plot$`Eigen Value Plot Test Procedure`
    no_of_factors_test <- test_res$test$no_of_factors
  }
  print(test_res)
  
  # Fit the dense model
  dense_fit <- TwoStepSDFM::nowcast(data = data_zoo, variables_of_interest = variables_of_interest,
                                    max_fcast_horizon = 2, delay = delay,
                                    selected = NULL, sparse = FALSE,
                                    frequency = frequency, no_of_factors = no_of_factors_test,
                                    decorr_errors = FALSE, max_ar_lag_order = 1, max_predictor_lag_order = 1
  )
  predictions_DFM_test <- na.omit(dense_fit$Forecasts$`Fcast GDP`)
  nowcasts$`DFM`[vintage_ind] <- centering[variables_of_interest] + predictions_DFM_test[1] * scaling[variables_of_interest]
  one_step$`DFM`[vintage_ind] <- centering[variables_of_interest] + predictions_DFM_test[2] * scaling[variables_of_interest]
  
  # SDFM model tuning nowcast and one step ahead
  # The model is only tuned if the number of factors changes outside of covid
  if(!covid_time){
    if(no_of_factors_old != no_of_factors_test){
      no_of_factors_old <- no_of_factors_test
      fcast_horizon <- ifelse(month(dates[length(dates)]) %% 3 != 0, 1, 0)
      cv_results_test <- crossVal(data = data_zoo, variable_of_interest = variables_of_interest,
                                  fcast_horizon = fcast_horizon, delay = delay, frequency = frequency,
                                  no_of_factors = no_of_factors_test, seed = 16102025, min_ridge_penalty = 100,
                                  max_ridge_penalty = 1000, cv_repetitions = 3, cv_size = cv_grid_size * no_of_factors_test,
                                  lasso_penalty_type = "selected", min_max_penalty = c(10, sum(frequency == 12)),
                                  parallel = TRUE, no_of_cores = floor(parallel::detectCores()/2), max_factor_lag_order = 10,
                                  comp_null = 1e-10, max_ar_lag_order = 1,
                                  max_predictor_lag_order = 1)
    }
  }
  
  # SDFM nowcasting
  nowcast_cv <- TwoStepSDFM::nowcast(data = data_zoo, variables_of_interest = variables_of_interest,
                                     max_fcast_horizon = 2, delay = delay,
                                     selected = cv_results_test$CV$`Min. CV`[3:(3 + no_of_factors_test - 1)],
                                     frequency = frequency, no_of_factors = no_of_factors_test,
                                     max_ar_lag_order = 1, max_predictor_lag_order = 1,
                                     ridge_penalty = cv_results_test$CV$`Min. CV`[2]
  )
  predictions_cv_test <- na.omit(nowcast_cv$Forecasts$`Fcast GDP`)
  nowcasts$`SDFM(CV)`[vintage_ind] <- centering[variables_of_interest] + predictions_cv_test[1] * scaling[variables_of_interest]
  one_step$`SDFM(CV)`[vintage_ind] <- centering[variables_of_interest] + predictions_cv_test[2] * scaling[variables_of_interest]
  
  nowcast_bic_test <- TwoStepSDFM::nowcast(data = data_zoo, variables_of_interest = variables_of_interest,
                                           max_fcast_horizon = 2, delay = delay,
                                           selected = cv_results_test$BIC$`Min. BIC`[3:(3 + no_of_factors_test - 1)],
                                           frequency = frequency, no_of_factors = no_of_factors_test,
                                           max_ar_lag_order = 1, max_predictor_lag_order = 1,
                                           ridge_penalty = cv_results_test$BIC$`Min. BIC`[2]
  )
  predictions_bic_test <- na.omit(nowcast_bic_test$Forecasts$`Fcast GDP`)
  nowcasts$`SDFM(BIC)`[vintage_ind] <- centering[variables_of_interest] + predictions_bic_test[1] * scaling[variables_of_interest]
  one_step$`SDFM(BIC)`[vintage_ind] <- centering[variables_of_interest] + predictions_bic_test[2] * scaling[variables_of_interest]

  # Dense EM nowcasting
  dense_em_fit <- sparseDFM::sparseDFM(data_zoo[, which(frequency == 12)], no_of_factors_test,
                                       standardize = FALSE, alg = "EM")
  factor_delay <- rep((3 - month(index(data_zoo)[nrow(data_zoo)]) %% 3) %% 3, no_of_factors_test)
  dense_em_forecast_data <- matrix(NA, dim(data_zoo)[1] + factor_delay[1], no_of_factors_test + 1)
  dense_em_forecast_data[1:dim(data_zoo)[1], 1] <- data_zoo[, c("GDP")]
  dense_em_forecast_data[1:dim(data_zoo)[1], 2:(no_of_factors_test + 1)] <- dense_em_fit$state$factors
  em_data_delay <- colSums(is.na(dense_em_forecast_data))
  pred <- MIDASForecaster(target_variables = t(dense_em_forecast_data[, 1, drop = FALSE]),
                          monthly_predictors = t(dense_em_forecast_data[, 2:(no_of_factors_test + 1)]),
                          target_variable_delay = em_data_delay[1],
                          predictors_delay = em_data_delay[2:(no_of_factors_test + 1)],
                          lag_estim_criterion = "BIC",
                          max_fcast_horizon = 2,
                          max_ar_lag_order = 1,
                          max_predictor_lag_order = 1)
  nowcasts$`DFM EM`[vintage_ind] <- centering[variables_of_interest] + pred$`Avg. Point Forecast`[1] * scaling[variables_of_interest]
  one_step$`DFM EM`[vintage_ind] <- centering[variables_of_interest] + pred$`Avg. Point Forecast`[2] * scaling[variables_of_interest]

  # Dense EM nowcasting
  sparse_em_fit <- sparseDFM::sparseDFM(data_zoo[, which(frequency == 12)], no_of_factors_test,
                                        standardize = FALSE, alg = "EM-sparse")
  factor_delay <- rep((3 - month(index(data_zoo)[nrow(data_zoo)]) %% 3) %% 3, no_of_factors_test)
  sparse_em_forecast_data <- matrix(NA, dim(data_zoo)[1] + factor_delay[1], no_of_factors_test + 1)
  sparse_em_forecast_data[1:dim(data_zoo)[1], 1] <- data_zoo[, c("GDP")]
  sparse_em_forecast_data[1:dim(data_zoo)[1], 2:(no_of_factors_test + 1)] <- sparse_em_fit$state$factors
  sparse_em_data_delay <- colSums(is.na(sparse_em_forecast_data))
  sparse_em_pred <- MIDASForecaster(target_variables = t(sparse_em_forecast_data[, 1, drop = FALSE]),
                                    monthly_predictors = t(sparse_em_forecast_data[, 2:(no_of_factors_test + 1)]),
                                    target_variable_delay = sparse_em_data_delay[1],
                                    predictors_delay = sparse_em_data_delay[2:(no_of_factors_test + 1)],
                                    lag_estim_criterion = "BIC",
                                    max_fcast_horizon = 2,
                                    max_ar_lag_order = 1,
                                    max_predictor_lag_order = 1)

  nowcasts$`SDFM EM`[vintage_ind] <- centering[variables_of_interest] + sparse_em_pred$`Avg. Point Forecast`[1] * scaling[variables_of_interest]
  one_step$`SDFM EM`[vintage_ind] <- centering[variables_of_interest] + sparse_em_pred$`Avg. Point Forecast`[2] * scaling[variables_of_interest]

  # Extract nowcast hyper-parameter distribution
  selected_cand <- matrix(round(runif(no_of_factors_test * cv_grid_size * density_reps, 10, sum(frequency == 12))), ncol = no_of_factors_test)
  ridge_cand <- exp(runif(cv_grid_size * density_reps, log(100), log(1000)))
  pb <- txtProgressBar(max = cv_grid_size * density_reps, style = 3)
  progressFunc <- function(n) setTxtProgressBar(pb, n)
  opts <- list(progress = progressFunc)
  cl <- makeCluster(parallel::detectCores() / 2 - 1)
  registerDoSNOW(cl)
  nowcast_dist <- foreach(sample = 1:(cv_grid_size * density_reps), .options.snow = opts, .combine = c) %dopar% {
    fcast_dist <- TwoStepSDFM::nowcast(data = data_zoo, variables_of_interest = variables_of_interest,
                                       max_fcast_horizon = 2, delay = delay,
                                       selected = selected_cand[sample, ],
                                       frequency = frequency, no_of_factors = no_of_factors_test,
                                       max_ar_lag_order = 1, max_predictor_lag_order = 1,
                                       ridge_penalty = ridge_cand[sample]
    )
    predictions_fcast_dist <- na.omit(fcast_dist$Forecasts$`Fcast GDP`)
    zoo::coredata(centering[variables_of_interest] + predictions_fcast_dist[1] * scaling[variables_of_interest])
  }
  close(pb)
  stopCluster(cl)
  marginal_density_nowcasts[vintage_ind, ] <- nowcast_dist
  plot(density(marginal_density_nowcasts[vintage_ind, ]))
  abline(v = nowcasts$`SDFM(CV)`[vintage_ind], col = "red")
  abline(v = nowcasts$`DFM`[vintage_ind], col = "green")
  abline(v = nowcasts$`SDFM EM`[vintage_ind], col = "blue")
  abline(v = nowcasts$`DFM EM`[vintage_ind], col = "black")

  # Print some stuff for tracking
  cat(paste0("\n##################################################################################################\n",
             "Executed vintage No.: ", vintage_ind, "; ", time(nowcasts)[vintage_ind], ".\n",
             "\n__________________________________________________________________________________________________\n",
             "__________________________________________________________________________________________________\n",
             "Model       SDFM(CV)      SDFM(BIC)      DFM         DFM EM             SDFM EM\n",
             "__________________________________________________________________________________________________\n",
             "MSNE    \n",
             " Month 1    ", sprintf("%07.4f", mean((nowcasts$`SDFM(CV)` - nowcasts$Realisation)[month_one_ind]^2, na.rm = TRUE)),
             "        ", sprintf("%07.4f", mean((nowcasts$`SDFM(BIC)` - nowcasts$Realisation)[month_one_ind]^2, na.rm = TRUE)),
             "        ", sprintf("%07.4f", mean((nowcasts$`DFM` - nowcasts$Realisation)[month_one_ind]^2, na.rm = TRUE)),
             "        ", sprintf("%07.4f", mean((nowcasts$`DFM EM` - nowcasts$Realisation)[month_one_ind]^2, na.rm = TRUE)),
             "        ", sprintf("%07.4f", mean((nowcasts$`SDFM EM` - nowcasts$Realisation)[month_one_ind]^2, na.rm = TRUE)), "\n",
             " Month 2    ", sprintf("%07.4f", mean((nowcasts$`SDFM(CV)` - nowcasts$Realisation)[month_two_ind]^2, na.rm = TRUE)),
             "        ", sprintf("%07.4f", mean((nowcasts$`SDFM(BIC)` - nowcasts$Realisation)[month_two_ind]^2, na.rm = TRUE)),
             "        ", sprintf("%07.4f", mean((nowcasts$`DFM` - nowcasts$Realisation)[month_two_ind]^2, na.rm = TRUE)),
             "        ", sprintf("%07.4f", mean((nowcasts$`DFM EM` - nowcasts$Realisation)[month_two_ind]^2, na.rm = TRUE)),
             "        ", sprintf("%07.4f", mean((nowcasts$`SDFM EM` - nowcasts$Realisation)[month_two_ind]^2, na.rm = TRUE)), "\n",
             " Month 3    ", sprintf("%07.4f", mean((nowcasts$`SDFM(CV)` - nowcasts$Realisation)[month_three_ind]^2, na.rm = TRUE)),
             "        ", sprintf("%07.4f", mean((nowcasts$`SDFM(BIC)` - nowcasts$Realisation)[month_three_ind]^2, na.rm = TRUE)),
             "        ", sprintf("%07.4f", mean((nowcasts$`DFM` - nowcasts$Realisation)[month_three_ind]^2, na.rm = TRUE)),
             "        ", sprintf("%07.4f", mean((nowcasts$`DFM EM` - nowcasts$Realisation)[month_three_ind]^2, na.rm = TRUE)),
             "        ", sprintf("%07.4f", mean((nowcasts$`SDFM EM` - nowcasts$Realisation)[month_three_ind]^2, na.rm = TRUE)), "\n",
             "__________________________________________________________________________________________________\n",
             "1-step-ahead MSFE    \n",
             " Month 1    ", sprintf("%07.4f", mean((one_step$`SDFM(CV)` - one_step$Realisation)[month_one_ind]^2, na.rm = TRUE)),
             "        ", sprintf("%07.4f", mean((one_step$`SDFM(BIC)` - one_step$Realisation)[month_one_ind]^2, na.rm = TRUE)),
             "        ", sprintf("%07.4f", mean((one_step$`DFM` - one_step$Realisation)[month_one_ind]^2, na.rm = TRUE)),
             "        ", sprintf("%07.4f", mean((one_step$`DFM EM` - one_step$Realisation)[month_one_ind]^2, na.rm = TRUE)),
             "        ", sprintf("%07.4f", mean((one_step$`SDFM EM` - one_step$Realisation)[month_one_ind]^2, na.rm = TRUE)), "\n",
             " Month 2    ", sprintf("%07.4f", mean((one_step$`SDFM(CV)` - one_step$Realisation)[month_two_ind]^2, na.rm = TRUE)),
             "        ", sprintf("%07.4f", mean((one_step$`SDFM(BIC)` - one_step$Realisation)[month_two_ind]^2, na.rm = TRUE)),
             "        ", sprintf("%07.4f", mean((one_step$`DFM` - one_step$Realisation)[month_two_ind]^2, na.rm = TRUE)),
             "        ", sprintf("%07.4f", mean((one_step$`DFM EM` - one_step$Realisation)[month_two_ind]^2, na.rm = TRUE)),
             "        ", sprintf("%07.4f", mean((one_step$`SDFM EM` - one_step$Realisation)[month_two_ind]^2, na.rm = TRUE)), "\n",
             " Month 3    ", sprintf("%07.4f", mean((one_step$`SDFM(CV)` - one_step$Realisation)[month_three_ind]^2, na.rm = TRUE)),
             "        ", sprintf("%07.4f", mean((one_step$`SDFM(BIC)` - one_step$Realisation)[month_three_ind]^2, na.rm = TRUE)),
             "        ", sprintf("%07.4f", mean((one_step$`DFM` - one_step$Realisation)[month_three_ind]^2, na.rm = TRUE)),
             "        ", sprintf("%07.4f", mean((one_step$`DFM EM` - one_step$Realisation)[month_three_ind]^2, na.rm = TRUE)),
             "        ", sprintf("%07.4f", mean((one_step$`SDFM EM` - one_step$Realisation)[month_three_ind]^2, na.rm = TRUE)), "\n",
             "__________________________________________________________________________________________________\n",
             "####################################################################################################\n"
  )
  )
  
  vintage_ind <- step_size + vintage_ind
  save.image("./safety_net_pkg_more_compet_09062026.RData")
}

# Nowcast evaluation #

# Create result plots #

# Create a long format data.table
results_combined_zoo <- merge(nowcasts$`SDFM(CV)`, nowcasts$`DFM`, nowcasts$`DFM EM`, nowcasts$`SDFM EM`, nowcasts$Realisation)
results_combined <- cbind(as.data.frame(as.yearqtr(time(results_combined_zoo)[month_three_ind])),
                          as.data.frame(coredata(results_combined_zoo[month_three_ind, ])))
colnames(results_combined) <- c("Dates", "SDFM", "DFM", "DFM EM", "SDFM EM", "Realisation")
results_long <- pivot_longer(results_combined[, 1:5], cols = c(SDFM, DFM, `DFM EM`, `SDFM EM`), names_to = "Model", values_to = "Value")

# Plot skelleton
nowcast_plot_base <- ggplot() +
  geom_hline(yintercept = 0) +
  geom_bar(data = results_combined, 
           aes(x = Dates, y = Realisation, fill = "Realisation"), 
           stat = "identity", alpha = 0.5, width = 0.1) +
  geom_line(data = results_long, 
            aes(x = Dates, y = Value, color = Model, group = Model, linetype = Model), 
            linewidth = 2) +
  labs(x = "Dates", y = expression(Delta ~ "log gdp")) +
  theme_minimal() 

# Full evaluation period plot
full_nowcast_plots <- nowcast_plot_base +
  scale_x_yearqtr(format = "%Y-Q%q", expand = c(0, 0)) +
  scale_fill_manual(values = c("Realisation" = "#85C0F9")) +
  scale_color_manual(values = c("SDFM" = "#000000", "DFM" = "#F5793A", "DFM EM" = "#0F9D58", "SDFM EM" = "#A95AA1")) +
  scale_linetype_manual(values = c("DFM" = "twodash", "SDFM EM" = "dotdash", "DFM EM" = "dotted", "SDFM" = "solid")) +
  theme(
    axis.text.x = element_text(angle = 30, vjust = 1, hjust = 1),
    text = element_text(size = 30),
    axis.line = element_line(color = "black"), 
    legend.title = element_blank()
  ) +
  theme()
full_nowcast_plots
ggsave("full_sample_nowcast_plots.pdf", full_nowcast_plots, width = 16, height = 9)

# Pre-Covid plot
nowcast_plot_pre_covid <- nowcast_plot_base +
  scale_x_yearqtr(format = "%Y-Q%q", expand = c(0, 0), 
                  limits = as.yearqtr(c("1999 Q3", "2019 Q4"))) +
  scale_y_continuous(limits = c(-2, 2)) +
  scale_fill_manual(values = c("Realisation" = "#85C0F9")) +
  scale_color_manual(values = c("SDFM" = "#000000", "DFM" = "#F5793A", "DFM EM" = "#0F9D58", "SDFM EM" = "#A95AA1")) +
  scale_linetype_manual(values = c("DFM" = "twodash", "SDFM EM" = "dotdash", "DFM EM" = "dotted", "SDFM" = "solid")) +
  theme(
    axis.text.x = element_text(angle = 30, vjust = 1, hjust = 1),
    text = element_text(size = 50),
    axis.line = element_line(color = "black"), 
    legend.title = element_blank()
  ) +
  theme(legend.position = "none")
nowcast_plot_pre_covid
ggsave("pre_covid_nowcasts_plot.pdf", nowcast_plot_pre_covid, width = 16, height = 9)

# Post-Covid plot
nowcast_plot_post_covid <- nowcast_plot_base +
  scale_x_yearqtr(format = "%Y-Q%q", expand = c(0, 0), 
                  limits = as.yearqtr(c("2020 Q4", "2025 Q4"))) +
  scale_y_continuous(limits = c(-2, 2)) +
  scale_fill_manual(values = c("Realisation" = "#85C0F9")) +
  scale_color_manual(values = c("SDFM" = "#000000", "DFM" = "#F5793A", "DFM EM" = "#0F9D58", "SDFM EM" = "#A95AA1")) +
  scale_linetype_manual(values = c("DFM" = "twodash", "SDFM EM" = "dotdash", "DFM EM" = "dotted", "SDFM" = "solid")) +
  theme(
    axis.text.x = element_text(angle = 30, vjust = 1, hjust = 1),
    text = element_text(size = 50),
    axis.line = element_line(color = "black"), 
    legend.title = element_blank()
  ) +
  theme(legend.position = "none")
nowcast_plot_post_covid
ggsave("post_covid_nowcast_plot.pdf", nowcast_plot_post_covid, width = 16, height = 9)

# Ribbon bands
nowcast_bands <- data.frame(
  Dates = as.yearqtr(index(nowcasts)),
  q005 = apply(marginal_density_nowcasts, 1, quantile, probs = 0.005, na.rm = TRUE),
  q025 = apply(marginal_density_nowcasts, 1, quantile, probs = 0.025, na.rm = TRUE),
  q25  = apply(marginal_density_nowcasts, 1, quantile, probs = 0.25, na.rm = TRUE),
  q75  = apply(marginal_density_nowcasts, 1, quantile, probs = 0.75, na.rm = TRUE),
  q975 = apply(marginal_density_nowcasts, 1, quantile, probs = 0.975, na.rm = TRUE),
  q995 = apply(marginal_density_nowcasts, 1, quantile, probs = 0.995, na.rm = TRUE)
)[month_three_ind, ]

# Plot skeleton
dist_plot_base <- ggplot() +
  geom_hline(yintercept = 0) +
  geom_line(data = results_combined, 
            aes(x = Dates, y = Realisation, fill = "Realisation"), 
            stat = "identity", linewidth = 1) +
  geom_ribbon(data = nowcast_bands,
              aes(x = Dates, ymin = q005, ymax = q995, fill = "99% interval"),
              alpha = 0.2) +
  geom_ribbon(data = nowcast_bands,
              aes(x = Dates, ymin = q025, ymax = q975, fill = "95% interval"),
              alpha = 0.35) +
  geom_ribbon(data = nowcast_bands,
              aes(x = Dates, ymin = q25, ymax = q75, fill = "75% interval"),
              alpha = 0.5) +
  labs(x = "Dates", y = expression(Delta ~ "log gdp")) +
  theme_minimal() 

# Full evaluation period density plot
full_dist_plots <- dist_plot_base +
  scale_x_yearqtr(format = "%Y-Q%q", expand = c(0, 0)) +
  scale_fill_manual(name = NULL, values = c("99% interval" = "#85C0F9", 
                                            "95% interval" = "#5B9BD5", "75% interval" = "#2F75B5")) +
  theme(
    axis.text.x = element_text(angle = 30, vjust = 1, hjust = 1),
    text = element_text(size = 30),
    axis.line = element_line(color = "black"), 
    legend.title = element_blank()
  ) +
  theme()
full_dist_plots
ggsave("full_dist_plots.pdf", full_dist_plots, width = 16, height = 9)

# Pre-Covid density plot
dist_plot_pre_covid <- dist_plot_base +
  scale_x_yearqtr(format = "%Y-Q%q", expand = c(0, 0), 
                  limits = as.yearqtr(c("1999 Q3", "2019 Q4"))) +
  scale_y_continuous(limits = c(-2, 2)) +
  scale_fill_manual(name = NULL, values = c("99% interval" = "#85C0F9", 
                                            "95% interval" = "#5B9BD5", "75% interval" = "#2F75B5")) +
  theme(
    axis.text.x = element_text(angle = 30, vjust = 1, hjust = 1),
    text = element_text(size = 50),
    axis.line = element_line(color = "black"), 
    legend.title = element_blank()
  ) +
  theme(legend.position = "none")
dist_plot_pre_covid
ggsave("pre_covid_dist_plot.pdf", dist_plot_pre_covid, width = 16, height = 9)

# Post-Covid density plot
dist_plot_post_covid <- dist_plot_base +
  scale_x_yearqtr(format = "%Y-Q%q", expand = c(0, 0), 
                  limits = as.yearqtr(c("2020 Q4", "2025 Q4"))) +
  scale_y_continuous(limits = c(-2, 2)) +
  scale_fill_manual(name = NULL, values = c("99% interval" = "#85C0F9", 
                                            "95% interval" = "#5B9BD5", "75% interval" = "#2F75B5")) +
  theme(
    axis.text.x = element_text(angle = 30, vjust = 1, hjust = 1),
    text = element_text(size = 50),
    axis.line = element_line(color = "black"), 
    legend.title = element_blank()
  ) +
  theme(legend.position = "none")
dist_plot_post_covid
ggsave("post_covid_dist_plot.pdf", dist_plot_post_covid, width = 16, height = 9)


# Create the result table #

# Period indicator
pre_corona_ind <- (round(as.numeric(time(nowcasts)), 3) < round(2020 + (3 - 1)/12, 3))
corona_ind <- (round(as.numeric(time(nowcasts)), 3) %in% round(seq(2020 + (3 - 1)/12, 2021 + (3 - 1)/12, by = 1/12), 3))
post_corona_ind <- (round(as.numeric(time(nowcasts)), 3) > round(2021 + (3 - 1)/12, 3))
dot_com_ind <- (round(as.numeric(time(nowcasts)), 3) %in% round(seq(2001 + (3 - 1)/12, 2001 + (12 - 1)/12, by = 1/12), 3))
financial_ind <- (round(as.numeric(time(nowcasts)), 3) %in% round(seq(2007 + (12 - 1)/12, 2009 + (6 - 1)/12, by = 1/12), 3))
no_crisis_ind <- !corona_ind & !post_corona_ind & !dot_com_ind & !financial_ind

# Create a relative and absolute MSNE results
createResultForMonth <- function(month_ind){
  msne_all_models <- (nowcasts$Realisation - nowcasts[, 1:5])[month_ind]^2
  results <- as.data.frame(matrix(NaN, 7, 5))
  colnames(results) <- colnames(one_step)[1:5]
  rownames(results) <- c("Full Sample w/o Covid", "Pre-Covid", "Covid", "Post-Covid", 
                         "Dotcom Crisis", "Financial Crisis", "No Crisis")
  results[1, ] <- colMeans(msne_all_models[!corona_ind[month_ind], ], na.rm = TRUE)
  results[2, ] <- colMeans(msne_all_models[pre_corona_ind[month_ind], ], na.rm = TRUE)
  results[3, ] <- colMeans(msne_all_models[corona_ind[month_ind], ], na.rm = TRUE)
  results[4, ] <- colMeans(msne_all_models[post_corona_ind[month_ind], ], na.rm = TRUE)
  results[5, ] <- colMeans(msne_all_models[dot_com_ind[month_ind], ], na.rm = TRUE)
  results[6, ] <- colMeans(msne_all_models[financial_ind[month_ind], ], na.rm = TRUE)
  results[7, ] <- colMeans(msne_all_models[no_crisis_ind[month_ind], ], na.rm = TRUE)
  return(results)
}

results_month_one <- createResultForMonth(month_one_ind)
rel_improv_res_month_one <- 100 * round(1 - results_month_one$`SDFM(CV)` / results_month_one, 4)
results_month_one
rel_improv_res_month_one

results_month_two <- createResultForMonth(month_two_ind)
rel_improv_res_month_two <- 100 * round(1 - results_month_two$`SDFM(CV)` / results_month_two, 4)
results_month_two
rel_improv_res_month_two

results_month_three <- createResultForMonth(month_three_ind)
rel_improv_res_month_three <- 100 * round(1 - results_month_three$`SDFM(CV)` / results_month_three, 4)
results_month_three
rel_improv_res_month_three

rel_imprv_res <- cbind(rel_improv_res_month_one[, c(-1, -2)],
                       rel_improv_res_month_two[, c(-1, -2)],
                       rel_improv_res_month_three[, c(-1, -2)])
rel_imprv_res

# Create absolute results tex table
full_results_all_month <- rbind(results_month_one,
                                results_month_two,
                                results_month_three)
full_results_all_month
print(xtable(full_results_all_month, digits = c(0, 3, 3, 3, 3, 3)), include.rownames = TRUE,
      floating = FALSE)

# Compute clark-west-type tests for the comparison of SDFM (CV) against all other models
clark_west_test <- matrix(NaN, 4, 9)
colnames(clark_west_test) <- colnames(rel_imprv_res)
rownames(clark_west_test) <- c("Regular S.E. $t$-Stat.", "Regular S.E. $p$-Val.",
                               "Newey-West S.E. $t$-Stat.", "Newey-West S.E. $p$-Val.")
clarkWestTest <- function(month_ind){
  results <- matrix(NaN, 4, 3)
  colnames(results) <- c("DFM", "DFM EM", "SDFM EM")
  for(bench in 1:3){
    # Compute the clark-west target
    sq_fcst_error_cv <- as.vector(nowcasts$Realisation - nowcasts$`SDFM(CV)`)[month_ind]^2
    sq_fcst_error_dense <- as.vector(nowcasts$Realisation - nowcasts[, colnames(results)[bench]])[month_ind]^2
    sq_fcst_diff <- as.vector(nowcasts[, colnames(results)[bench]] - nowcasts$`SDFM(CV)`)[month_ind]^2
    adj_loss_differential <- (sq_fcst_error_dense - sq_fcst_error_cv + sq_fcst_diff)[!corona_ind[month_ind]]
    
    # Clark-West approximate normal test for equal predicitive ability of nested models (Clark, T. E., & West, K. D. (2007). Approximately normal tests for equal predictive accuracy in nested models. Journal of econometrics, 138(1), 291-311.)
    clark_west_fit <- lm(adj_loss_differential ~ 1)
    model_summary <- summary(clark_west_fit)
    results[1, bench] <- model_summary$coefficients[3]
    results[2, bench] <- pt(results[1, bench], df = clark_west_fit$df.residual, lower.tail = FALSE)
    nw_var_cov <- NeweyWest(clark_west_fit)
    results[3, bench] <- model_summary$coefficients[1] / sqrt(nw_var_cov)
    results[4, bench] <- pt(results[3, bench], df = clark_west_fit$df.residual, lower.tail = FALSE)
  }
  
  return(round(results, 3))
}

clark_west_test[, 1:3] <- clarkWestTest(month_one_ind)
clark_west_test[, 4:6] <- clarkWestTest(month_two_ind)
clark_west_test[, 7:9] <- clarkWestTest(month_three_ind)
clark_west_test

# Compute Giacomini-White-type tests for the cond. comparison of SDFM (CV) against all other models
make_GW_test <- function(month_ind){
  
  # Create containers
  coef_estim_results <- as.data.frame(matrix(NaN, 7, 3))
  rownames(coef_estim_results) <- c("No Crisis", "Covid", "Post-Covid",
                                    "Dotcom Crisis", "Financial Crisis",
                                    "$\\Delta L_{t-1}$", "$\\Delta L_{t-2}$")
  colnames(coef_estim_results) <- c("DFM", "DFM EM", "SDFM EM")
  
  se_estim_results <- as.data.frame(matrix(NaN, 7, 3))
  rownames(se_estim_results) <- c("No Crisis", "Covid", "Post-Covid",
                                  "Dotcom Crisis", "Financial Crisis",
                                  "$\\Delta L_{t-1}$", "$\\Delta L_{t-2}$")
  colnames(se_estim_results) <- c("DFM", "DFM EM", "SDFM EM")
  
  
  p_val_results <- as.data.frame(matrix(NaN, 7, 3))
  rownames(p_val_results) <- c("No Crisis", "Covid", "Post-Covid",
                               "Dotcom Crisis", "Financial Crisis",
                               "$\\Delta L_{t-1}$", "$\\Delta L_{t-2}$")
  colnames(p_val_results) <- c("DFM", "DFM EM", "SDFM EM")
  
  # Create period dummies
  corona_dummy <- as.numeric(corona_ind[month_ind])
  post_corona_dummy <- as.numeric(post_corona_ind[month_ind])
  dot_com_dummy <- as.numeric(dot_com_ind[month_ind])
  financial_dummy <- as.numeric(financial_ind[month_ind])
  no_crisis_dummy <- as.numeric(no_crisis_ind[month_ind])
  
  for(bench in 1:dim(p_val_results)[2]) {
    sq_fcst_error_cv <- as.vector(nowcasts$Realisation - nowcasts$`SDFM(CV)`)[month_ind]^2
    sq_fcst_error_bench <- as.vector(nowcasts$Realisation - nowcasts[, colnames(p_val_results)[bench]])[month_ind]^2
    loss_differential <- ts(sq_fcst_error_bench - sq_fcst_error_cv)
    fit <- dynlm(loss_differential ~ -1 + no_crisis_dummy + corona_dummy
                 + post_corona_dummy + dot_com_dummy + financial_dummy
                 # + L(loss_differential, 1) + L(loss_differential, 2) # Uncomment for appendix robustness results
                 )
    nw_var_cov <- NeweyWest(fit)
    test_res <- coeftest(fit, df = fit$df.residual, vcov. = nw_var_cov)
    coef_estim_results[1:5, bench]  <- round(test_res[, 1], 3)
    se_estim_results[1:5, bench]  <- round(test_res[, 2], 3)
    p_val_results[1:5, bench] <- round(test_res[, 4], 3)
  }
  return(list("coef_estim_results" = coef_estim_results, 
              "se_estim_results" = se_estim_results, 
              "p_val_results" = p_val_results))
}

month_one_res <- make_GW_test(month_one_ind)
coef_estim_results_month_one <- month_one_res$coef_estim_results
coef_estim_results_month_one
se_estim_results_month_one <- month_one_res$se_estim_results
se_estim_results_month_one
p_val_results_month_one <- month_one_res$p_val_results
p_val_results_month_one

month_two_res <- make_GW_test(month_two_ind)
coef_estim_results_month_two <- month_two_res$coef_estim_results
coef_estim_results_month_two
se_estim_results_month_two <- month_two_res$se_estim_results
se_estim_results_month_two
p_val_results_month_two <- month_two_res$p_val_results
p_val_results_month_two

month_three_res <- make_GW_test(month_three_ind)
coef_estim_results_month_three <- month_three_res$coef_estim_results
coef_estim_results_month_three
se_estim_results_month_three <- month_three_res$se_estim_results
se_estim_results_month_three
p_val_results_month_three <- month_three_res$p_val_results
p_val_results_month_three

coef_estim_results <- cbind(coef_estim_results_month_one,
                            coef_estim_results_month_two,
                            coef_estim_results_month_three)
coef_estim_results

se_estim_results <- cbind(se_estim_results_month_one,
                          se_estim_results_month_two,
                          se_estim_results_month_three)
se_estim_results

p_val_results <- cbind(p_val_results_month_one,
                       p_val_results_month_two,
                       p_val_results_month_three)
p_val_results

# Recreate final table
coeffWithStars <- function(coef, p_val){
  stars <- ifelse(p_val < 0.01, "^{***}",
                  ifelse(p_val < 0.05, "^{**}",
                         ifelse(p_val < 0.1, "^{*}", "")))
  paste0("$", sprintf("%.3f", coef), stars, "$")
}

makeTexTable <- function(){
  cat(
    "% ", as.character(Sys.time()), "
  \\begin{table}
  \\begin{center}
  \\caption{Relative nowcasting results with statistical tests disaggregated by subperiods}\\label{tab::rel_results}
  \\resizebox{\\textwidth}{!}{\\begin{tabular}{l|ccc|ccc|ccc}
  \\hline
  \\hline 
   & \\multicolumn{3}{c|}{First Month of Quarter} & \\multicolumn{3}{c|}{Second Month of Quarter} & \\multicolumn{3}{c}{Third Month of Quarter} \\\\
  Time Period & DFM & DFM EM & SDFM EM & DFM & DFM EM & SDFM EM & DFM & DFM EM & SDFM EM \\\\
  \\hline
  \\multicolumn{10}{c}{Relative Improvement (\\%)}\\\\
  \\hline\n")
  for(row in 1:dim(rel_imprv_res)[1]){
    cat(paste(rownames(rel_imprv_res)[row], paste0(sprintf("%.2f", rel_imprv_res[row, ]), collapse = " & "), sep = " & "), "\\\\ \n")
  }
  cat("\\hline
  \\hline
  \\multicolumn{10}{c}{Clark-West-Type Test for Full Sample w/o Covid}\\\\
  \\hline\n")
  for(row in 1:dim(clark_west_test)[1]){
    cat(paste(rownames(clark_west_test)[row], paste0(sprintf("%.3f", clark_west_test[row, ]), collapse = " & "), sep = " & "), "\\\\ \n")
  }
  cat("\\hline
  \\hline
  \\multicolumn{10}{c}{Giacomini-White-Type Regression}\\\\
  \\hline\n")
  for(row in 1:dim(coef_estim_results)[1]){
    coeff_wit_stars <- c()
    for(col in 1:dim(coef_estim_results)[2]){
      coeff_wit_stars[col] <- coeffWithStars(coef_estim_results[row, col], p_val_results[row, col])
    }
    cat(paste(rownames(coef_estim_results)[row], paste0(coeff_wit_stars, collapse = " & "), sep = " & "), "\\\\\n",
        " &", paste0("(", paste0(sprintf("%.3f", se_estim_results[row, ]), collapse = ") & ("), ")"), " \\\\\n")
  }
  cat("\\hline
  \\hline
  \\end{tabular}}
  \\end{center}
  \\begin{tablenotes}
  \\scriptsize\\item $^*$: $p < 0.10$, $^{**}$: $p < 0.05$, $^{***}$: $p < 0.01$
  \\item The crisis periods are constructed according to \\citet{hamilton2026usrecessions} as: 2020Q1 to 2021Q1 for ``Covid'', 2001Q1 to 2001Q4 for ``Dotcom Crisis'', 2007Q4 to 2009Q2 for ``Financial Crisis''. ``Full Sample w/o Covid'' refers to the entire out-of-sample period from 1999Q3 to 2025Q3 excluding the covid period from 2020Q1 to 2021Q1. ``Pre-Covid'' refers to the period of 1999Q3 to 2019Q4. ``Post-Covid'' refers to the period from 2021Q2 to 2025Q3. All periods outside of the periods of ``Covid'', ``Post-Covid'', ``Dotcom Crisis'', and ``Financial Crisis'' are referred to as ``No Crisis''. The MSNE Reduction is computed as in \\autoref{sec::sim}.
  \\item ``DFM'' denotes the two-step estimator for DFMs by \\citet{Giannone2008Nowcasting}. ``DFM EM'' denotes the expectation-maximisation estimator for DFMs by \\citet{banbura2014maximum}. ``SDFM EM'' denotes the expectation-maximisation estimator for SDFMs by \\citet{mosley2023sparse} with hyper-parameters validated via the BIC.
  \\item $\\Delta L_{t-1}$ and $\\Delta L_{t-2}$ denote the first and second lag of the loss differential between the two models, i.e., $\\Delta L_t := L_{t,\\text{SDFM}} - L_{t,i}$, where, for all $t$, $L_{t,i} = (\\widehat{x}_{t,i} - x_{t,i})^2$ for nowcast $\\widehat{x}_{t,i}$ retrieved via model $i\\in\\{\\text{DFM}, \\text{DFM EM}, \\text{SDFM EM}\\}$ of observation $x_t$. The individual coefficient significant tests where conducted using heteroscedasticity- and autocorrelation-robust variance covariance matrix estimators. The same variance-covariance estimators have been employed in the ``Robust S.E.'' case of the Clark-West-type test.
  \\end{tablenotes} 
  \\end{table}")
}
makeTexTable()