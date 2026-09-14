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
load("./current_results_final_29062026.RData")
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
date_ind <- as.numeric(as.yearqtr(as.Date(Realisation[, 2])))

all_names <- unlist(read.table("./clean_data/fredgdp_real_time/vintage_2025-09-01/names.csv", sep = ",", header = FALSE))

# Individual variables
nowcasts <- as.zoo(ts(matrix(NA, no_of_predictions, length(all_names) + 3 - 1),
                      start = c(year(start_date), month(start_date)),
                      end = c(year(end_date), month(end_date)),
                      frequency = 12))
colnames(nowcasts) <- c("Full Model", all_names[-1], "Realisation", "Target Date")
nowcasts$`Target Date` <- as.yearqtr(time(nowcasts))
nowcast_ind <- which(nowcasts$`Target Date` %in%  as.numeric(as.yearqtr(as.Date(Realisation$V2))))
reali_nowcast_ind <- which(as.numeric(as.yearqtr(as.Date(Realisation$V2))) %in% nowcasts$`Target Date`)
nowcasts$Realisation[nowcast_ind] <- Realisation[reali_nowcast_ind, 1]
month_three_ind <- which(month(time(nowcasts)) %% 3 == 0)
month_two_ind <- which(month(time(nowcasts)) %% 3 == 2)
month_one_ind <- which(month(time(nowcasts)) %% 3 == 1)
nowcasts$Realisation[c(month_one_ind, month_two_ind)] <- NA
nowcasts$Realisation <- na.locf(nowcasts$Realisation, fromLast = TRUE, na.rm = FALSE)

one_step <- as.zoo(ts(matrix(NA, no_of_predictions, length(all_names) + 3 - 1),
                      start = c(year(start_date), month(start_date)),
                      end = c(year(end_date), month(end_date)),
                      frequency = 12))
colnames(one_step) <- c("Full Model", all_names[-1], "Realisation", "Target Date")
one_step$`Target Date` <- nowcasts$`Target Date` + 0.25
osh_ind <- which(one_step$`Target Date` %in%  as.numeric(as.yearqtr(as.Date(Realisation$V2))))
reali_osh_ind <- which(as.numeric(as.yearqtr(as.Date(Realisation$V2))) %in% one_step$`Target Date`)
one_step$Realisation[osh_ind] <- Realisation[reali_osh_ind, 1]
one_step$Realisation[c(month_one_ind, month_two_ind)] <- NA
one_step$Realisation <- na.locf(one_step$Realisation, fromLast = TRUE, na.rm = FALSE)

# Groups
output_income <- c("RPI", "W875RX1", "INDPRO", "IPFPNSS", "IPFINAL", "IPCONGD", "IPDCONGD",
                   "IPNCONGD", "IPBUSEQ", "IPMAT", "IPDMAT", "IPNMAT", "IPMANSICS",
                   "IPB51222S", "IPFUELS", "NAPMPI", "CUMFNS")

labour_market <- c("HWI", "HWIURATIO", "CLF16OV", "CE16OV", "UNRATE", "UEMPMEAN", 
                   "UEMPLT5", "UEMP5TO14", "UEMP15OV", "UEMP15T26", "UEMP27OV",
                   "CLAIMSx", "PAYEMS", "USGOOD", "CES1021000001", "USCONS",
                   "MANEMP", "DMANEMP", "NDMANEMP", "SRVPRD", "USTPU", "USWTRADE",
                   "USTRADE", "USFIRE", "USGOVT", "CES0600000007", "AWOTMAN", "AWHMAN",
                   "NAPMEI", "CES0600000008", "CES2000000008", "CES3000000008")

consumption_and_orders <- c("HOUST", "HOUSTNE", "HOUSTMW", "HOUSTS", "HOUSTW", "PERMIT",
                            "PERMITNE", "PERMITMW", "PERMITS", "PERMITW")

orders_and_invent <- c("DPCERA3M086SBEA", "CMRMTSPLx", "RETAILx", "NAPM", "NAPMNOI",
                       "NAPMSDI", "NAPMII", "ACOGNO", "AMDMNOx", "ANDENOx", "AMDMUOx",
                       "BUSINVx", "ISRATIOx", "UMCSENTx")

money_and_credit <- c("M1SL", "M2SL", "M2REAL", "BOGMBASE", "TOTRESNS", "NONBORRES",
                      "BUSLOANS", "BUSLOANS", "REALLN", "NONREVSL", "CONSPI", "MZMSL",
                      "DTCOLNVHFNM", "DTCTHFNM", "INVEST", "AMBSL",
                      "MZMSL")

int_and_ex <- c("FEDFUNDS", "CP3Mx", "TB3MS", "TB6MS", "GS1", "GS5", "GS10", "AAA", 
                "BAA", "COMPAPFFx", "TB3SMFFM", "TB6SMFFM", "T1YFFM", "T5YFFM",
                "T10YFFM", "AAAFFM", "BAAFFM", "TWEXAFEGSMTHx", "EXSZUSx", "EXJPUSx",
                "EXUSUKx", "EXCAUSx", "TWEXMMTH")

prices <- c("WPSFD49207", "WPSFD49502", "WPSID61", "WPSID62", "OILPRICEx", "PPICMM", "NAPMPRI",
            "CPIAUCSL", "CPIAPPSL", "CPITRNSL", "CPIMEDSL", "CUSR0000SAC", "CUSR0000SAD",
            "CUSR0000SAS", "CPIULFSL", "CUSR0000SA0L2", "CUSR0000SA0L5", "PCEPI",
            "DDURRG3M086SBEA", "DNDGRG3M086SBEA", "DSERRG3M086SBEA",  "PPIFGS", "PPIFCG",
            "PPIITM", "PPICRM", "CUUR0000SAD", "CUUR0000SA0L2")


stock_market <- c("S.P.500", "S.P.indust", "S.P.div.yield", "S.P.PE.ratio", "VIXCLSx", "S.P..indust",
                  "VXOCLSx")

# By groups
group_nowcasts <- as.zoo(ts(matrix(NA, no_of_predictions, 11),
                            start = c(year(start_date), month(start_date)),
                            end = c(year(end_date), month(end_date)),
                            frequency = 12))
colnames(group_nowcasts) <- c("Full Model", "Output & Income", "Labour Market",
                              "Consump. & Orders", "Orders & Inventories", "Money & Credit",
                              "Interest Rates & FX", "Prices", "Stock Market", "Realisation", "Target Date")
group_nowcasts$`Target Date` <- as.yearqtr(time(group_nowcasts))
group_nowcasts$Realisation[nowcast_ind] <- Realisation[reali_nowcast_ind, 1]
group_nowcasts$Realisation[c(month_one_ind, month_two_ind)] <- NA
group_nowcasts$Realisation <- na.locf(group_nowcasts$Realisation, fromLast = TRUE, na.rm = FALSE)

group_one_step <- as.zoo(ts(matrix(NA, no_of_predictions, 11),
                            start = c(year(start_date), month(start_date)),
                            end = c(year(end_date), month(end_date)),
                            frequency = 12))
colnames(group_one_step) <- c("Full Model", "Output & Income", "Labour Market",
                              "Consump. & Orders", "Orders & Inventories", "Money & Credit",
                              "Interest Rates & FX", "Prices", "Stock Market", "Realisation", "Target Date")
group_one_step$`Target Date` <- nowcasts$`Target Date` + 0.25
group_one_step$Realisation[osh_ind] <- Realisation[reali_osh_ind, 1]
group_one_step$Realisation[c(month_one_ind, month_two_ind)] <- NA
group_one_step$Realisation <- na.locf(group_one_step$Realisation, fromLast = TRUE, na.rm = FALSE)

# start prediction loop over all vintages #
vintage_ind <- 3  # Start at the end of the first quarter
vintage_name <- vintage_file_list[[vintage_ind]]
step_size <- 1
no_of_factors_old <- 0
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
  
  # SDFM nowcasting
  full_nowcast_cv <- TwoStepSDFM::nowcast(data = data_zoo, variables_of_interest = variables_of_interest,
                                          max_fcast_horizon = 2, delay = delay,
                                          selected = cv_results_test$CV$`Min. CV`[3:(3 + no_of_factors_test - 1)],
                                          frequency = frequency, no_of_factors = no_of_factors_test,
                                          max_ar_lag_order = 1, max_predictor_lag_order = 1,
                                          ridge_penalty = cv_results_test$CV$`Min. CV`[2]
  )
  predictions_cv_full <- na.omit(full_nowcast_cv$Forecasts$`Fcast GDP`)
  nowcasts[vintage_ind, 1] <- centering[variables_of_interest] + predictions_cv_full[1] * scaling[variables_of_interest]
  one_step[vintage_ind, 1] <- centering[variables_of_interest] + predictions_cv_full[2] * scaling[variables_of_interest]
  group_nowcasts[vintage_ind, 1] <- centering[variables_of_interest] + predictions_cv_full[1] * scaling[variables_of_interest]
  group_one_step[vintage_ind, 1] <- centering[variables_of_interest] + predictions_cv_full[2] * scaling[variables_of_interest]
  
  full_fit <- twoStepSDFM(data = data_zoo[, which(frequency == 12)], delay = delay[which(frequency == 12)],
                          selected = cv_results_test$CV$`Min. CV`[3:(3 + no_of_factors_test - 1)],
                          no_of_factors = no_of_factors_test, ridge_penalty = cv_results_test$CV$`Min. CV`[2])
  full_loading <- full_fit$loading_matrix_estimate
  
  # By group #
  output_income_ind <- which(colnames(data_zoo) %in% output_income)
  labour_market_ind <- which(colnames(data_zoo) %in% labour_market)
  consumption_and_orders_ind <- which(colnames(data_zoo) %in% consumption_and_orders)
  orders_and_invent_ind <- which(colnames(data_zoo) %in% orders_and_invent)
  money_and_credit_ind <- which(colnames(data_zoo) %in% money_and_credit)
  int_and_ex_ind <- which(colnames(data_zoo) %in% int_and_ex)
  prices_ind <- which(colnames(data_zoo) %in% prices)
  stock_market_ind <- which(colnames(data_zoo) %in% stock_market)
  all_names_sofar <- c(output_income, labour_market, consumption_and_orders,
                       orders_and_invent, money_and_credit, int_and_ex,
                       prices, stock_market)
  
  all_ind <- list(output_income_ind, labour_market_ind, consumption_and_orders_ind,
                  orders_and_invent_ind, money_and_credit_ind, int_and_ex_ind,
                  prices_ind, stock_market_ind)
  
  for(group in 1:length(all_ind)){
    curr_data <- data_zoo[, c(-1, -all_ind[[group]])]
    var_ind <- which(colnames(data_zoo[, frequency == 12]) %in% colnames(data_zoo)[all_ind[[group]]])
    if(all(full_loading[var_ind, ] == 0)){
      group_nowcasts[vintage_ind, group + 1] <- nowcasts[vintage_ind, 1]
      group_one_step[vintage_ind, group + 1] <- one_step[vintage_ind, 1]
      next
    }
    curr_loading_matrix <- full_loading[-var_ind, ]
    curr_factor_mat <- curr_data %*% curr_loading_matrix %*% solve(t(curr_loading_matrix) %*% curr_loading_matrix)
    curr_factor_fit <- na.omit(zoo(curr_factor_mat, order.by = index(data_zoo)))
    lag_order <- full_fit$factor_var_lag_order
    stacked_factors <- embed(curr_factor_fit, lag_order + 1)
    lags <- stacked_factors[, -(1:ncol(curr_factor_fit)), drop = FALSE]
    current <- stacked_factors[, 1:ncol(curr_factor_fit), drop = FALSE]
    var_coeff <- t(solve(t(lags) %*% lags) %*% t(lags) %*% current)
    curr_trans_var_coeff <- rbind(var_coeff, cbind(diag(no_of_factors_test * (lag_order - 1)), 
                                                   matrix(0, no_of_factors_test * (lag_order - 1), no_of_factors_test))
    )
    factor_errors <- current - lags %*% t(var_coeff)
    curr_trans_error_var_cov <- (1 / (dim(factor_errors)[1] - 1)) * t(factor_errors) %*% factor_errors
    meas_errors <- na.omit(curr_data) - curr_factor_fit %*% t(curr_loading_matrix)
    curr_meas_error_var_cov <- (1 / (dim(factor_errors)[1] - 1)) * t(meas_errors) %*% meas_errors
    curr_companion_loading_matrix <- cbind(curr_loading_matrix, matrix(0, dim(curr_data)[2],no_of_factors_test * lag_order))
    
    # Filter and smooth factors
    kalman_fit <- kalmanFilterSmoother(curr_data, delay[-c(1, all_ind[[group]])], 
                                       no_of_factors = no_of_factors_test,
                                       loading_matrix = curr_loading_matrix, 
                                       meas_error_var_cov = curr_meas_error_var_cov,
                                       trans_error_var_cov = curr_trans_error_var_cov,
                                       trans_var_coeff = var_coeff,
                                       factor_lag_order = lag_order)
    curr_factors_kfs <- kalman_fit$smoothed_factors
    
    curr_pred <- MIDASForecaster(target_variables = t(coredata(data_zoo[, 1, drop = FALSE])),
                                 monthly_predictors = t(coredata(curr_factors_kfs)),
                                 target_variable_delay = delay[1],
                                 predictors_delay = rep(0, no_of_factors_test),
                                 lag_estim_criterion = "BIC",
                                 max_fcast_horizon = 2,
                                 max_ar_lag_order = 1,
                                 max_predictor_lag_order = 1)
    
    group_nowcasts[vintage_ind, group + 1] <- centering[variables_of_interest] + curr_pred$`Avg. Point Forecast`[1] * scaling[variables_of_interest]
    group_one_step[vintage_ind, group + 1] <- centering[variables_of_interest] + curr_pred$`Avg. Point Forecast`[2] * scaling[variables_of_interest]
  }
  
  print(group_nowcasts[vintage_ind, ])
  vintage_ind <- vintage_ind + step_size
}

corona_ind <- (round(as.numeric(time(group_nowcasts)), 3) %in% round(seq(2020 + (3 - 1)/12, 2021 + (3 - 1)/12, by = 1/12), 3))
month_one_res <- ((group_nowcasts[, 1:9] - group_nowcasts$Realisation)^2)[month_one_ind]
month_one_msne <- colMeans(month_one_res[!corona_ind[month_one_ind], ], na.rm = TRUE)
month_one_msne

month_two_res <- ((group_nowcasts[, 1:9] - group_nowcasts$Realisation)^2)[month_two_ind]
month_two_msne <- colMeans(month_two_res[!corona_ind[month_two_ind], ], na.rm = TRUE)
month_two_msne

month_three_res <- ((group_nowcasts[, 1:9] - group_nowcasts$Realisation)^2)[month_three_ind]
month_three_msne <- colMeans(month_three_res[!corona_ind[month_three_ind], ], na.rm = TRUE)
month_three_msne

var_log_gdp <- var(group_nowcasts$Realisation[month_three_ind][!corona_ind[month_three_ind], ], na.rm = TRUE)

deco_msne <- rbind(month_one_msne, month_two_msne, month_three_msne)
rownames(deco_msne) <- c("Month 1", "Month 2", "Month 3")

deco_msne_prepared <- deco_msne %>%
  as.data.frame() %>%
  mutate(Month = rownames(.)) %>%
  pivot_longer(cols = -Month, names_to = "Group", values_to = "MSNE")
deco_msne_prepared$Group <- factor(deco_msne_prepared$Group, levels = colnames(deco_msne))

cumm_importance_plot <- ggplot(deco_msne_prepared, aes(x = Group, y = MSNE, fill = Group)) +
  geom_col(width = 0.8) +
  geom_hline(aes(yintercept = var_log_gdp, linetype = "Var"),
             colour = "black") +
  facet_wrap(~ Month, nrow = 1, scales = "free_x", axes = "all_y") +
  scale_fill_manual(values = c("Full Model" = "#000000", "Output & Income" = "#0072B2",
                               "Labour Market" = "#E69F00", "Consump. & Orders" = "#009E73",
                               "Orders & Inventories" = "#D55E00", "Money & Credit" = "#CC79A7",
                               "Interest Rates and FX" = "#56B4E9", "Prices" = "#F0E442",
                               "Stock Market" = "#999999")) +
  scale_linetype_manual(values = c("Var" = "dashed"), labels = expression(Var(Delta*log~gdp))) +
  guides(fill = "none", linetype = guide_legend(title = NULL)) +
  labs(x = NULL, y = "MSNE") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 75, vjust = 1, hjust = 1),
    text = element_text(size = 30),
    axis.line = element_line(color = "black"),
    legend.position = "top",
    panel.spacing = unit(2, "lines")
  )
cumm_importance_plot
ggsave("cumm_importance_plot.pdf", cumm_importance_plot, width = 16, height = 9)

# Plot msne deco for post covid
post_corona_ind <- (round(as.numeric(time(nowcasts)), 3) > round(2021 + (3 - 1)/12, 3))
post_covid_group_deco <- fortify(na.omit(group_nowcasts[month_three_ind, -11][post_corona_ind[month_three_ind], ]))
colnames(post_covid_group_deco)[1] <- "Date"
post_covid_group_deco$Date <- as.yearqtr(post_covid_group_deco$Date)
full_model <- post_covid_group_deco$`Full Model`
post_covid_group_deco_long <- post_covid_group_deco[, -c(2, 11)] %>%
  pivot_longer(cols = -Date, names_to = "Group", values_to = "Value") %>%
  mutate(Base = full_model[match(Date, post_covid_group_deco$Date)],
         Delta = Value - Base)  %>%
  group_by(Date) %>%
  arrange(desc(abs(Delta)), .by_group = TRUE) %>%
  ungroup()
post_covid_group_deco_long$Group <- factor(post_covid_group_deco_long$Group, levels = colnames(deco_msne)[-1])

group_impoirtance_plot_post_covid <- ggplot() +
  scale_x_yearqtr(format = "%Y-Q%q") +
  geom_segment(data = post_covid_group_deco_long, aes(x = Date, xend = Date, y = Base,
                                                      yend = Base + Delta, color = Group),
               alpha = 1, linewidth = 8, lineend = "butt", group = seq_len(nrow(post_covid_group_deco_long))) +
  geom_line(data = post_covid_group_deco, aes(x = Date, y = Realisation, group = 1), color = "black", linewidth = 2) +
  geom_point(data = post_covid_group_deco, aes(x = Date, y = `Full Model`), color = "#000000", size = 5, shape = 18) +
  scale_color_manual(values = c("Output & Income" = "#0072B2",
                                "Labour Market" = "#E69F00", "Consump. & Orders" = "#009E73",
                                "Orders & Inventories" = "#D55E00", "Money & Credit" = "#CC79A7",
                                "Interest Rates and FX" = "#56B4E9", "Prices" = "#F0E442",
                                "Stock Market" = "#999999")) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 30, vjust = 1, hjust = 1),
    text = element_text(size = 30),
    axis.line = element_line(color = "black"),
    legend.title = element_blank(),
    legend.position = "top",
    legend.box = "horizontal"
  ) +
  labs(x = "Date", y = expression(Delta*log~gdp))
group_impoirtance_plot_post_covid
ggsave("group_impoirtance_plot_post_covid.pdf", group_impoirtance_plot_post_covid, width = 16, height = 9)

# Plot msne deco for pre covid
pre_corona_ind <- index(nowcasts)[(round(as.numeric(time(nowcasts)), 3) < round(2020 + (3 - 1)/12, 3))]
pre_covid_group_deco <- fortify(na.omit(group_nowcasts[month_three_ind, -11][pre_corona_ind[month_three_ind], ]))
colnames(pre_covid_group_deco)[1] <- "Date"
pre_covid_group_deco$Date <- as.yearqtr(pre_covid_group_deco$Date)
full_model <- pre_covid_group_deco$`Full Model`
pre_covid_group_deco_long <- pre_covid_group_deco[, -c(2, 11)] %>%
  pivot_longer(cols = -Date, names_to = "Group", values_to = "Value") %>%
  mutate(Base = full_model[match(Date, pre_covid_group_deco$Date)],
         Delta = Value - Base)  %>%
  group_by(Date) %>%
  arrange(desc(abs(Delta)), .by_group = TRUE) %>%
  ungroup()
pre_covid_group_deco_long$Group <- factor(pre_covid_group_deco_long$Group, levels = colnames(deco_msne)[-1])

group_impoirtance_plot_pre_covid <- ggplot() +
  scale_x_yearqtr(format = "%Y-Q%q") +
  geom_segment(data = pre_covid_group_deco_long, aes(x = Date, xend = Date, y = Base,
                                                     yend = Base + Delta, color = Group),
               alpha = 1, linewidth = 5, lineend = "butt", group = seq_len(nrow(pre_covid_group_deco_long))) +
  geom_line(data = pre_covid_group_deco, aes(x = Date, y = Realisation, group = 1), color = "black", linewidth = 2) +
  geom_point(data = pre_covid_group_deco, aes(x = Date, y = `Full Model`), color = "#000000", size = 5, shape = 18) +
  scale_color_manual(values = c("Output & Income" = "#0072B2",
                                "Labour Market" = "#E69F00", "Consump. & Orders" = "#009E73",
                                "Orders & Inventories" = "#D55E00", "Money & Credit" = "#CC79A7",
                                "Interest Rates and FX" = "#56B4E9", "Prices" = "#F0E442",
                                "Stock Market" = "#999999")) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 30, vjust = 1, hjust = 1),
    text = element_text(size = 30),
    axis.line = element_line(color = "black"),
    legend.title = element_blank(),
    legend.position = "top",
    legend.box = "horizontal"
  ) +
  labs(x = "Date", y = expression(Delta*log~gdp))
group_impoirtance_plot_pre_covid
ggsave("group_impoirtance_plot_pre_covid.pdf", group_impoirtance_plot_pre_covid, width = 16, height = 9)


library(knitr)
library(kableExtra)
group_nowcasts_moht_three <- fortify(na.omit(group_nowcasts[month_three_ind, 1:10]))
colnames(group_nowcasts_moht_three)[1] <- "Date"
kable(group_nowcasts_moht_three, format = "latex", booktabs = TRUE, digits = 3, 
      caption = "Group omission impact on the full model nowcast compared to the realisations
from 1999Q3 to 2025Q2") %>%
  kable_styling(
    latex_options = c("landscape", "scale_down", "hold_position")
  )

group_nowcasts_moht_three %>%
  kable(
    format = "latex",
    booktabs = TRUE,
    longtable = TRUE,
    align = c("l", rep("c", ncol(group_nowcasts_moht_three)-1)),
    digits = 3,
    caption = "Group omission impact on the full model nowcast compared to the realisations from 1999Q3 to 2025Q2"
  ) %>%
  kable_styling(
    latex_options = c("repeat_header"),
    position = "center"
  )
