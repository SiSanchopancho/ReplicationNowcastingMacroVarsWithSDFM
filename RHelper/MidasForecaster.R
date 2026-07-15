# SPDX-License-Identifier: GPL-3.0-or-later
#
#  Copyright (C) 2024 Domenic Franjic
#
#  This file is part of TwoStepSDFM.
#
#  TwoStepSDFM is free software: you can redistribute
#  it and/or modify it under the terms of the GNU General Public License as
#  published by the Free Software Foundation, either version 3 of the License,
#  or (at your option) any later version.
#
#  TwoStepSDFM is distributed in the hope that it
#  will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
#  of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with TwoStepSDFM. If not, see <https://www.gnu.org/licenses/>.

MIDASForecaster <- function(target_variables,
                            monthly_predictors,
                            target_variable_delay,
                            predictors_delay,
                            lag_estim_criterion,
                            max_fcast_horizon,
                            max_ar_lag_order,
                            max_predictor_lag_order
)
{
  target_variables <- cbind(matrix(NaN, dim(target_variables)[1], 2), target_variables)
  monthly_predictors <- cbind(matrix(NaN, dim(monthly_predictors)[1], 2), monthly_predictors)
  no_of_target_vars <- dim(target_variables)[1]
  no_of_predictors <- dim(monthly_predictors)[1]
  no_of_vars  <- no_of_predictors + no_of_target_vars
  no_of_qrtly_obs <- (dim(monthly_predictors)[2] - 2) / 3
  min_fcast_horizons <- ifelse(target_variable_delay == 0, 
                               1,
                               -floor(target_variable_delay / 3) + 1)
  return_object <- list()
  all_qtrly_data_delay <- c(floor(c(target_variable_delay, predictors_delay) / 3))
  
  # Start quarterfication loop according to Mariano and Murasawa #
  
  # Note: It is implicitly assumed that the stationary monthly data starts at the 
  #   second month of the first quarter. Further, it is assumed that the data set
  #   ends with an observation in the last month of the quarter.
  
  all_qrtly_data <- matrix(NaN, no_of_vars, no_of_qrtly_obs)
  for(t in seq(5, dim(monthly_predictors)[2], 3)){
    all_qrtly_data[1:no_of_target_vars, (t - 2)/3] <- target_variables[, t]
    all_qrtly_data[(no_of_target_vars + 1):(no_of_vars), (t - 2)/3] <- rowSums(
      cbind(1/3 * monthly_predictors[, t, drop = FALSE],
            2/3 * monthly_predictors[, t - 1, drop = FALSE],
            1 * monthly_predictors[, t - 2, drop = FALSE],
            2/3 * monthly_predictors[, t - 3, drop = FALSE],
            1/3 * monthly_predictors[, t - 4, drop = FALSE]),
      na.rm = TRUE)
  }
  
  # End quarterfication loop according to Mariano and Murasawa #
  
  # Start ARDL estimation loop over the target variables #
  
  fcasts <- matrix(NaN, no_of_target_vars, max_fcast_horizon - min(min_fcast_horizons) + 1)
  for(current_target in 1){
    
    # Start ARDL estimation loop over the predictor #
    
    current_fcasts <- matrix(NaN, no_of_vars, max_fcast_horizon - min_fcast_horizons[current_target] + 1)
    for(current_predictor in 1:(no_of_vars)){
      if(current_target == current_predictor){
        next # Skip using the current target_variable as a single predictor as its always included as predictor
      }      
      if(all_qtrly_data_delay[current_target] < all_qtrly_data_delay[current_predictor]){
        next # Skip a predictor if it is dalyed further back compared to the target variable (we do not expect forecasting gains from using variables that are further behind then the target)
      }
      
      rel_fcast_horizons <- min_fcast_horizons[current_target]:max_fcast_horizon + all_qtrly_data_delay[current_predictor]
      for(h in rel_fcast_horizons){
        
        # Fit the model for the specific forecasting horizon
        horizon_adjustment <- which(rel_fcast_horizons == h)
        
        horizon_specific_target <- matrix(
          all_qrtly_data[current_target, 
                         (horizon_adjustment + 1):(no_of_qrtly_obs - all_qtrly_data_delay[current_target])],
          ncol = 1)
        
        horizon_specific_ar_lag <- matrix(
          all_qrtly_data[current_target, 
                         1:(no_of_qrtly_obs - all_qtrly_data_delay[current_target] - horizon_adjustment)],
          ncol = 1)
        
        horizon_specific_predictor <- matrix(
          all_qrtly_data[current_predictor, 
                         (horizon_adjustment + 1 - h):(no_of_qrtly_obs - all_qtrly_data_delay[current_target] - h)],
          ncol = 1)
        
        if(max_ar_lag_order != 0){
          ardl_fit <- runARDL(horizon_specific_target,
                              horizon_specific_ar_lag,
                              horizon_specific_predictor,
                              max(max_ar_lag_order - max(h, 0), 1), 
                              max(max_predictor_lag_order - max(h, 0), 1),
                              lag_estim_criterion)
          
          # Forecast
          forecast_predictors <- matrix(1, sum(ardl_fit$optimL_lag_order) + 3, 1) # Add three for the intercept and the "contemporaenous" observations
          
          forecast_predictors[2:(ardl_fit$optimL_lag_order[1] + 2), ] <- 
            head(all_qrtly_data[current_target, 
                                (no_of_qrtly_obs - all_qtrly_data_delay[current_target]):1],
                 ardl_fit$optimL_lag_order[1] + 1)
          
          forecast_predictors[(ardl_fit$optimL_lag_order[1] + 3):(ardl_fit$optimL_lag_order[1] + ardl_fit$optimL_lag_order[2] + 3), ] <- 
            head(all_qrtly_data[current_predictor, 
                                (no_of_qrtly_obs - all_qtrly_data_delay[current_predictor]):1],
                 ardl_fit$optimL_lag_order[2] + 1)
          
          current_fcasts[current_predictor, which(rel_fcast_horizons == h)] <-
            matrix(ardl_fit$coefficients, nrow = 1) %*% forecast_predictors
        }else{
          ardl_fit <- runDL(horizon_specific_target,
                            horizon_specific_predictor,
                            max(max_predictor_lag_order - max(h, 0), 1),
                            lag_estim_criterion)
          
          # Forecast
          forecast_predictors <- matrix(1, ardl_fit$optimL_lag_order + 2, 1) # Add two for the intercept and the "contemporaenous" observations
          
          forecast_predictors[2:(ardl_fit$optimL_lag_order[1] + 2), ] <- 
            head(all_qrtly_data[current_predictor, 
                                (no_of_qrtly_obs - all_qtrly_data_delay[current_predictor]):1],
                 ardl_fit$optimL_lag_order + 1)
          
          current_fcasts[current_predictor, which(rel_fcast_horizons == h)] <-
            matrix(ardl_fit$coefficients, nrow = 1) %*% forecast_predictors
        }
        
      }
      
      # End loop over the forecasting horizons #
      
    }
    
    # Store the final point forecast using simple forecast averaging for each target
    rownames(current_fcasts) <- c(rownames(target_variables), rownames(monthly_predictors))
    return_object[[current_target]] <- current_fcasts
    names(return_object)[current_target] <- paste0("Single Predictor Forecasts ", rownames(target_variables)[current_target], collapse = "")
    fcasts[current_target, (max_fcast_horizon - min(min_fcast_horizons) - length(rel_fcast_horizons) + 2):(max_fcast_horizon - min(min_fcast_horizons) + 1)] <-
      colMeans(current_fcasts, na.rm = TRUE)
    
    # Start ARDL estimation loop over the predictor #
    
  }
  
  # End ARDL estimation loop over the target variables #
  
  rownames(fcasts) <- rownames(target_variables)
  return_object[[no_of_target_vars + 1]] <- fcasts
  names(return_object)[no_of_target_vars + 1] <- "Avg. Point Forecast"
  
  return(return_object)
  
}


makeMarianoMurasawaData <- function(target_variables,
                                    monthly_predictors,
                                    target_variable_delay,
                                    predictors_delay
)
{
  monthly_index <- index(target_variables)
  target_variables <- cbind(matrix(NaN, dim(target_variables)[2], 2), t(coredata(target_variables)))
  monthly_predictors <- cbind(matrix(NaN, dim(monthly_predictors)[2], 2), t(coredata(monthly_predictors)))
  no_of_target_vars <- dim(target_variables)[1]
  no_of_predictors <- dim(monthly_predictors)[1]
  no_of_vars  <- no_of_predictors + no_of_target_vars
  no_of_qrtly_obs <- (dim(monthly_predictors)[2] - 2) / 3
  min_fcast_horizons <- ifelse(target_variable_delay == 0, 
                               1,
                               -floor(target_variable_delay / 3) + 1)
  return_object <- list()
  all_qtrly_data_delay <- c(floor(c(target_variable_delay, predictors_delay) / 3))
  
  # Start quarterfication loop according to Mariano and Murasawa #
  
  # Note: It is implicitly assumed that the stationary monthly data starts at the 
  #   second month of the first quarter. Further, it is assumed that the data set
  #   ends with an observation in the last month of the quarter.
  
  all_qrtly_data <- matrix(NaN, no_of_vars, no_of_qrtly_obs)
  for(t in seq(5, dim(monthly_predictors)[2], 3)){
    all_qrtly_data[1:no_of_target_vars, (t - 2)/3] <- target_variables[, t]
    all_qrtly_data[(no_of_target_vars + 1):(no_of_vars), (t - 2)/3] <- rowSums(
      cbind(1/3 * monthly_predictors[, t, drop = FALSE],
            2/3 * monthly_predictors[, t - 1, drop = FALSE],
            1 * monthly_predictors[, t - 2, drop = FALSE],
            2/3 * monthly_predictors[, t - 3, drop = FALSE],
            1/3 * monthly_predictors[, t - 4, drop = FALSE]),
      na.rm = TRUE)
  }
  
  return(zoo(t(all_qrtly_data), order.by = monthly_index[which(month(monthly_index) %% 3 == 0)]))
  
}
