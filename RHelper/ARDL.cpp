/* SPDX-License-Identifier: GPL-3.0-or-later */
/*
 * Copyright (C) 2024 Domenic Franjic
 *
 * This file is part of TwoStepSDFM.
 *
 * TwoStepSDFM is free software: you can redistribute
 * it and/or modify it under the terms of the GNU General Public License as
 * published by the Free Software Foundation, either version 3 of the License,
 * or (at your option) any later version.

 * TwoStepSDFM is distributed in the hope that it
 * will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
 * of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.

 * You should have received a copy of the GNU General Public License
 * along with TwoStepSDFM. If not, see <https://www.gnu.org/licenses/>.
 */

 // [[Rcpp::depends(RcppEigen)]]
#include "ARDL.h"

 //' @description
 //' This function is for internal use only and may change in future releases
 //' without notice. 
 //'
 // [[Rcpp::export]]
Rcpp::List runARDL(
	Rcpp::NumericVector target_variable,
	Rcpp::NumericVector target_variable_predictor,
	Rcpp::NumericVector predictor_variable,
	const unsigned max_target_lags,
	const unsigned max_predictor_lags,
	const std::string crit
)
{

	// Disable Eigen parallelisation as this function is regularly employed in aparallelised environment
	Eigen::setNbThreads(1);

	Eigen::Map<Eigen::VectorXd> target_variable_eigen(Rcpp::as<Eigen::Map<Eigen::VectorXd>>(target_variable));
	Eigen::Map<Eigen::VectorXd> target_variable_predictor_eigen(Rcpp::as<Eigen::Map<Eigen::VectorXd>>(target_variable_predictor));
	Eigen::Map<Eigen::VectorXd> predictor_variable_eigen(Rcpp::as<Eigen::Map<Eigen::VectorXd>>(predictor_variable));
	const unsigned no_of_target_observations = target_variable_predictor.size();
	const unsigned no_of_variables = max_target_lags + 1 // max_target_lags number of lags plus the "contemporaneous" observation
		+ max_predictor_lags + 1 // max_predictor_lags number of lags plus the "contemporaneous" observation
		+ 1; // Add one for the intercept
	const unsigned max_lag = std::max<unsigned>(max_target_lags, max_predictor_lags);
	const unsigned effective_no_of_obs = no_of_target_observations - max_lag;

	// Build full lag predictor matrix
	Eigen::MatrixXd predictor_matrix = Eigen::MatrixXd::Constant(effective_no_of_obs, no_of_variables, 1.0);
	Eigen::Vector<Eigen::Index, Eigen::Dynamic> target_lag_sequence = Eigen::Vector<Eigen::Index, Eigen::Dynamic>::LinSpaced(max_target_lags + 1, max_target_lags, 0);
	Eigen::Vector<Eigen::Index, Eigen::Dynamic> predictor_lag_sequence = Eigen::Vector<Eigen::Index, Eigen::Dynamic>::LinSpaced(max_predictor_lags + 1, max_predictor_lags, 0);
	for (unsigned lag = 0; lag <= max_lag; ++lag) {
		if (lag <= max_target_lags) {
			predictor_matrix.col(lag + 1) = target_variable_predictor_eigen.segment(target_lag_sequence(lag), effective_no_of_obs);
		}
		if (lag <= max_predictor_lags) {
			predictor_matrix.col(max_target_lags + 2 + lag) = predictor_variable_eigen.segment(predictor_lag_sequence(lag), effective_no_of_obs);
		}
	}

	// Fit ARDL model
	Eigen::LLT<Eigen::MatrixXd> llt;
	double residuals_sum_of_squares = 0.0;
	double current_ic = 0.0;
	double best_ic = DBL_MAX;
	Eigen::VectorXd beta_optimal = Eigen::VectorXd::Zero(0);
	Eigen::VectorXd current_beta = Eigen::VectorXd::Zero(0);
	Eigen::Vector2i optimal_lags = Eigen::Vector2i::Zero(2);
	Eigen::MatrixXd current_predictor_matrix = Eigen::MatrixXd::Zero(0, 0);
	for (unsigned target_lag = 0; target_lag <= max_target_lags; ++target_lag) {
		for (unsigned predictor_lag = 0; predictor_lag <= max_predictor_lags; ++predictor_lag) {
			current_predictor_matrix = Eigen::MatrixXd::Zero(effective_no_of_obs, target_lag + 1 + predictor_lag + 1 + 1);
			current_predictor_matrix.leftCols(target_lag + 2) = predictor_matrix.leftCols(target_lag + 2);
			current_predictor_matrix.rightCols(predictor_lag + 1) = predictor_matrix.block(0, max_target_lags + 2, effective_no_of_obs, predictor_lag + 1);

			llt.compute(current_predictor_matrix.transpose() * current_predictor_matrix);

			if (llt.info() != Eigen::Success) {
				Rcpp::Rcerr << "\nWARNING: LLT Decomposition FAILED in runARDL.\n";
				current_beta = Eigen::VectorXd::Zero(current_predictor_matrix.cols());
				residuals_sum_of_squares = DBL_MAX;
				continue;
			}

			current_beta = llt.solve(current_predictor_matrix.transpose() * target_variable_eigen.tail(effective_no_of_obs));

			// Compute RSS and evaluat AIC/BIC
			residuals_sum_of_squares = (target_variable_eigen.tail(effective_no_of_obs) - current_predictor_matrix * current_beta).squaredNorm();
			double log_like_equivalent = static_cast<double>(effective_no_of_obs) * std::log(residuals_sum_of_squares / static_cast<double>(effective_no_of_obs));
			if (crit == "AIC") {
				current_ic = log_like_equivalent + 2 * static_cast<double>(target_lag + predictor_lag + 2);
				if (current_ic < best_ic) {
					best_ic = current_ic;
					beta_optimal = current_beta;
					optimal_lags(0) = target_lag;
					optimal_lags(1) = predictor_lag;
				}
			}
			else if (crit == "BIC") {
				current_ic = log_like_equivalent + std::log(static_cast<double>(effective_no_of_obs)) * static_cast<double>(target_lag + predictor_lag + 2);
				if (current_ic < best_ic) {
					best_ic = current_ic;
					beta_optimal = current_beta;
					optimal_lags(0) = target_lag;
					optimal_lags(1) = predictor_lag;
				}
			}
		}
	}

	// Enable Eigen parallelisation for proper clean up
	Eigen::setNbThreads(0);

	// Convert the results back to Rcpp types and return
	return Rcpp::List::create(Rcpp::Named("coefficients") = Rcpp::wrap(beta_optimal),
		Rcpp::Named("optimL_lag_order") = Rcpp::wrap(optimal_lags)
	);

}

//' @description
//' This function is for internal use only and may change in future releases
//' without notice. 
//'
// [[Rcpp::export]]
Rcpp::List runDL(
	Rcpp::NumericVector target_variable,
	Rcpp::NumericVector predictor_variable,
	const unsigned max_predictor_lags,
	const std::string crit
)
{

	// Disable Eigen parallelisation as this function is regularly employed in aparallelised environment
	Eigen::setNbThreads(1);

	Eigen::Map<Eigen::VectorXd> target_variable_eigen(Rcpp::as<Eigen::Map<Eigen::VectorXd>>(target_variable));
	Eigen::Map<Eigen::VectorXd> predictor_variable_eigen(Rcpp::as<Eigen::Map<Eigen::VectorXd>>(predictor_variable));
	const unsigned no_of_target_observations = target_variable.size();
	const unsigned no_of_variables = max_predictor_lags + 1 + 1; // Add one for the intercept and one for the contemporaneous effect
	const unsigned effective_no_of_obs = no_of_target_observations - max_predictor_lags;

	// Build full lag predictor matrix
	Eigen::MatrixXd predictor_matrix = Eigen::MatrixXd::Constant(effective_no_of_obs, no_of_variables, 1.0);
	Eigen::Vector<Eigen::Index, Eigen::Dynamic> predictor_lag_sequence = Eigen::Vector<Eigen::Index, Eigen::Dynamic>::LinSpaced(max_predictor_lags + 1, max_predictor_lags, 0);
	for (unsigned lag = 0; lag <= max_predictor_lags; ++lag) {
		predictor_matrix.col(1 + lag) = predictor_variable_eigen.segment(predictor_lag_sequence(lag), effective_no_of_obs);
	}

	// Fit distributed lag model
	Eigen::LLT<Eigen::MatrixXd> llt;
	double residuals_sum_of_squares = 0.0;
	double current_ic = 0.0;
	double best_ic = DBL_MAX;
	Eigen::VectorXd beta_optimal = Eigen::VectorXd::Zero(0);
	Eigen::VectorXd current_beta = Eigen::VectorXd::Zero(0);
	int optimal_lag = 0;
	Eigen::MatrixXd current_predictor_matrix = Eigen::MatrixXd::Zero(0, 0);
	for (unsigned predictor_lag = 0; predictor_lag <= max_predictor_lags; ++predictor_lag) {
		current_predictor_matrix = predictor_matrix.leftCols(predictor_lag + 2);

		llt.compute(current_predictor_matrix.transpose() * current_predictor_matrix);

		if (llt.info() != Eigen::Success) {
			Rcpp::Rcerr << "\nWARNING: LLT Decomposition FAILED in runARDL.\n";
			current_beta = Eigen::VectorXd::Zero(current_predictor_matrix.cols());
			residuals_sum_of_squares = DBL_MAX;
			continue;
		}

		current_beta = llt.solve(current_predictor_matrix.transpose() * target_variable_eigen.tail(effective_no_of_obs));

		// Compute RSS and evaluat AIC/BIC
		residuals_sum_of_squares = (target_variable_eigen.tail(effective_no_of_obs) - current_predictor_matrix * current_beta).squaredNorm();
		double log_like_equivalent = static_cast<double>(effective_no_of_obs) * std::log(residuals_sum_of_squares / static_cast<double>(effective_no_of_obs));
		if (crit == "AIC") {
			current_ic =  log_like_equivalent + 2 * static_cast<double>(predictor_lag + 2);
			if (current_ic < best_ic) {
				best_ic = current_ic;
				beta_optimal = current_beta;
				optimal_lag = predictor_lag;
			}
		}
		else if (crit == "BIC") {
			current_ic = log_like_equivalent + std::log(static_cast<double>(effective_no_of_obs)) * static_cast<double>(predictor_lag + 2);
			if (current_ic < best_ic) {
				best_ic = current_ic;
				beta_optimal = current_beta;
				optimal_lag = predictor_lag;
			}
		}
	}

	// Enable Eigen parallelisation for proper clean up
	Eigen::setNbThreads(0);

	// Convert the results back to Rcpp types and return
	return Rcpp::List::create(Rcpp::Named("coefficients") = Rcpp::wrap(beta_optimal),
		Rcpp::Named("optimL_lag_order") = Rcpp::wrap(optimal_lag)
	);

}