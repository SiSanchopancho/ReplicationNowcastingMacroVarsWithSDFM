#define _USE_MATH_DEFINES // If you need some math constants

// Defin PI when using g++
#ifndef M_PI
#define M_PI (3.14159265358979323846)
#endif

#define NLOPT_DLL

// Externakl includes

#include <stdlib.h>
#include <cfloat>
#include <iostream>
#include <Eigen/Eigen>
#include <math.h>
#include <RcppCommon.h>
#include <Rcpp.h>
#include <RcppEigen.h>
#include "TwoStepSDFM_types.h"

// Internal Incldues

#include "Internals/CrossVal.h" // Cross-Validation Wrappers
#include "Internals/DataGen.h" // Data generation
#include "Internals/DataHandle.h" // Load and save data
#include "Internals/Developer.h" // DEvelopment helper
#include "Internals/ElNetSolve.h" // Elastic-Net solvers
#include "Internals/EmpiricsSetup.h" // Set-Up for handling an empiric data set
#include "Internals/Filtering.h" // Univariate representataion of multivariate KFS according to Koopman and Durbin (2010)
#include "Internals/Forecast.h" // Forecasting function
#include "Internals/Orders.h" // Functions to infer on "orders" of the proces (number of factors, VAR order, etc.)
#include "Internals/SparseDFM.h" // Wrapper for the Sparse DFM estimation
#include "Internals/SparsePCA.h" // Sparse Principal Components Analysis
#include "Internals/TargetedPredictors.h" // Targeted Predictors 

using namespace CrossVal;
using namespace DataGen;
using namespace DataHandle;
using namespace Developer;
using namespace ElNetSolve;
using namespace EmpiricsSetup;
using namespace Filtering;
using namespace Forecast;
using namespace Orders;
using namespace SparseDFM;
using namespace SparsePCA;
using namespace TargetedPredictors;
using namespace Rcpp;
using namespace RcppEigen;

/* Function to run the two-step SDFM estimation procedure*/
/*Source:
-  Franjic, Domenic and Schweikert, Karsten, Nowcasting Macroeconomic Variables with a Sparse Mixed Frequency Dynamic Factor Model (February 21, 2024). Available at SSRN: https://ssrn.com/abstract=4733872 or http://dx.doi.org/10.2139/ssrn.4733872 
*/

//' @description
//' This function is for internal use only and may change in future releases
//' without notice. Users should use `SimFM()` instead for a stable and
//' supported interface.
//'
// [[Rcpp::export]]
List runSDFMKFS(
	NumericMatrix X_in,
	IntegerVector delay,
	IntegerVector selected,
	int K,
	int order,
	bool decorr_errors,
	const char* crit,
	const char* method,
	double l2,
	NumericVector l1,
	double alpha,
	double l1_start,
	double ratio,
	int grid_size,
	int max_iterations,
	int steps,
	double comp_null,
	bool check_rank,
	double conv_crit,
	double conv_threshold,
	bool log,
	int KFS_conv_crit
)
{

	// Initialise the result object
	KFS_fit results;

	// Map the numeric matrices and vectors to eigen objects
	Eigen::Map<Eigen::MatrixXd> X_in_eigen(Rcpp::as<Eigen::Map<Eigen::MatrixXd>>(X_in));
	Eigen::Map<Eigen::VectorXi> delay_eigen(Rcpp::as<Eigen::Map<Eigen::VectorXi>>(delay));
	Eigen::Map<Eigen::VectorXi> selected_eigen(Rcpp::as<Eigen::Map<Eigen::VectorXi>>(selected));
	Eigen::Map<Eigen::VectorXd> l1_eigen(Rcpp::as<Eigen::Map<Eigen::VectorXd>>(l1));

	// Handle the case where l1, l1_start and or steps is not provided
	if (steps == -2147483647)
	{
		steps = INT_MIN;
	}

	if (selected == -2147483647)
	{
		steps = INT_MAX;
	}

	if (l1_start == -std::numeric_limits<double>::infinity())
	{
		l1_start = NAN;
	}

	// Estimate the sparse DFM
	SDFMKFS(results, X_in_eigen, delay_eigen, selected_eigen, K, order, decorr_errors, crit, method, l2,
		l1_eigen, alpha, l1_start, ratio, grid_size, max_iterations, steps, comp_null, check_rank,
		conv_crit, conv_threshold, log, KFS_conv_crit);

	// Re-correlate the loadings fit if necessary
	if (decorr_errors)
	{

		results.Lambda_hat = results.C.triangularView<Eigen::Lower>().solve(MatrixXd::Identity(X_in_eigen.cols(), X_in_eigen.cols())) * results.Lambda_hat;

		for (int n = 0; n < results.Lambda_hat.rows(); ++n) {
			for (int k = 0; k < results.Lambda_hat.cols(); ++k) {
				if (results.Zero_Indeces(n, k) == 1) {

					results.Lambda_hat(n, k) = 0.;

				}
			}
		}

	}

	// Convert the results back to Rcpp types and return
	return List::create(Named("Lambda_hat") = wrap(results.Lambda_hat),
		Named("Pt") = wrap(results.Pt),
		Named("F") = wrap(results.F),
		Named("Wt") = wrap(results.Wt),
		Named("C") = wrap(results.C),
		Named("P") = results.order);

}

/* Simulate an approximate DFM */

//' @description
//' This function is for internal use only and may change in future releases
//' without notice. Users should use `SimFM()` instead for a stable and
//' supported interface.
//'
// [[Rcpp::export]]
List runStaticFM(
	int T,
	const int& N,
	NumericMatrix S,
	NumericMatrix Lambda,
	NumericVector mu_e,
	NumericMatrix Sigma_e,
	NumericMatrix A,
	int order,
	bool quarterfy,
	bool corr,
	double beta_param,
	double m,
	int seed,
	int K,
	int burn_in,
	bool rescale
)
{

	// Map the numeric matrices and vectors to eigen objects

	Eigen::Map<Eigen::MatrixXd> S_eigen(Rcpp::as<Eigen::Map<Eigen::MatrixXd>>(S));
	Eigen::Map<Eigen::MatrixXd> Lambda_eigen(Rcpp::as<Eigen::Map<Eigen::MatrixXd>>(Lambda));
	Eigen::Map<Eigen::VectorXd> mu_e_eigen(Rcpp::as<Eigen::Map<Eigen::VectorXd>>(mu_e));
	Eigen::Map<Eigen::MatrixXd> Sigma_e_eigen(Rcpp::as<Eigen::Map<Eigen::MatrixXd>>(Sigma_e));
	Eigen::Map<Eigen::MatrixXd> A_eigen(Rcpp::as<Eigen::Map<Eigen::MatrixXd>>(A));

	FM results;
	StaticFM(results, T, N, S_eigen, Lambda_eigen, mu_e_eigen, Sigma_e_eigen, A_eigen, order, K, 1, quarterfy, corr, beta_param, m, seed, K, burn_in, rescale);

	return List::create(Named("F") = wrap(results.F),
		Named("Phi") = wrap(results.Phi),
		Named("Lambda") = wrap(results.Lambda),
		Named("Sigma_xi") = wrap(results.Sigma_e),
		Named("Sigma_epsilon") = wrap(results.Sigma_epsilon),
		Named("Xi") = wrap(results.e),
		Named("X") = wrap(results.X.transpose()),
		Named("frequency") = wrap(results.frequency));

}