/* SPDX-License-Identifier: GPL-3.0-or-later */
/*
 * Copyright © 2024-2026 Domenic Franjic
 *
 * This file is part of ReplicationNowcastingMacroVarsWithSDFM.
 *
 * ReplicationNowcastingMacroVarsWithSDFM is free software: you can redistribute
 * it and/or modify it under the terms of the GNU General Public License as
 * published by the Free Software Foundation, either version 3 of the License,
 * or (at your option) any later version.

 * ReplicationNowcastingMacroVarsWithSDFM is distributed in the hope that it
 * will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
 * of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.

 * You should have received a copy of the GNU General Public License
 * along with ReplicationNowcastingMacroVarsWithSDFM. If not, see <https://www.gnu.org/licenses/>.
 */

#include "Filtering.h"

/* Univariate representation of the Multivariate Kalman Filter */
/* Sources:
Koopman, S. J., & Durbin, J. (2000). Fast filtering and smoothing for multivariate state space models. Journal of time series analysis, 21(3), 281-296. (https://doi.org/10.1111/1467-9892.00186)
*/
void Filtering::UVMVKalmanFilter(
    Filtering::KFS_fit& results, // Object for storing the results
    Eigen::MatrixXd& Vt_inv, // Inverse of the variance matrix V
    Eigen::MatrixXd& Kt, // Matrix storing the Kalman gain at time t
    Eigen::MatrixXd& e, //Matrix storing the errors
    const int& N, // Number of variables
    const int& T, // Number of observations
    const int& K, // Number of factors
    const Eigen::MatrixXd& X_in, // Data matrix of dimensions (NxT) !
    const Eigen::MatrixXd& Lambda, // Factor Loading Matrix
    const Eigen::MatrixXd& Phi, // VAR coefficient matrix for factor process
    const Eigen::MatrixXd& Sigma_e, // Idiosyncratic error variance-covariance-matrix
    const Eigen::MatrixXd& Sigma_epsilon, // Variance-covariance matrix of the factor VAR process
    const Eigen::MatrixXd& Pt0, // Initial conditional variance-covariance matrix of the factors 
    const Eigen::VectorXd& Ft0, // Initial value of the factors
    const Eigen::VectorXi& missings, // Vector indicating the obseravtion index for missing observations for each variable x_n in X
    const int& h // Forecasting horizon
)
{

    /* Dummies */

    // Reals
    double Vt = 0.;

    // Vectors
    Eigen::VectorXd Ftt = Eigen::VectorXd::Zero(K);

    // Matrices
    Eigen::MatrixXd Ptt = Eigen::MatrixXd::Zero(K, K), Lambda_T = Lambda.transpose(), X = X_in;

    // Initialisation
    results.F.col(0) = Ft0;
    results.Pt.block(0, 0, K, K) = Pt0;

    /* Filter loop */
    for (int t = 0; t < T; ++t)
    {
        Ftt = results.F.col(t);
        Ptt = results.Pt.block(0, (t) * K, K, K);
        for (int n = 0; n < N; ++n)
        {
            if (T - missings(n) <= t)
            {
                // Skip missing observations
                continue;
            }

            // Update
            Vt = Lambda.row(n) * Ptt * Lambda_T.col(n) + Sigma_e.diagonal()(n);
            e(n, t) = X(n, t) - Lambda.row(n) * Ftt;
            Vt_inv(n, t) = 1. / Vt;
            Kt.col(t * N + n) = Ptt * Lambda_T.col(n);
            Ftt += (Vt_inv(n, t)) * Kt.col(t * N + n) * e(n, t);
            Ptt -= Kt.col(t * N + n) * (Vt_inv(n, t)) * Kt.col(t * N + n).transpose();

        }

        // Store
        results.F.col(t) = Ftt;
        results.Pt.block(0, (t)*K, K, K) = Ptt;

        // Propagate
        results.F.col(t + 1) = Phi * Ftt;
        results.Pt.block(0, (t + 1) * K, K, K) = Phi * Ptt * Phi.transpose() + Sigma_epsilon;
    }

    /* End filter loop */

    // Forecasting
    if (0 < h)
    {
        for (int t = T; t < T + h; ++t)
        {
            results.F.col(t + 1) = Phi * results.F.col(t);
            results.Pt.block(0, (t + 1) * K, K, K) = Phi * results.Pt.block(0, t * K, K, K) * Phi.transpose() + Sigma_epsilon;
        }
    }
    return;
}

/* Univariate representation of the Multivariate Kalman Smoother */
/* Sources:
Koopman, S. J., & Durbin, J. (2000). Fast filtering and smoothing for multivariate state space models. Journal of time series analysis, 21(3), 281-296. (https://doi.org/10.1111/1467-9892.00186)
*/
void Filtering::UVMVKalmanSmoother(
    Filtering::KFS_fit& results, // Object for storing the results
    Eigen::MatrixXd& Kt, // Matrix storing the Kalman gain at time t
    Eigen::MatrixXd& Vt_inv, // Inverse of the variance matrix V
    Eigen::MatrixXd& e, //Matrix storing the errors
    const int& N, // Number of variables
    const int& T, // Number of observations
    const int& K, // Number of factors
    const Eigen::MatrixXd& X, // Data matrix of dimensions (NxT) !
    const Eigen::MatrixXd& Lambda, // Factor Loading Matrix
    const Eigen::MatrixXd& Phi, // VAR coefficient matrix for factor process
    const Eigen::MatrixXd& Sigma_e, // Idiosyncratic error variance-covariance-matrix
    const Eigen::MatrixXd& Sigma_epsilon, // Variance-covariance matrix of the factor VAR process
    const Eigen::VectorXi& missings // Vector indicating the obseravtion index for missing observations for each variable x_n in X
)
{

    /* Dummies */
     
    // Vectors
    Eigen::VectorXd rtt = Eigen::VectorXd::Zero(K);

    // Matrices
    Eigen::MatrixXd Lt = Eigen::MatrixXd::Zero(K, K), Ntt = Eigen::MatrixXd::Zero(K, K), Lambda_T = Lambda.transpose(), IdentK = Eigen::MatrixXd::Identity(K, K);


    /* Smoother loop */
    for (int t = T - 1; 0 <= t; --t)
    {

        for (int n = N - 1; 0 <= n; --n)
        {
            if (T - missings(n) <= t)
            {
                // Skip missing observations
                continue;
            }

            // Update
            Lt = IdentK - Kt.col(t * N + n) * Lambda.row(n) * Vt_inv(n, t);
            rtt = Lambda_T.col(n) * Vt_inv(n, t) * e(n, t) + Lt.transpose() * rtt;
            Ntt = Lambda_T.col(n) * Vt_inv(n, t) * Lambda.row(n) + Lt.transpose() * Ntt * Lt;
        
        }
        
        // Smooth
        results.F.col(t) += results.Pt.block(0, t * K, K, K) * rtt;
        results.Wt.block(0, t * K, K, K) = results.Pt.block(0, t * K, K, K) - results.Pt.block(0, t * K, K, K) * Ntt * results.Pt.block(0, t * K, K, K);
        
        // Propagate
        rtt = Phi.transpose() * rtt;
        Ntt = Phi.transpose() * Ntt * Phi;
    
    }

    /* End smoother loop */

    return;
}

