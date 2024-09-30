/* SPDX-License-Identifier: GPL-3.0-or-later */
/*
 * Copyright © 2024 Domenic Franjic
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

#include "SparsePCA.h"

/* Sparse Principal Components Analysis */
/* Source:
- Zou, H., Hastie, T., & Tibshirani, R. (2006). Sparse principal component analysis. Journal of computational and graphical statistics, 15(2), 265-286. (https://doi.org/10.1198/106186006X113430)
- Zou, H., Hastie, T., & Zou, M. H. (2016). Package ‘elasticnet’. (https://cran.r-project.org/web/packages/elasticnet/index.html)
*/
void SparsePCA::SparsePC(
    SPC_fit& results, // Object for saving the results
    const Eigen::MatrixXd& X_in, // Data matrix
    const Eigen::VectorXi& selected, // Vector with number of variables aimed to be selected for each factor
    const int& K, // Number of factors
    const double& l2, // l2 penalty
    Eigen::VectorXd l1, // Vector of l1 penalties for each factor
    const double& alpha, // Weighting coefficient between LASSO and RIDGE (Only "CD"; currently disabled)
    double l1_start, // Starting value of the l1 penalty sequence (Only "CD"; currently disabled)
    double ratio, // Ratio betwen the smalles and largest value of the l1-pebnalty sequence (Only "CD"; currently disabled)
    const int& grid_size, // Size of the penalty grid (Only "CD"; currently disabled)
    const int& max_iterations, // Maximum number of iterations (Only "CD"; currently disabled)
    int steps, // Maximum nzumber number of steps ("LARS" only)
    const double& comp_null, // Computational zero
    const bool& check_rank, // Check rank
    const double& conv_crit, // Conversion criterion for the Loadings matrix
    const double& conv_threshold, // Conversion criterion (Only "CD"; currently disabled)
    const bool& normalise // indicate whether PCs should be normalised
)
{
 
    /* Singular value decomposition of the data */
    
    // Demean the data
    Eigen::MatrixXd X = X_in;

    /* Dummies */ 

    // Integers
    const int N = X.cols(), T = X.rows();

    // Reals
    double tot_var = 0.;

    // Decompositions
    Eigen::JacobiSVD<Eigen::MatrixXd> SVD_i, X_SVD(X, Eigen::ComputeThinU | Eigen::ComputeThinV);
    Eigen::ColPivHouseholderQR<Eigen::MatrixXd> QR;

    // Vectors
    Eigen::VectorXd y_j = Eigen::VectorXd::Zero(T), var_exp = Eigen::VectorXd::Zero(K);

    // Matrices
    Eigen::MatrixXd R = Eigen::MatrixXd::Zero(T, K), F_hat = Eigen::MatrixXd::Zero(T, K), GB = Eigen::MatrixXd::Zero(N, N), Gram = Eigen::MatrixXd::Zero(N, N), A = Eigen::MatrixXd::Zero(N, K), B = Eigen::MatrixXd::Zero(N, K), U = X_SVD.matrixU(),
        V = X_SVD.matrixV(), Sigma = X_SVD.singularValues().asDiagonal(), Lambda_hat = Eigen::MatrixXd::Zero(N, K);

    // Save total variance for later
    tot_var = (X_SVD.singularValues().array() * X_SVD.singularValues().array()).sum();

    // Calculate Gram matrix
    Gram = X.transpose() * X;

    // Check rank of the matrix
    if (check_rank)
    {
        FullPivLU<Eigen::MatrixXd> Gram_LU(Gram);
        if (Gram_LU.rank() < N)
        {
            std::cout << '\n' << "Error! The Matrix product XT * X does not have full rank." << '\n';
            return;
        }
        else
        {
            std::cout << '\n' << "The Matrix product XT * X has full rank." << '\n';
        }
    }

    /* Peuso EM loop */

    // Initialise A = V(,0:k-1)
    A = V(all, seq(0, K - 1));

    // Solve ^b_j = arg min_{b_j} || X * a_j - X * b_j||^2 + l||b_2||^2 + l_{1,j}||b_j||_1 for j = 1,...,k
    for (int k = 0; k < K; ++k)
    {
        // Create target vector
        y_j = X * A.col(k);

        B.col(k) = LARS<false>(y_j, X, l2, l1(k), selected(k), steps, comp_null);
    }

    // Update the current loadings matrix
    Lambda_hat = B;

    /* Rince and repeat 1 and 2 until stopping rule or convergence */

    for (int i = 0; i < max_iterations; ++i)
    {
        // Update A
        GB = Gram * B;
        SVD_i.compute(GB, Eigen::ComputeThinU | Eigen::ComputeThinV);
        A = SVD_i.matrixU() * SVD_i.matrixV().transpose();

        for (int k = 0; k < K; ++k)
        {
            // create target vector
            y_j = X * A.col(k);
            B.col(k) = LARS<false>(y_j, X, l2, l1(k), selected(k), steps, comp_null);
        }
        // Check for convergence of the Loadings matrix
        if ((Lambda_hat - B).squaredNorm() <= conv_crit)
        {
            Lambda_hat = B;
            break;
        }
        else
        {
            Lambda_hat = B;
        }
    }

    if (normalise)
    {
        // 4 Step: Normalise V
        for (int k = 0; k < K; ++k)
        {
            Lambda_hat.col(k).normalize();
        }
    }

    // Calculate the factors
    F_hat = X * Lambda_hat;
    F_hat.transposeInPlace();

    // Calculate variance explained
    R = QR.compute(F_hat.transpose()).matrixQR().triangularView<Eigen::Upper>();
    var_exp = 1 / tot_var * (R.diagonal().array() * R.diagonal().array()).matrix();

    // Save the results as an SPC_fit object
    results.F_hat = F_hat;
    results.Lambda_hat = Lambda_hat;
    results.pct_var_expl = var_exp;
    results.tot_var = tot_var;

    return;
}