/* SPDX-License-Identifier: GPL-3.0-or-later */
/*
 * Copyright © 2026 Domenic Franjic
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

#include "DataGen.h"

/* Draw random correlation matrices */
Eigen::MatrixXd DataGen::rndCorrMat(std::mt19937& gen, const double& beta_param, const int& N)
{
    /* Dummies */

    // Matrices
    Eigen::MatrixXd P = Eigen::MatrixXd::Zero(N, N); // Permutation matrix
    Eigen::MatrixXd S = Eigen::MatrixXd::Identity(N, N);

    // Distributions
    std::gamma_distribution<> dist1(beta_param, 1.0);
    std::gamma_distribution<> dist2(beta_param, 1.0);

    // Reals
    double pki, pli, plk;

    // Construct off diagonal entries of S using P
    for (int k = 0; k < N - 1; ++k)
    {
        for (int i = k + 1; i < N; ++i)
        {
            double gamma1 = dist1(gen);
            double gamma2 = dist2(gen);
            P(k, i) = 2 * (gamma1 / (gamma1 + gamma2)) - 1;
            pki = P(k, i);
            if (k != 0)
            {
                int l = (k - 1);
                while (l >= 0)
                {
                    pli = P(l, i);
                    plk = P(l, k);
                    pki = pki * sqrt((1 - pow(pli, 2)) * (1 - pow(plk, 2))) + pli * plk;
                    --l;
                }
            }
            S(k, i) = pki;
            S(i, k) = pki;
        }
    }

    // Permute S
    Eigen::PermutationMatrix<Dynamic, Dynamic> Pi(N);
    Pi.setIdentity();
    std::shuffle(Pi.indices().data(), Pi.indices().data() + Pi.indices().size(), gen);

    return S(Pi.indices(), Pi.indices());
}

/* Draw Standard Multivariate Normal Data */
Eigen::MatrixXd DataGen::MVData(std::mt19937& gen, const int& T, const int& N, const bool& corr, const double& beta_param)
{

    // draw Matrix with standard normal entries

    // Dummies

    // Matrices

    Eigen::MatrixXd X = Eigen::MatrixXd::Zero(T, N);

    // Distribution
    std::normal_distribution<double> dist(0., 1.);


    for (int t = 0; t < T; ++t)
    {
        for (int n = 0; n < N; ++n)
        {
            X(t, n) = dist(gen);
        }
    }

    // If data should be correlated: corralte data using L * X' where L stems from R = L * L' and R is the correlation matrix
    if (corr == 1)
    {
        Eigen::MatrixXd R = rndCorrMat(gen, beta_param, N);
        Eigen::LLT<Eigen::MatrixXd> Rllt(R);
        Eigen::MatrixXd U = Rllt.matrixU();
        X = X * U;
    }

    return X;
}

/* Draw Multivariate Normal Data*/
Eigen::MatrixXd DataGen::MVData(std::mt19937& gen, const int& T, const int& N, const VectorXd& m, const Eigen::MatrixXd& S, const bool& corr, const double& beta_param)
{
    /* Dummies */

    // Matrices
    Eigen::MatrixXd X = Eigen::MatrixXd::Zero(T, N);

    // Distributions
    std::normal_distribution<> dist(0., 1.);

    for (int t = 0; t < T; ++t)
    {
        for (int n = 0; n < N; ++n)
        {
            X(t, n) = dist(gen);
        }
    }

    if (corr == 1)
    {
        // If data should be correlated: corralte data using L * X' where L stems from R = L * L' and R is the correlation matrix
        Eigen::MatrixXd R = rndCorrMat(gen, beta_param, N);
        Eigen::LLT<Eigen::MatrixXd> Rllt(R);
        Eigen::MatrixXd U = Rllt.matrixU();
        X = X * U;

    }
    if (!S.isDiagonal())
    {

        // If correlation matrix is proviede externally via S: Left-multiply data by the upper triangular of the LLT decomposition of S
        LLT<Eigen::MatrixXd> Rllt(S);
        Eigen::MatrixXd U = Rllt.matrixU();
        X = X * U;

    }
    else
    {
        // If data should have non unit variance: Multiply data by diag(s1, ..., sn) from the right
        X = X * S.cwiseSqrt();

    }

    // If data should have non-zero mean: Add mean to each column of the data
    for (int n = 0; n < N; ++n)
    {
        for (int t = 0; t < T; ++t)
        {
            X(t, n) += m(n);
        }
    }

    return X;
}


/* Generate a FM with potentially correlated factors following a VAR(order) process using the factor process VAR companion form */
void DataGen::staticFM(
    DataGen::FM& Results, // FM class object storing the results
    int T, // Time Dimension
    const int& N, // Number of variables
    const Eigen::MatrixXd& S, // Diagonal 2*R x x2*R variance matrix of the transition errors
    const Eigen::MatrixXd& Lambda, // Factor loadings Matrix
    const VectorXd& mu_e, // Mean of factor model innovations 
    const Eigen::MatrixXd& Sigma_e, // Variance of factor model innovations
    const Eigen::MatrixXd& A, // Correlation matrix of the static factors
    std::mt19937& gen, // RNG
    const int& order, // Order of the factor VAR proces
    const bool& quarterfy, // Indicate whether to convert the first m * N observations into quarterly data
    const bool& corr, // Indicator whether measurement error of FM should be cross correlated
    const double& beta_param, // Degree of correlation of idiosyncratic error (Elements of the parameters are beta distributed where the parameter is in (0, 1])
    const double& m, // Fraction of quarterly observations (only if quarterfy == 1)
    const int& R, // Number of factors
    const int& burn_in, // burning in period
    const bool& rescale // internaly rescaling the error variance
)
{

    /* Adjust the number of observations */
    
    // Add additional observations to weaken the influence of the starting values
    T += burn_in;

    // If observations will be quarterfied, two additional observations are needed to end up at T observations
    if (quarterfy)
    {
        T += 2;
    }

    /* Dummies */

    // Integers
    int factors = R * order;

    // Matrices
    Eigen::MatrixXd U = Eigen::MatrixXd::Zero(factors, T + 1), F = Eigen::MatrixXd::Zero(factors, T + 1), 
      Phi = Eigen::MatrixXd::Zero(factors, factors), L = Eigen::MatrixXd::Zero(N, factors), X = Eigen::MatrixXd::Zero(T, N);
    
    /* Create factors */

    // Fill in placeholders
    U.topLeftCorner(R, T + 1) = MVData(gen, T + 1, R, VectorXd::Zero(R), S).transpose();
    Results.Sigma_epsilon = S;
    F.col(0) = MVData(gen, factors, 1);
    Phi.topLeftCorner(R, factors) = A;

    if (1 < order)
    {
        Phi.bottomLeftCorner(factors - R, factors - R) = Eigen::MatrixXd::Identity(factors - R, factors - R);
    }

    for (int t = 1; t < T + 1; ++t)
    {
        F.col(t) = Phi * F.col(t - 1) + U.col(t);
    }

    // Remove the first observation which
    DataHandle::removeCol(F, 0, 0);

    /* Create the data */

    L.topLeftCorner(N, R) = Lambda;

    // Create random vector for innovations and rescale if wished
    Eigen::MatrixXd VarCov;

    if (corr)
    {
        Eigen::MatrixXd R = rndCorrMat(gen, beta_param, N);
        if (!R.isApprox(R.transpose()))
        {
            std::cout << '\n' << "Error! Random Covariance-Matrix is not symmetric." << '\n';
            return;
        }
        VarCov = Sigma_e.cwiseSqrt() * R * Sigma_e.cwiseSqrt();
    }
    else
    {
        VarCov = Sigma_e;
    }
    if (rescale)
    {
        Eigen::MatrixXd Common = L * F;
        Eigen::MatrixXd Scaler = Eigen::MatrixXd::Identity(N, N);
        Scaler.diagonal() = DataHandle::cov(Common.transpose()).diagonal().cwiseSqrt();
        for (int n = 0; n < N; ++n)
        {
            if (Scaler.diagonal()(n) == 0)
            {
                Scaler.diagonal()(n) = 1;
            }
        }
        VarCov = Scaler * VarCov * Scaler;
    }

    Eigen::MatrixXd e = MVData(gen, T, N, mu_e, VarCov, false, beta_param);
    e.transposeInPlace();
    Results.Sigma_e = VarCov;

    // Create observations
    for (int t = 0; t < T; ++t)
    {
        X(t, seq(0, N - 1)) = L * F.col(t) + e.col(t);
    }

    /* Quarterfication */

    if (quarterfy)
    {

        // Quarterfy the first floor(m*N) observations according to Mariano-Murasawa
        for (int q = 0; q < std::floor(m * double(N)); ++q)
        {
            VectorXd xq = VectorXd::Zero(T);
            int t = 4;
            while (t < T)
            {
                xq(t) = 1. / 3. * (X(t, q) + 2. * X(t - 1, q) + 3. * X(t - 2, q) + 2. * X(t - 3, q) + X(t - 4, q));
                t += 3;
            }

            t = 4;
            while (t < T)
            {
                for (int tt = t - 2; tt <= t; ++tt)
                {
                    X(tt, q) = xq(t);
                }
                t += 3;
            }
        }

        // Remove the two additional observations that where used for aggregation
        DataHandle::removeRow(X, 0, 0);
        DataHandle::removeRow(X, 0, 0);
        DataHandle::removeCol(F, 0, 0);
        DataHandle::removeCol(F, 0, 0);
        DataHandle::removeCol(e, 0, 0);
        DataHandle::removeCol(e, 0, 0);
    }

    // Binding the results

    Results.F = F(seq(0, R - 1), seq(burn_in - 1, last));
    Results.Lambda = Lambda;
    Results.Phi = A;
    Results.X = X(seq(burn_in - 1, last), all);
    Results.frequency = VectorXi::Constant(N, 12);
    Results.e = e(all, seq(burn_in - 1, last));

    if (quarterfy)
    {
        Results.frequency(seq(0, std::floor(m * double(N)) - 1)).setConstant(4);
    }

    // Pseudo-date vector which might be necessary later
    Results.date = VectorXi::LinSpaced(X.rows(), 0, X.rows() - 1);

    return;
}
