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

#pragma once

#ifndef DATA_GENERATOR
#define DATA_GENERATOR

// Including external libraries
#include <cfloat>
#include <random>
#include <stdlib.h>
#include <Eigen>

// Including internal libraries
#include "DataHandle.h"

namespace DataGen {

    class FM // Class storing the results and parameters of the data generating process
    {
    public:
        Eigen::MatrixXd F; // Factor matrix
        Eigen::MatrixXd Phi; // Matrix governing the VAR processes of the factors themselves
        Eigen::MatrixXd Lambda; // Loading matrix
        Eigen::MatrixXd Sigma_e; // Correlation matrix of the error in the state equation
        Eigen::MatrixXd Sigma_epsilon; // Correlation matrix of the error term in the transition equation
        Eigen::MatrixXd e; // Idiosyncratic error term of the relevant variables
        Eigen::MatrixXd X; // Data matrix
        Eigen::VectorXi frequency; // Frequency vector
        Eigen::VectorXi date; // Date vector

        // Default Constructor
        FM()
        {
            F = Eigen::MatrixXd::Zero(0, 0);
            Phi = Eigen::MatrixXd::Zero(0, 0);
            Lambda = Eigen::MatrixXd::Zero(0, 0);
            Sigma_e = Eigen::MatrixXd::Zero(0, 0);
            Sigma_epsilon = Eigen::MatrixXd::Zero(0, 0);
            e = Eigen::MatrixXd::Zero(0, 0);
            X = Eigen::MatrixXd::Zero(0, 0);
            frequency = Eigen::VectorXi::Zero(0);
            date = Eigen::VectorXi::Zero(0);
        }

        // Constructor
        FM(int R, int T, int N)
        {
            F = Eigen::MatrixXd::Zero(R, T);
            Phi = Eigen::MatrixXd::Zero(R, R);
            Lambda = Eigen::MatrixXd::Zero(N, R);
            Sigma_e = Eigen::MatrixXd::Zero(N, N);
            Sigma_epsilon = Eigen::MatrixXd::Zero(R, R);
            e = Eigen::MatrixXd::Zero(N, T);
            X = Eigen::MatrixXd::Zero(T, N);
            frequency = Eigen::VectorXi::Zero(N);
            date = Eigen::VectorXi::Zero(N);
        }

        // Copy constructor
        FM(const FM& other) = default;

        // Move constructor
        FM(FM&& other) noexcept = default;

        // Copy assignment operator
        FM& operator=(const FM& other) = default;

        // Move assignment operator
        FM& operator=(FM&& other) noexcept = default;

        // Destructor
        ~FM() = default;
    };

    MatrixXd rndCorrMat(std::mt19937& gen, const double& beta_param = 1, const int& N = 100);
    MatrixXd MVData(std::mt19937& gen, const int& T, const int& N, const Eigen::VectorXd& m, const Eigen::MatrixXd& S, const bool& corr = false, const double& beta_param = 1.0);
    MatrixXd MVData(std::mt19937& gen, const int& T, const int& N, const bool& corr = false, const double& beta_param = 1.0);

    void staticFM(
        DataGen::FM& Results,
        int T,
        const int& N,
        const Eigen::MatrixXd& S,
        const Eigen::MatrixXd& Lambda,
        const Eigen::VectorXd& mu_e,
        const Eigen::MatrixXd& Sigma_e,
        const Eigen::MatrixXd& A,
        std::mt19937& gen,
        const int& order = 1,
        const bool& quarterfy = true,
        const bool& corr = true,
        const double& beta_param = 1.0,
        const double& m = 0.0,
        const int& R = 2,
        const int& burn_in = 1000,
        const bool& rescale = true
    );
};
#endif /* defined(DATA_GENERATOR) */
