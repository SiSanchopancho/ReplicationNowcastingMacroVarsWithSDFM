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

#pragma once

#ifndef CV
#define CV

// Including third party libraries
#include <Eigen>
#include <stdlib.h>
#include <math.h>
#include <random>
#include <chrono>
#include <omp.h>

// Including internal libraries
#include "SparseDFM.h"
#include "Forecast.h"

namespace CrossVal {

    /* Class storing the results of the cross-validation procedure */
    class CV_res {
    public:
        Eigen::VectorXi sel; // Number of variables aimed to be extracted for each factor
        Eigen::VectorXd l1; // Factor specific l1 penalty
        double l2; // Global l2 penalty
        double min_mse; // Minimum nowcasting error (MSE)

        // Default constructor
        CV_res()
        {
            sel = Eigen::VectorXi::Zero(0);
            l1 = Eigen::VectorXd::Zero(0);
            l2 = 0.;
            min_mse = 0.;
        }

        // Constructor
        CV_res(int R)
        {
            sel = Eigen::VectorXi::Zero(0);
            l1 = Eigen::VectorXd::Zero(0);
            l2 = 0.;
            min_mse = 0.;

        }

        // Copy constructor
        CV_res(const CV_res& other) = default;

        // Move constructor
        CV_res(CV_res&& other) = default;

        // Copy assignment operator
        CV_res& operator=(const CV_res& other) = default;

        // Move assignment operator
        CV_res& operator=(CV_res&& other) = default;

        // Destructor
        ~CV_res() = default;
    };

    extern

        CrossVal::CV_res randCrossVal(
            const int& x0_ind,
            const Eigen::MatrixXd& X,
            const int& H_in,
            const int& R,
            const int& max_models,
            const Eigen::VectorXi& delay,
            const Eigen::VectorXi& date,
            const Eigen::VectorXi& frequency,
            std::mt19937& gen,
            const bool& decorr_errors = false,
            const int& order = 10,
            const double& l2_min = -6,
            const double& l2_max = +6,
            const int& sel_min = 5,
            int sel_max = 10e+6,
            const char* crit = "BIC",
            const double& comp_null = 10e-15,
            const double& conv_threshold = 10e-6,
            const bool& timer = true,
            const bool& log = false
        );

    CrossVal::CV_res parallelRandCrossVal(
        const int& x0_ind,
        const Eigen::MatrixXd& X,
        const int& H_in,
        const int& R,
        const int& max_models,
        const Eigen::VectorXi& delay,
        const Eigen::VectorXi& date,
        const Eigen::VectorXi& frequency,
        std::mt19937& gen,
        const bool& decorr_errors = false,
        const int& order = 10,
        const double& l2_min = -6,
        const double& l2_max = +6,
        const int& sel_min = 5,
        int sel_max = 10e+6,
        const char* crit = "BIC",
        const double& comp_null = 10e-15,
        const double& conv_threshold = 10e-6,
        const bool& timer = true,
        const bool& log = false
    );

    CrossVal::CV_res randBIC(
        const Eigen::MatrixXd& X_in,
        const int& R,
        const int& max_models,
        const Eigen::VectorXi& delay,
        const Eigen::VectorXi& date,
        const Eigen::VectorXi& frequency,
        std::mt19937& gen,
        const bool& decorr_errors = false,
        const int& order = 10,
        const double& l2_min = -6,
        const double& l2_max = +6,
        const int& sel_min = 5,
        int sel_max = 10e+6,
        const char* crit = "BIC",
        const double& comp_null = 10e-15,
        const double& conv_threshold = 10e-6,
        const bool& timer = true,
        const bool& log = false
    );

    CrossVal::CV_res parallelRandBIC(
        const Eigen::MatrixXd& X_in,
        const int& R,
        const int& max_models,
        const Eigen::VectorXi& delay,
        const Eigen::VectorXi& date,
        const Eigen::VectorXi& frequency,
        std::mt19937& gen,
        const bool& decorr_errors = false,
        const int& order = 10,
        const double& l2_min = -6,
        const double& l2_max = +6,
        const int& sel_min = 5,
        int sel_max = 10e+6,
        const char* crit = "BIC",
        const double& comp_null = 10e-15,
        const double& conv_threshold = 10e-6,
        const bool& timer = true,
        const bool& log = false
    );

    CrossVal::CV_res parallelRandSinglePenaltyBIC(
        const Eigen::MatrixXd& X_in,
        const int& R,
        const int& max_models,
        const Eigen::VectorXi& delay,
        const Eigen::VectorXi& date,
        const Eigen::VectorXi& frequency,
        std::mt19937& gen,
        const bool& decorr_errors = false,
        const int& order = 10,
        const double& l1_max = NAN,
        const double& l1_min = NAN,
        const double& l2_min = -6,
        const double& l2_max = +6,
        const char* crit = "BIC",
        const double& comp_null = 10e-15,
        const double& conv_threshold = 10e-6,
        const bool& timer = true,
        const bool& log = false
    );

};
#endif /* defined(CV) */

