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

#ifndef SPCA
#define SPCA

// Including external libraries
#include <Eigen>
#include <stdlib.h>
#include <math.h>

// Including internal libraries
#include "Developer.h"
#include "ElNetSolve.h"

namespace SparsePCA {

	class SPC_fit
	{
	public:
		Eigen::MatrixXd F_hat;
		Eigen::MatrixXd Lambda_hat;
		Eigen::VectorXd pct_var_expl;
		double tot_var;

        // Default Constructor
        SPC_fit()
        {
            F_hat = MatrixXd::Zero(0, 0);
            Lambda_hat = MatrixXd::Zero(0, 0);
            pct_var_expl = VectorXd::Zero(0);
            tot_var = 0.;
        }

        // Constructor
        SPC_fit(int K, int T) :
            tot_var(0)
        {
            F_hat = MatrixXd::Zero(K, T);
            Lambda_hat = MatrixXd::Zero(K, K);
            pct_var_expl = VectorXd::Zero(K);
        }

        // Copy constructor
        SPC_fit(const SPC_fit& other) = default;

        // Move constructor
        SPC_fit(SPC_fit&& other) noexcept = default;

        // Copy assignment operator
        SPC_fit& operator=(const SPC_fit& other) = default;

        // Move assignment operator
        SPC_fit& operator=(SPC_fit&& other) noexcept = default;

        // Destructor
        ~SPC_fit() = default;
	};

    extern
        void SparsePC(
            SPC_fit& results,
            const Eigen::MatrixXd& X_in,
            const Eigen::VectorXi& selected,
            const int& K,
            const double& l2,
            Eigen::VectorXd l1,
            const double& alpha,
            double l1_start,
            double ratio,
            const int& grid_size,
            const int& max_iterations,
            int steps,
            const double& comp_null,
            const bool& check_rank,
            const double& conv_crit,
            const double& conv_threshold,
            const bool& normalise = 0
            );

};
#endif /* defined(SPCA) */
