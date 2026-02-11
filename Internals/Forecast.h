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

#ifndef FC
#define FC

// Including external libraries
#include <Eigen>
#include <stdlib.h>
#include <math.h>
#include <stdlib.h>


// Include internal libraries
#include "DataHandle.h"
#include "Filtering.h"

namespace Forecast {

    double factorForecaster(const Eigen::MatrixXd& X, const Filtering::KFS_fit& Fit, const int& K, const int& frequency, const int& delay, const int& h, const int& x0_ind);
    double factorDistance(const Eigen::MatrixXd& F_hat, const Eigen::MatrixXd& F);

};
#endif /* defined(FC) */

