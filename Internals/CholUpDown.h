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

#pragma once

#ifndef CHOL_UP_DOWN
#define CHOL_UP_DOWN

// Including external libraries
#include <Eigen>
#include <math.h>
#include <cfloat>

// Including internal libraries
#include "DataHandle.h"

namespace CholUpDown {
    void cholUpdate(Eigen::MatrixXd& L, const Eigen::VectorXd& c, const Eigen::MatrixXd& X, const int& index, const double& l2);
    void cholDowndate(Eigen::MatrixXd& L, const Eigen::VectorXd& c, const Eigen::MatrixXd& X, const int& index, const double& l2);
}


#endif /* defined(CHOL_UP_DOWN) */

