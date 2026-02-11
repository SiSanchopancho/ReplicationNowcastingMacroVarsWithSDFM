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

#ifndef ORDERS
#define ORDERS

#define _USE_MATH_DEFINES

// Including external libraries
#include <Eigen>
#include <stdlib.h>
#include <math.h>

#include <iostream>

/* Infer the order of a VAR process using IC criteria */
enum Criterion {
    AIC,
    BIC,
    HIC
};
inline void computeVAROrder(
    Eigen::MatrixXd& ICs_res, // Matrix to store the results for each order and criterion
    const Eigen::MatrixXd& F, // Data matrix
    const int& O, // Maximum order to be tested
    const double& comp_null // Comptational zero
)
{

    /* Dummies */

    // Integers
    int K = static_cast<int>(F.rows()), T = static_cast<int>(F.cols());

    /* Loop over orders */
    for (int o = 1; o <= O; ++o)
    {

        // Companion form
        Eigen::MatrixXd F_curr(K * o, T);
        for (int oo = 0; oo < o; ++oo)
        {
            F_curr.block(K * oo, oo, K, T - oo) = F.block(0, 0, K, T - oo);
        }

        // OLS
        Eigen::MatrixXd F_t = F_curr(Eigen::all, Eigen::seq(o + 1, Eigen::last)).transpose();
        Eigen::MatrixXd F_t_lag = F_curr(Eigen::all, Eigen::seq(o, Eigen::last - 1)).transpose();
        Eigen::MatrixXd Phi = ((F_t_lag.transpose() * F_t_lag).llt().solve(Eigen::MatrixXd::Identity(K * o, K * o)) * F_t_lag.transpose() * F_t).transpose();
        Phi = Phi.unaryExpr([comp_null](double x) {return (comp_null < std::abs(x)) ? x : 0.; }).eval();
        Eigen::MatrixXd F_hat = (Phi.topLeftCorner(K, K * o) * F_curr(Eigen::all, Eigen::seq(o, Eigen::last - 1))).transpose();

        // Calculate the residuals variance-covariance-matrix and its determinant
        Eigen::MatrixXd Res_Var_Cov = (F_t(0, Eigen::seq(0, K - 1)) - F_hat.row(0)).transpose() * (F_t(0, Eigen::seq(0, K - 1)) - F_hat.row(0));
        for (int t = 1; t < T - 1 - o; ++t)
        {
            Res_Var_Cov += (F_t(t, Eigen::seq(0, K - 1)) - F_hat.row(t)).transpose() * (F_t(t, Eigen::seq(0, K - 1)) - F_hat.row(t));
        }
        Res_Var_Cov /= double(T);
        Eigen::MatrixXd C = Res_Var_Cov.llt().matrixL();
        double det = 2 * C.diagonal().array().log().sum();

        // Calculate the information criterion
        ICs_res(o - 1, 0) = o;
        ICs_res(o - 1, 1) = det + 2. * double(K * K * o) / double(T); // AIC
        ICs_res(o - 1, 2) = det + double(K * K * o) * std::log(double(T)) / double(T); // BIC
        ICs_res(o - 1, 3) = det + 2. * double(K * K * o) * std::log(std::log(double(T))) / double(T); // HIC
    }

    /* End of order loop */

    return;
}
template<Criterion criterion>
int VARorder(const Eigen::MatrixXd& F, const int& O, const double& comp_null);
template<>
inline int VARorder<AIC>(const Eigen::MatrixXd& F, const int& O, const double& comp_null) {
    Eigen::MatrixXd ICs_res = Eigen::MatrixXd::Constant(O, 4, DBL_MAX);
    computeVAROrder(ICs_res, F, O, comp_null);
    int order = 0;
    ICs_res.col(1).minCoeff(&order);
    return order + 1;

}
template<>
inline int VARorder<BIC>(const Eigen::MatrixXd& F, const int& O, const double& comp_null) {
    Eigen::MatrixXd ICs_res = Eigen::MatrixXd::Constant(O, 4, DBL_MAX);
    computeVAROrder(ICs_res, F, O, comp_null);
    int order = 0;
    ICs_res.col(2).minCoeff(&order);
    return order + 1;
}
template<>
inline int VARorder<HIC>(const Eigen::MatrixXd& F, const int& O, const double& comp_null) {
    Eigen::MatrixXd ICs_res = Eigen::MatrixXd::Constant(O, 4, DBL_MAX);
    computeVAROrder(ICs_res, F, O, comp_null);
    int order = 0;
    ICs_res.col(3).minCoeff(&order);
    return order + 1;
}

#endif /* defined(ORDERS) */

